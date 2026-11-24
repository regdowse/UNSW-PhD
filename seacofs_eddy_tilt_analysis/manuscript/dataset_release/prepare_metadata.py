"""Build local release metadata/static grid; never modify or copy source tables.

Run with the SEACOFS Python environment. Outputs are in outputs/review/.
README.md and variables.csv are maintained beside this script for review.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import zipfile
import yaml
import numpy as np
import netCDF4 as nc
import pyarrow.parquet as pq

HERE=Path(__file__).resolve().parent
DATA=Path('/srv/scratch/z5297792/SEACOFS_26yr_eddy_dataset_modular')
MODEL=Path('/srv/scratch/z3533156/26year_BRAN2020/outer_avg_01461.nc')
ARCHIVE=Path('/tmp/eddy-dataset-pipeline-v0.1.1.zip')
OUT=HERE/'outputs/review'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def distance(lat1,lon1,lat2,lon2):
    lat1,lon1,lat2,lon2=map(np.radians,[lat1,lon1,lat2,lon2])
    a=np.sin((lat2-lat1)/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2
    return 6357.*2*np.arctan2(np.sqrt(a),np.sqrt(1-a))


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    # Use the retained audit copy of the exact archived YAML, without altering paths.
    config=HERE.parent/'context/revision/outputs/dataset_release_audit/generation_config_as_archived.yaml'
    if ARCHIVE.exists():
        z=zipfile.ZipFile(ARCHIVE);base=z.namelist()[0].split('/')[0]+'/'
        assert config.read_bytes()==z.read(base+'config/local.yaml')
    (HERE/'generation_config.yaml').write_bytes(config.read_bytes())
    with nc.Dataset(MODEL) as f:
        arrays={name:np.asarray(f[name][:]).T.copy() for name in ['lon_rho','lat_rho','h','mask_rho','angle']}
        reference=float(f['angle'][0,0])
    lon=arrays['lon_rho'];lat=arrays['lat_rho'];nx,ny=lon.shape
    x=np.r_[0,np.cumsum(distance(lat[:-1,ny//2],lon[:-1,ny//2],lat[1:,ny//2],lon[1:,ny//2]))]
    y=np.r_[0,np.cumsum(distance(lat[nx//2,:-1],lon[nx//2,:-1],lat[nx//2,1:],lon[nx//2,1:]))]
    zr=np.load(DATA/'z_r.npy',mmap_mode='r')
    assert zr.shape[:2]==lon.shape
    depths=np.r_[0,np.abs(zr[150,150,1:])]
    grid=OUT/'static_grid.nc'
    definitions={
        'x':('km','Cumulative model-x coordinate; ic indexes this dimension'),
        'y':('km','Cumulative model-y coordinate; jc indexes this dimension'),
        'lon_rho':('degrees_east','Longitude at original model rho points, transposed'),
        'lat_rho':('degrees_north','Latitude at original model rho points, transposed'),
        'h':('m','Model bathymetric depth, positive downward'),
        'mask_rho':('1','Model rho-point mask: 1 water, 0 land'),
        'angle':('radian','Original model-grid angle, transposed'),
        'pipeline_reference_angle':('radian','Single angle at original [0,0] used by extraction'),
        'profile_depth':('m','Surface plus positive-down reference profile depths'),
        'profile_level':('1','Profile z index: zero surface; subsequent fixed target levels')}
    with nc.Dataset(grid,'w',format='NETCDF4') as f:
        f.createDimension('x',nx);f.createDimension('y',ny);f.createDimension('profile_level',len(depths))
        f.title='Static SEACOFS geometry for the eddy dataset'
        f.source_model_file=MODEL.name
        f.source_model_doi='10.5281/zenodo.12703063'
        f.pipeline_doi='10.5281/zenodo.24457378'
        f.array_order='x,y; transposed from original model y,x'
        f.coordinate_construction='Middle-row/column cumulative spherical distances; Earth radius 6357 km; not a standard map projection'
        f.z_r_array_order='Separate z_r.npy uses x,y,native_level, uppermost level first; signed negative depths in metres'
        f.history='Release-time extraction of existing static arrays; no model/eddy analysis rerun'
        values={**arrays,'x':x,'y':y,'profile_level':np.arange(len(depths)),'profile_depth':depths,'pipeline_reference_angle':reference}
        for name,a in values.items():
            dims=() if name=='pipeline_reference_angle' else ('profile_level',) if name in ['profile_level','profile_depth'] else (name,) if name in ['x','y'] else ('x','y')
            dtype='i4' if name=='profile_level' else np.asarray(a).dtype
            v=f.createVariable(name,dtype,dims,zlib=bool(dims),complevel=4)
            v[:]=a;v.units,v.long_name=definitions[name]
            if name in ['profile_depth','h']:v.positive='down'
            if name=='mask_rho':v.flag_values=np.array([0,1],dtype=np.asarray(a).dtype);v.flag_meanings='land water'
    with nc.Dataset(grid) as f:
        for name,a in arrays.items():np.testing.assert_array_equal(f[name][:],a)
        np.testing.assert_array_equal(f['profile_depth'][:],depths)
        np.testing.assert_array_equal(f['x'][:],x);np.testing.assert_array_equal(f['y'][:],y)
        assert float(f['pipeline_reference_angle'][:])==reference
    source_rows=[]
    for p in sorted(DATA.rglob('*')):
        if not p.is_file():continue
        row=dict(path=str(p.relative_to(DATA)),bytes=p.stat().st_size,sha256=sha(p))
        if p.suffix=='.parquet':row['rows']=pq.ParquetFile(p).metadata.num_rows
        source_rows.append(row)
    audit=HERE.parent/'context/revision/outputs/dataset_release_audit/file_inventory.csv'
    import csv
    with audit.open(newline='') as f:
        expected={r['path']:r['sha256'] for r in csv.DictReader(f)}
    if expected != {r['path']:r['sha256'] for r in source_rows}:
        raise ValueError('Source snapshot changed since audit; review it before preparing a release')
    citation=yaml.safe_load((HERE/'CITATION.cff').read_text())
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='prepared locally; dataset DOI reserved, publication not verified',
        dataset_doi=citation.get('doi'),
        dataset_version=str(citation['version']),creators=citation['authors'],license=citation['license'],
        pipeline_doi='10.5281/zenodo.24457378',pipeline_commit='143b9a31ec0586206cf440df3da129fb886bc4b8',
        generation_esp_version='v1.1',compatible_esp_version='v1.2 (byte-identical computational functions.py)',
        source_dataset=str(DATA),source_files=source_rows,
        static_grid=dict(file=grid.name,source_model_file=str(MODEL),shape=[nx,ny],reference_angle_rad=reference,sha256=sha(grid)),
        z_r_sha256=sha(DATA/'z_r.npy'),configuration_sha256=sha(HERE/'generation_config.yaml'),
        prepared_files={p.name:sha(p) for p in [HERE/'README.md',HERE/'variables.csv',HERE/'generation_config.yaml',HERE/'CITATION.cff',HERE/'LICENSE.txt',grid]},
        preparation_code_sha256=sha(__file__))
    (OUT/'release_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with (OUT/'checksums.sha256').open('w') as f:
        for row in source_rows:f.write(f'{row["sha256"]}  {row["path"]}\n')
        for p in [HERE/'README.md',HERE/'variables.csv',HERE/'generation_config.yaml',HERE/'CITATION.cff',HERE/'LICENSE.txt',grid]:f.write(f'{sha(p)}  {p.name}\n')
        f.write(f'{sha(OUT/"release_manifest.json")}  release_manifest.json\n')
    print('Prepared',grid,'bytes',grid.stat().st_size,flush=True)
    print('Source snapshot:',len(source_rows),'files;',sum(r['bytes'] for r in source_rows),'bytes',flush=True)
    print('Documentation checksums and release manifest:',OUT,flush=True)


if __name__=='__main__':run()
