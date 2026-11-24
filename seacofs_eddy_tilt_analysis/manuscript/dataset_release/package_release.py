"""Assemble and SHA256/CRC-verify the local review ZIP; never upload or edit data.

Run from manuscript/: python dataset_release/package_release.py
Run prepare_metadata.py first if rebuilding the static grid or source inventory.
"""
from pathlib import Path, PurePosixPath
from datetime import datetime, timezone
import hashlib
import io
import json
import os
import zipfile

HERE=Path(__file__).resolve().parent
REVIEW=HERE/'outputs/review'
OUTPUT=HERE/'outputs/seacofs_eddy_dataset_v1.0.0.zip'


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):result.update(chunk)
    return result.hexdigest()


def run():
    manifest_path=REVIEW/'release_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    source=Path(manifest['source_dataset'])
    expected={row['path']:row['sha256'] for row in manifest['source_files']}
    actual={str(p.relative_to(source)) for p in source.rglob('*') if p.is_file()}
    if actual != set(expected):
        raise ValueError('Source inventory changed since audit; review before packaging')
    files={name:source/name for name in expected}
    for name in manifest['prepared_files']:
        p=REVIEW/name if name=='static_grid.nc' else HERE/name
        if sha(p)!=manifest['prepared_files'][name]:
            raise ValueError(f'{name} changed since preparation; refresh reviewed metadata first')
        files[name]=p;expected[name]=manifest['prepared_files'][name]
    files['release_manifest.json']=manifest_path
    expected['release_manifest.json']=sha(manifest_path)
    files['checksums.sha256']=REVIEW/'checksums.sha256'
    listed={name:digest for digest,name in (line.split('  ',1) for line in files['checksums.sha256'].read_text().splitlines())}
    if listed!=expected:raise ValueError('Checksum list is inconsistent with the release manifest')
    expected['checksums.sha256']=sha(files['checksums.sha256'])
    for name in files:
        if PurePosixPath(name).is_absolute() or '..' in PurePosixPath(name).parts:
            raise ValueError(f'Unsafe archive name: {name}')
    temporary=OUTPUT.with_suffix('.zip.partial')
    try:
        with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=1,allowZip64=True) as z:
            for index,(name,p) in enumerate(files.items(),1):
                state=p.stat();digest=hashlib.sha256()
                with p.open('rb') as src,z.open(name,'w',force_zip64=True) as dst:
                    for chunk in iter(lambda:src.read(8*1024*1024),b''):
                        digest.update(chunk);dst.write(chunk)
                after=p.stat()
                if digest.hexdigest()!=expected[name] or (state.st_size,state.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):
                    raise ValueError(f'Source changed or checksum differs: {name}')
                if index%150==0:print(f'Packaged {index}/{len(files)} files',flush=True)
        print('Checking archived SHA256 values and ZIP CRCs...',flush=True)
        with zipfile.ZipFile(temporary) as z:
            if len(z.namelist())!=len(files) or set(z.namelist())!=set(files):
                raise ValueError('Archive members differ from expected inventory')
            for name in z.namelist():
                digest=hashlib.sha256()
                with z.open(name) as f:
                    for chunk in iter(lambda:f.read(8*1024*1024),b''):digest.update(chunk)
                if digest.hexdigest()!=expected[name]:raise ValueError(f'Archive verification failed: {name}')
            # Exercise the public file formats, beyond byte/CRC integrity.
            import numpy as np
            import netCDF4 as nc
            import pyarrow.parquet as pq
            for name in ['processed/eddy_dataset_processed.parquet',
                         'tilt/tilt_dataset.parquet',
                         'analysis/tilt_summary.parquet']:
                pq.read_table(io.BytesIO(z.read(name)))
            depth_array=np.load(io.BytesIO(z.read('z_r.npy')),allow_pickle=False)
            with nc.Dataset('release_grid',memory=z.read('static_grid.nc')) as grid:
                assert depth_array.shape[:2] == (len(grid.dimensions['x']),len(grid.dimensions['y']))
                assert len(grid.dimensions['profile_level'])==depth_array.shape[2]
        os.replace(temporary,OUTPUT)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    digest=sha(OUTPUT)
    OUTPUT.with_suffix('.zip.sha256').write_text(f'{digest}  {OUTPUT.name}\n')
    # Provide separately downloadable documentation for the eventual Zenodo page.
    for name in ['README.md','variables.csv','generation_config.yaml','CITATION.cff','LICENSE.txt']:
        (REVIEW/name).write_bytes((HERE/name).read_bytes())
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),archive=OUTPUT.name,
        sha256=digest,bytes=OUTPUT.stat().st_size,archive_members=len(files),
        original_source_files=len(manifest['source_files']),
        all_member_sha256_and_crc_verified=True,
        primary_parquet_numpy_netcdf_readability_verified=True,
        dataset_doi=manifest.get('dataset_doi'),
        publication_status='Prepared locally; dataset DOI reserved; not uploaded/published by this script')
    (REVIEW/'package_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':run()
