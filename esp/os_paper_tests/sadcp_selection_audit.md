# S-ADCP selection audit: SOLO Figure 6 and LATTE Figure 12

Audited 30 October 2026 against `/home/z5297792/ESP_zonodo`, commit `a5ce3a7bb9fc21560c3963fc5024db81d6a24184`. Results below come from executing the current notebook selection cells and tracing the actual fitting functions with the local data. No source code, notebook, or data was modified. Cell numbers are **one-based positions including Markdown cells**, not execution counts.

## Time basis, source data, and units

Both notebooks, cell 3, read `data/sadcp_data.mat`: **2,174 time records × 79 depths**, with depths 27–807 in steps of 10. The source spans **2023-10-09 11:56:08 to 2023-11-01 04:37:59**, an elapsed **22 d 16 h 41 min 50 s**. Both first retain `Lat <= -37` (degrees north): **433 records** spanning **2023-10-17 10:27:31 to 2023-10-21 22:27:30**, elapsed **4 d 11 h 59 min 59 s**. These spans describe the first and last timestamps, not continuous coverage.

All displayed timestamps and durations are rounded independently to the nearest second; fractional endpoints are supplied below for reproducibility. Times are reported on the requested **UTC convention**. Satellite UTC is explicit in its NetCDF metadata. **The MATLAB file has no time-zone or unit metadata**: its S-ADCP times are converted exactly as the notebooks do, from Julian days using `2000-01-01 12:00 + (tim_sadcp - 2451545) days`, yielding timezone-naive timestamps. Their UTC provenance cannot be independently verified from this file. Similarly, S-ADCP depth in metres and velocity in m/s follow the notebook interpretation; the MATLAB arrays have no unit attributes. Satellite velocities explicitly have units m/s.

Both notebooks use Cartesian coordinates in metres:

```python
x = (Lon - 154.3) * 111000 * cos(deg2rad(-37.9))
y = (Lat + 37.9) * 111000
```

## SOLO — `SOLO_transect_data.ipynb`

Cell 3 specifies `t0 = 2023-10-18 01:30:00`, `dt = 0.5 days`, and the inclusive selection `t0-dt <= time <= t0+dt`. Thus **the nominal window is 24 hours**, **2023-10-17 13:30:00–2023-10-18 13:30:00**, rather than a 12-hour total window. Nearest-depth selection `.sel(z=37, method="nearest")` selects **37 m** exactly.

The **96 retained records** actually span **2023-10-17 13:42:31–2023-10-18 13:27:30**, elapsed **23 h 45 min 00 s**. Unrounded endpoints are `2023-10-17 13:42:30.666659` and `2023-10-18 13:27:30.333340` (elapsed `23:44:59.666681`). Their xarray mean timestamp is `2023-10-18 01:35:00.562501248`, or **01:35:01 UTC** to the nearest second.

| Stage | Count | Actual selection |
|---|---:|---|
| Source | 2,174 | All S-ADCP timestamps, all 79 depths |
| Latitude-filtered | 433 | `Lat <= -37` |
| Time/depth-filtered | 96 | Inclusive nominal window; z = 37 m |
| Valid original vectors | 96 | Finite `x, y, u, v` |
| Initial along-transect normal-velocity fit | 96 | Finite projected `l, VT, VN` |
| Inner-core refit | **30** | `abs(l - x0_initial) <= 30000` m |
| Outer spatial candidates | 93 | `rho² <= 100000²` m² about the fitted centre |
| Outer fit after sign and validity filters | **93** | All 93 candidates pass |

**Inner-core sample meaning:** cell 7 calls `project_sadcp_to_transect` (`functions.py:109–148`), which drops nonfinite vectors, fits a straight line, sorts the projected positions, and **interpolates velocities onto 96 evenly spaced transect positions**. These are resampled points, not a subset of 96 original timestamps. `solo` (`functions.py:23–69`) first fits the normal-velocity cubic to all 96 points, obtaining `x0_initial = 90200.693249767 m`; it then selects 30 points within ±30 km of that initial zero crossing. The retained `l` range is `62105.092994063–118387.833519933 m`. It refits the normal velocity and fits the tangential velocity on those same 30 points, yielding `x0_refined = 84534.743152263 m`. It does not reselect after this refit.

**Outer-core sample meaning:** cell 7 uses the original 96 vectors and `Q = I`, with centre `(xc, yc) = (-9434.639930002, 2507.191383929) m`. Consequently `rho² = (x-xc)² + (y-yc)²`, and 93 original vectors satisfy the 100 km cutoff. `out_core_param_fit` (`functions.py:478–501`) additionally retains `vt <= 0` for this negative initial angular velocity, finite `rho², Qr, vt`, `rho² >= 0`, and `Qr != 0`, where `Qr = |Q(x-xc, y-yc)|`. All 93 pass. Velocities enter in m/s. **The outer fit includes inner-region vectors; it is not an annulus.**

Cell 13's magenta map overlay uses the different rule `hypot(x-xc, y-yc) <= 30000`, selecting **20 original vectors**. That display count must not be reported as the 30-point inner-core fit count.

## LATTE — `LATTE_app1_multiplatform_data.ipynb`

Cell 3 retains the full **433-record** latitude-filtered set stated above, without an additional S-ADCP time window. Its unrounded endpoints are `2023-10-17 10:27:31.000004` and `2023-10-21 22:27:30.333340` (elapsed `4 days 11:59:59.333336`). `latte_source_selector` selects its default `z_target=37.0` by nearest depth, obtaining **37 m**, and all **433 vectors have finite x, y, u, v** (`functions.py:691–697`).

Cell 5 reads `data/satellite_data.nc`, a **41 × 46** grid (1,886 locations), covering latitude −42° to −34° and longitude 149° to 158° at 0.2° spacing. It replaces each velocity component failing `abs(component) < 4.9 m/s` with NaN; it interpolates onto daily **06:00 UTC** fields, retaining **2023-10-19 06:00:00–2024-01-27 06:00:00** (101 fields).

The selector then linearly interpolates these satellite fields to the **mean timestamp of all 433 latitude-filtered S-ADCP records**, calculated **before** depth/finite-vector screening: `2023-10-19 16:27:30.509622172`, or **2023-10-19 16:27:31 UTC** to the nearest second (`functions.py:680–689`). The bracketing daily fields are 19 and 20 October at 06:00 UTC. At that time, **1,747 satellite vectors** have finite x, y, u, v. Thus the combined input has **433 S-ADCP + 1,747 satellite = 2,180 vectors**.

### Combined (`source="multi"`) selection counts

| Stage | S-ADCP | Satellite | Total |
|---|---:|---:|---:|
| Finite combined input | 433 | 1,747 | **2,180** |
| First inner-core selection | 173 | 12 | **185** |
| Refined inner-core selection / second LATTE fit | 174 | 12 | **186** |
| Outer spatial selection, including inner points | 395 | 80 | **475** |
| Outer fit after rotation-sign filter | 393 | 80 | **473** |
| Outer fit after remaining validity filters | **393** | **80** | **473** |
| Display rows labelled `core="outer"` | 221 | 68 | **289** |

Cell 7 calls the selector without radius overrides. The precise masks (`functions.py:711–771`) are:

1. **First core:** `x² + y² <= 35000²` m², using the default prior centre `(0,0)` and `Q=I`. The first `latte` fit uses all 185 vectors.
2. **Refined core:** evaluate `rho² = q11*dx² + 2*q12*dx*dy + q22*dy²` over **the entire combined input**, about the first fitted centre `(-345.156650844, -422.599277080) m`, with `Q1 = [[1.137989611957, -0.044924660013], [-0.044924660013, 0.880516144039]]`; retain `rho² <= 35000²` m². This is a fresh selection, not an intersection with the first mask. The second LATTE fit uses 186 vectors; its own finite-vector check drops none.
3. **Outer spatial mask:** recompute `rho²` about the refined centre `(-372.710282652, -429.926758603) m`, using `Q = [[1.144555154095, -0.046655243357], [-0.046655243357, 0.875603685980]]`; retain `rho² <= 100000²` m². These 475 vectors include the inner-region data.
4. **Actual outer fit:** `out_core_param_fit` retains `vt <= 0` for the negative fitted angular velocity, removing **two S-ADCP vectors**; its finite/nonnegative-radius/nonzero-`Qr` checks remove no more. The fit therefore uses **473**, not 475 or 433, vectors. Here `rho` and `Qr` are in metres and `vt` in m/s; the fitted normalised velocity is `vt * rho / Qr`.

The selector's returned `df_xyuv` labels refined-core rows as `inner` and `outer & ~core2` rows as `outer`. These **289 outer-labelled plotting rows are not the outer-fit count**. The returned `source="multi"` label identifies the combined experiment, not each row's physical platform; counts above track the S-ADCP-first/satellite-second concatenation. The stored `land_mask` is not an additional selection mask in the selector; finite satellite velocities determine eligibility. Plot limits and quiver spacing also do not filter the fit.

For completeness, cell 7 also runs the existing single-platform comparisons:

| Run | Finite input (S-ADCP / satellite / total) | First core | Refined core | Outer spatial | Actual outer fit |
|---|---|---|---|---|---|
| `sadcp` | 433 / 0 / 433 | 173 / 0 / 173 | 174 / 0 / 174 | 395 / 0 / 395 | 393 / 0 / 393 |
| `sat` | 0 / 1,747 / 1,747 | 0 / 12 / 12 | 0 / 9 / 9 | 0 / 82 / 82 | 0 / 82 / 82 |

### Radius mismatch

**The implemented LATTE defaults are 35 km / 100 km**, explicitly `rho_core=35000.0`, `rho_outer=100000.0` (`functions.py:661–663`). Cell 7 supplies neither override. Therefore the **30 km inner-core manuscript claim supplied in the audit brief conflicts with this implementation**; the 100 km outer limit agrees. `out_core_param_fit` uses a separate 30 km split for colouring its plotted points (`functions.py:561–564`), which does not change the selector's 35 km inner-core fits. The manuscript itself was not supplied and its wording has not been independently checked. No 30 km alternative experiment was run.

## Saved notebook state versus current execution

SOLO's saved cell 3 count (96), fitted centre, `Rc = 67805.071560 m`, mean timestamp, and cell 14 `R² = 0.948525` agree with execution to displayed precision. The absolute satellite path in SOLO cell 5 exists; its file is byte-identical to the bundled `data/satellite_data.nc`. This satellite field is not an input to SOLO's cell 7 fit. Cell 5 also overwrites `t0` with the satellite dataset's first time, so its final nearest-time satellite selection is not controlled by the transect midpoint variable originally defined in cell 3.

LATTE's saved count (433), interpolation time, centres, Q matrices, axis ratios, and outer-fit parameters agree with current execution to displayed precision. **Its saved R² values do not:**

| Run | Saved cells 8 / 12 R² | Current cell 8 R² |
|---|---:|---:|
| Combined | 0.907048 | 0.907558923 |
| S-ADCP | 0.902745 | 0.899771753 |
| Satellite | 0.969165 | 0.967839083 |

The cause of that saved-output difference is **not verifiable from the saved notebook alone**. Saved execution counts are also nonsequential (plot cells 13/14 precede a results cell saved as execution 9). Saved figures are not evidence of exact selection counts; all counts here use current execution.

## Reproduction

Executed with `/home/z5297792/myenv310/bin/python -B`: Python 3.10.8, NumPy 2.2.6, pandas 2.3.3, SciPy 1.15.3, xarray 2025.6.1, netCDF4 1.7.3. The following self-contained script executes the original selection cells without saving notebooks and captures the fitting masks directly. Save it under `/tmp` and run with that interpreter and `-B` to avoid bytecode writes. Matplotlib output/cache stays in `/tmp`. It prints per-platform counts in **S-ADCP, satellite, total** order.

```python
import os
os.environ.update(MPLBACKEND="Agg", MPLCONFIGDIR="/tmp/sadcp-audit-mpl")
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
root = Path("/home/z5297792/ESP_zonodo")
os.chdir(root)
sys.path.insert(0, str(root))
import functions
seen = []
def trace(frame, event, arg):
    if (event == "return" and frame.f_code.co_filename == str(root / "functions.py")
        and frame.f_code.co_name in {"solo", "fit_x0_from_vn",
                                    "latte_source_selector", "out_core_param_fit"}):
        seen.append((frame.f_code.co_name, frame.f_locals.copy(), arg))
    return trace
for name in ["SOLO_transect_data", "LATTE_app1_multiplatform_data"]:
    nb = json.loads((root / (name + ".ipynb")).read_text())
    ns = {}
    def cell(n):
        exec(compile("".join(nb["cells"][n-1]["source"]),
                     f"{name}:cell{n}", "exec"), ns)
    seen.clear()
    for n in [1, 3, 5]:
        cell(n)
    t = ns["ds_sadcp"].time
    print(name, "records", t.size, "endpoints", t.min().values, t.max().values,
          "elapsed", pd.Timedelta(t.max().values - t.min().values),
          "mean", t.mean().values)
    sys.settrace(trace)
    try:
        cell(7)
    finally:
        sys.settrace(None)
    outer_calls = [d for n, d, r in seen if n == "out_core_param_fit"]
    if name.startswith("SOLO"):
        print("depth", ns["ds_sadcp"].z.values,
              "valid projected points", len(ns["df_proj"]))
        for n, d, r in seen:
            if n == "fit_x0_from_vn":
                print("VN fit points / x0", len(d["l"]), r[0])
        o = outer_calls[0]
        print("outer spatial / sign / final", ns["outer"].sum(),
              o["sign_mask"].sum(), o["m"].sum())
        for n in [9, 11]:
            cell(n)
    else:
        selectors = [d for n, d, r in seen if n == "latte_source_selector"]
        for d, o in zip(selectors, outer_calls):
            labels = (np.r_[np.zeros(len(d["xd"]), int),
                            np.ones(len(d["xg"]), int)] if d["source"] == "multi"
                      else np.full(len(d["xi"]), int(d["source"] == "sat")))
            def count(indices):
                a = labels[indices]
                return [int((a == 0).sum()), int((a == 1).sum()), len(a)]
            print(d["source"], "depth/time", d["z0"], d["time"])
            for key in ["core1", "core2", "outer", "outer_only"]:
                print(key, count(d[key]))
            idx = np.flatnonzero(d["outer"])[o["sign_mask"]]
            print("finite input", count(np.arange(len(labels))),
                  "outer sign", count(idx), "outer final", count(idx[o["m"]]))
        cell(8)
    print(ns["df_res"].to_string(index=False))
    plt.close("all")
```

SHA-256 of the actual inputs/source used (relative to the repository):

```text
functions.py                         10208a788eea38e9451edc3a2f31aaa351f043e3cb6ed5a07c52b28b5bc05655
SOLO_transect_data.ipynb             b8b99957d5aacf60c223642c0414d46e26f513270f17a08a78e516770016f14f
LATTE_app1_multiplatform_data.ipynb   63078df9c019cd9a4f68fdc5f2de13843d2647f4f9dbf45a6a54a343d9896752
data/sadcp_data.mat                  9fd84d62f636ae0f5e90488a48af717b7764333c8e5fcb9605a918539f3620a8
data/satellite_data.nc               1c0f9625f946b1ccfc1ea53895b1d21ed7e57dbbed48570e49cfa28742adf053
```

## Manuscript-ready factual sentences

SOLO retained 96 S-ADCP records at 37 m depth, used 30 resampled transect points for inner-core estimation, and fitted 93 original vectors within 100 km of the estimated centre.

LATTE used 174 S-ADCP and 12 satellite vectors for its refined 35 km elliptical inner-core fit, and 393 S-ADCP and 80 satellite vectors for its sign-filtered fit within an elliptical radius of 100 km.
