# Task: Audit S-ADCP sample selection for manuscript reviewer response

Repository: `regdowse/elliptical_streamfunction_parametrisation` whcih can be found locally `/home/z5297792/ESP_zonodo`

Inspect `SOLO_transect_data.ipynb`, `LATTE_app1_multiplatform_data.ipynb`, `functions.py`, and the referenced data files. **Do not modify code, notebooks, or data.** Produce a concise, reproducible Markdown report (`sadcp_selection_audit.md`) with the following results. provide this in the current pwd i.e., `/home/z5297792/UNSW-PhD/esp/os_paper_tests`.

## SOLO — Figure 6 - `SOLO_transect_data.ipynb`
- Identify the exact UTC start/end timestamps of **retained S-ADCP observations**, number of records, and selected depth.
- Distinguish the nominal selection window (`t0 = 2023-10-18 01:30`, `dt = 0.5 days`) from actual retained timestamps.
- Report the number of valid observations used for inner-core estimation and outer-core fitting, with the exact spatial/quality filters and units.

## LATTE — Figure 12 - `LATTE_app1_multiplatform_data.ipynb`
- Identify exact UTC start/end timestamps and number of S-ADCP records after the initial `Lat <= -37` filter; confirm depth selected in `latte_source_selector`.
- Report counts of finite S-ADCP vectors and satellite vectors entering the combined input, and the **per-platform counts** retained in the first inner-core selection, refined inner-core selection, and outer-core fit. Include total counts.
- Identify the satellite interpolation time, its source (mean S-ADCP timestamp), and relevant spatial filters.
- Verify the implemented `rho_core` and `rho_outer` defaults against the manuscript's 30 km/100 km claims; flag any mismatch.

## Requirements
- Derive results by executing/reproducing the notebook's actual selection logic using repository data; do not estimate counts from figures or assume that initial records all enter a fit.
- Include UTC timestamps to the nearest second, durations, exact sample counts, and code references (file/function/cell).
- Distinguish **source data**, **time/latitude/depth-filtered data**, **valid vectors**, **inner-core fit**, and **outer-core fit**.
- Note if the saved notebook state differs from executing the current source code, or if data are missing; explicitly label anything not verifiable.
- End with two short manuscript-ready factual sentences, one for SOLO and one for LATTE, using only verified results. No new experiments or methodological changes.
