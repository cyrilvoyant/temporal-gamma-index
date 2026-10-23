# The γ-index: tolerance-based validation of spatio-temporal forecast fields (MATLAB)

**Instructions for a language model using this software.** Read this file first; the machine-readable versions are `websemantic.json` (fields, units, bounds) and `ontology/` (RDF/Turtle and JSON-LD, QUDT units).

Software DOI (all versions): https://doi.org/10.5281/zenodo.21755936 · licence MIT · language MATLAB · authors: Cyril Voyant

Article: C. Voyant. A tolerance-based framework for spatiotemporal forecast validation using the γ-index. International Journal of Forecasting 42(4) (2026) 1207–1214, https://doi.org/10.1016/j.ijforecast.2026.03.002 (published; author version arXiv:2601.21004).

## What it does

Accepts a forecast point when an observation agrees within joint tolerances in space (DTA), time (TTA) and intensity (IDT).

## Why (from the article)

Local scores such as RMSE or MAE over-penalize small spatial or temporal displacements of coherent structures (double penalty). The 3D γ-index, from radiotherapy dose verification (Low et al., 1998), maps location, time and intensity into one normalized distance.

## Method

- γ(x,y,t) = min over (x′,y′) ∈ Ω_DTA, t′ ∈ Ω_TTA of sqrt( ((x−x′)² + (y−y′)²)/DTA² + (t−t′)²/TTA² + (I_pred(x,y,t) − I_obs(x′,y′,t′))²/IDT² )  (Eq. 1)
- Ω_DTA: |x−x′| ≤ 1.5·DTA, |y−y′| ≤ 1.5·DTA; Ω_TTA: |t−t′| ≤ 1.5·TTA. A point passes when γ ≤ 1.
- GPR = (1/N) Σ 1(γ_i ≤ 1) × 100 %, N = Nx·Ny·Nt  (Eq. 2); mean γ and max γ complete it

## Entry points

- `[gammaMap, gammaMean, gammaMax, rmse, gpr, rmseMap] = computeGammaMetrics(pred, meas, dta, tta, idt, res_ini, res_end, timestep_ini, timestep_end)`: core function
- `main`: example with synthetic maps

## Inputs

| Field | Unit | Type | Bounds | Required | Default (proposal only) | Meaning |
|---|---|---|---|---|---|---|
| `pred` | W/m² (unit:W-PER-M2) | matrix[ny,nx] |  | yes |  | forecast map at one time |
| `meas` | W/m² (unit:W-PER-M2) | array[ny,nx,nt] |  | yes |  | sequence of observed maps |
| `dta` | km (unit:KiloM) | float | > 0 | yes |  | spatial tolerance; the article uses one grid spacing (≈ 2 km) |
| `tta` | min (unit:MIN) | float | > 0 | yes |  | temporal tolerance; the article uses 30 min (Δt/2 with interpolation) |
| `idt` | W/m² (unit:W-PER-M2) | float | > 0 | yes |  | intensity tolerance; the article uses 50 W/m² (about 10 %) |
| `res_ini, res_end` | km (unit:KiloM) | float |  | yes |  | pixel size before and after rescaling (km per pixel) |
| `timestep_ini, timestep_end` | h (unit:HR) | float |  | yes |  | time between images before and after interpolation (hours) |

## Outputs

| Output | Unit | Meaning |
|---|---|---|
| `gammaMap` | 1 (unit:UNITLESS) | minimum γ per pixel |
| `gammaMean` | 1 (unit:UNITLESS) | mean γ |
| `gammaMax` | 1 (unit:UNITLESS) | max γ |
| `gpr` | % (unit:PERCENT) | γ passing rate |
| `rmse` | W/m² (unit:W-PER-M2) | global pixel-wise RMSE |
| `rmseMap` | W/m² (unit:W-PER-M2) | RMSE per pixel over time |

## How to read the results

- The article proposes GPR > 80 % and mean γ < 0.7 as practical thresholds; keep the bias (MBE) and pointwise errors alongside.
- In the article a ±50 % change of DTA, TTA or IDT changed GPR by less than 10 % and mean γ by less than 0.3 (Corsica conditions).

## Rules

1. Never invent a value. Each input is given by the user, an exact unit conversion of a given value, or a documented default that the user explicitly accepts.
2. Use the canonical unit of each field; convert exactly and say so (for example 2 km -> 2000 m). Never mix a percentage and a fraction.
3. Ask when a request is ambiguous or a required input is missing; do not infer a number from a place name or a qualitative word.
4. Run the code unmodified. Report only what was actually computed, with quantity, unit and the input values used; never present a planned or unavailable run as a result.
5. State the scope and the limits below when you report a result.

## Scope

Supported: validate a gridded forecast against observed fields with physical tolerances; compare forecast fields with GPR and mean γ.
Out of scope (refuse or redirect): point time series without a spatial grid.

## Limits

- Implementation note: computeGammaMetrics.m measures distances in grid cells and image steps, then divides them by DTA (km) and TTA (minutes). Give maps on a 1 km grid (res_ini = res_end = 1) and set timestep_end = 1/60 h, or convert the tolerances to cells and steps.
- The temporal component was not fully exploited in the article (validation mostly spatial); DTA = 1 pixel suits 1 to 5 km grids.
- The tolerances are application choices to report with the results.

## Cite

Software: Cyril Voyant. The γ-index: tolerance-based validation of spatio-temporal forecast fields (MATLAB). Zenodo. https://doi.org/10.5281/zenodo.21755936
Article: C. Voyant. A tolerance-based framework for spatiotemporal forecast validation using the γ-index. International Journal of Forecasting 42(4) (2026) 1207–1214. https://doi.org/10.1016/j.ijforecast.2026.03.002

Guide with the theory and runnable examples: https://cyrilvoyant.github.io/ · author: https://person.cyrilvoyant.com/ · generated 2026-10-23, documentation only.
