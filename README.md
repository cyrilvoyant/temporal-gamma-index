# Gamma Index Evaluation for Solar Irradiance Forecasts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21755937.svg)](https://doi.org/10.5281/zenodo.21755937)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:9c56a454a1d7e6f5efa05c3652fbf4e15faf9fc9/)](https://archive.softwareheritage.org/swh:1:dir:9c56a454a1d7e6f5efa05c3652fbf4e15faf9fc9;origin=https://github.com/cyrilvoyant/temporal-gamma-index;visit=swh:1:snp:9794f4366ae14dfb227ebba68aedd2faf21d2c8e;anchor=swh:1:rev:e0c451738c92325dc90f877a84b7fb2d8a6ceabe;path=/)
A MATLAB implementation of a **tolerance-based, 3D Gamma Index framework** for validating gridded (spatio-temporal) forecasts, adapting a method originally developed for radiotherapy dose verification (Low et al., 1998) to the double-penalty problem in map forecast verification: classical pixel-wise scores such as RMSE/MAE over-penalize forecasts that are structurally correct but slightly displaced in space or time. The Gamma Index instead asks whether a prediction falls within explicit physical tolerances — spatial (DTA), temporal (TTA), and intensity (IDT) — of *some* nearby observation, rather than the exact same pixel at the exact same instant. Companion code for the manuscript "A Tolerance-Based Framework for Spatio-Temporal Forecast Validation Using the Γ-Index" (see [`Gamma_Index_Cyril_Voyant.pdf`](Gamma_Index_Cyril_Voyant.pdf), currently anonymized for peer review — see [Citation](#citation)).

## Overview

This package implements a comprehensive Gamma Index analysis tool specifically designed for evaluating solar irradiance forecasts against measurements. The Gamma Index is widely used in validation studies to assess prediction accuracy by considering spatial, temporal, and intensity deviations simultaneously.

### What is the Gamma Index?

The Gamma Index (Γ) is a composite metric that evaluates how well a prediction matches measurements by searching for the closest match within defined tolerance criteria:

- **DTA (Distance-To-Agreement)**: Spatial tolerance in kilometers
- **TTA (Time-To-Agreement)**: Temporal tolerance in minutes  
- **IDT (Intensity Deviation Tolerance)**: Allowable intensity difference in W/m²

For a prediction at grid coordinates $(x,y,t)$, the 3D Gamma Index is the minimum, over a spatio-temporal search neighborhood, of the normalized joint distance to the observed field:

$$\Gamma(x,y,t)=\min_{(x',y',t')\in\mathcal N}\sqrt{\frac{(x-x')^2+(y-y')^2}{\mathrm{DTA}^2}+\frac{(t-t')^2}{\mathrm{TTA}^2}+\frac{\big(\mathrm{pred}(x,y,t)-\mathrm{obs}(x',y',t')\big)^2}{\mathrm{IDT}^2}}$$

where the search neighborhood $\mathcal N$ extends slightly beyond the nominal tolerances ($\tfrac32\,\mathrm{DTA}$ spatially, $\tfrac32\,\mathrm{TTA}$ temporally) to guarantee a feasible match is found. A pixel "passes" the test if $\Gamma\le1$. Two aggregate statistics summarize the whole field:

$$\mathrm{GPR}=\frac{1}{N}\sum_{i=1}^{N}\mathbb 1(\Gamma_i\le1)\times100\%,\qquad \bar\Gamma=\frac{1}{N}\sum_{i=1}^N \Gamma_i,\qquad \Gamma_{\max}=\max_i \Gamma_i$$

GPR (Gamma Pass Rate) is the percentage of the field within tolerance; the mean and max Gamma complement it by quantifying how far outlier pixels fall from the acceptance boundary. Sensitivity analysis in the companion paper shows that a ±50% change in DTA, TTA, or IDT changes GPR by less than 10 percentage points and mean Γ by less than 0.3, on the HelioClim-3/Corsica dataset studied — i.e. the diagnosis is not overly sensitive to the exact tolerance choice.

## Files

- **`computeGammaMetrics.m`**: Core function that calculates Gamma Index and RMSE maps
- **`main.m`**: Example script demonstrating usage with synthetic data
- **`Gamma_Index_Cyril_Voyant.pdf`**: Companion manuscript with the full derivation, sensitivity analysis, and the HelioClim-3/CAMS satellite validation over Corsica (see [Citation](#citation))
- **`README.md`**: This documentation file

## Key Features

✓ **3D Gamma Index calculation** with spatial, temporal, and intensity components  
✓ **Spatial and temporal rescaling** to match different data resolutions  
✓ **Multiple output metrics**: Gamma maps, RMSE maps, GPR (Gamma Pass Rate), mean/max gamma  
✓ **Visualization tools** with publication-quality figures  
✓ **Flexible tolerance parameters** adaptable to different forecast horizons

## Installation

1. Download all `.m` files to your MATLAB working directory
2. Ensure you have MATLAB R2020a or later with Image Processing Toolbox
3. Run `main.m` to verify installation with synthetic data

## Quick Start

```matlab
% Load or prepare your data
pred_map = ... % Your predicted SSI map [ny, nx]
meas_maps = ... % Your measured SSI maps [ny, nx, nt]

% Set tolerance parameters (example for short-term forecast)
dta = 50;    % 50 km spatial tolerance
tta = 30;    % 30 minutes temporal tolerance
idt = 50;    % 50 W/m² intensity tolerance

% Set resolution parameters
res_ini = 1;         % Initial resolution: 1 km/pixel
res_end = 1;         % Target resolution: 1 km/pixel
timestep_ini = 1;    % Initial timestep: 1 hour
timestep_end = 1;    % Target timestep: 1 hour

% Compute Gamma Index
[gammaMap, gammaMean, gammaMax, rmse, gpr, rmseMap] = ...
    computeGammaMetrics(pred_map, meas_maps, dta, tta, idt, ...
                        res_ini, res_end, timestep_ini, timestep_end);

% Display results
fprintf('Mean Gamma: %.3f\n', gammaMean);
fprintf('Max Gamma: %.3f\n', gammaMax);
fprintf('GPR (Gamma Pass Rate): %.1f%%\n', gpr);
fprintf('RMSE: %.2f W/m²\n', rmse);
```

## Function Reference

### computeGammaMetrics

```matlab
[gammaMap, gammaMean, gammaMax, rmse, gpr, rmseMap] = ...
    computeGammaMetrics(pred, meas, dta, tta, idt, res_ini, res_end, timestep_ini, timestep_end)
```

**Inputs:**
- `pred` [ny × nx]: Predicted SSI map at a single time point
- `meas` [ny × nx × nt]: Temporal sequence of measured SSI maps
- `dta` (scalar): Spatial tolerance in km
- `tta` (scalar): Temporal tolerance in minutes
- `idt` (scalar): Intensity deviation tolerance in W/m²
- `res_ini` (scalar): Initial pixel resolution in km
- `res_end` (scalar): Target pixel resolution in km (for rescaling)
- `timestep_ini` (scalar): Initial time interval between images in hours
- `timestep_end` (scalar): Target time interval in hours (for rescaling)

**Outputs:**
- `gammaMap` [ny × nx]: Map of minimum gamma values for each pixel
- `gammaMean` (scalar): Mean gamma value across all pixels
- `gammaMax` (scalar): Maximum gamma value observed
- `rmse` (scalar): Global RMSE across all pixels
- `gpr` (scalar): Gamma Pass Rate - percentage of pixels with γ ≤ 1
- `rmseMap` [ny × nx]: RMSE map for each pixel across time

**Processing Steps:**
1. Applies spatial rescaling if `res_ini ≠ res_end`
2. Applies temporal rescaling via interpolation if `timestep_ini ≠ timestep_end`
3. For each pixel in the prediction:
   - Defines a 3D search neighborhood based on DTA and TTA
   - Calculates spatial, temporal, and intensity distances
   - Computes local gamma values: γ = √(d_spatial + d_temporal + d_intensity)
   - Records the minimum gamma value found
4. Computes RMSE for each pixel across all time steps
5. Calculates global statistics (mean γ, max γ, GPR)

## Understanding the Example (main.m)

The example script uses **synthetic data** for demonstration purposes. This is intentional:

### Why Synthetic Data?

Real solar irradiance data is proprietary and cannot be freely distributed. The synthetic data in `main.m` is designed to:

1. **Demonstrate proper data formats** required by `computeGammaMetrics`
2. **Show realistic value ranges** (700-1000 W/m² for SSI)
3. **Illustrate the workflow** from data preparation to visualization
4. **Serve as a template** you can adapt for your own data

### What the Example Shows

```matlab
% Section 1-2: Generate measurement maps with noise
meas_maps = zeros(50, 50, 24);  % 50×50 pixels, 24 time steps
for t = 1:24
    meas_maps(:,:,t) = base_map + noise;  % Add Gaussian noise
end

% Section 3: Generate prediction with smoothing and uncertainty
pred_map = filtered(base_map) + reduced_noise;

% Section 4: Select temporal window (12 before/after target time)
meas_selected = meas_maps(:, :, t_center-12:t_center+12);

% Section 5: Call computeGammaMetrics (THIS IS THE KEY PART)
[gammaMap, gammaMean, gammaMax, rmse, gpr, rmseMap] = ...
    computeGammaMetrics(pred_map, meas_selected, dta, tta, idt, ...
                        res_ini, res_end, timestep_ini, timestep_end);

% Section 6: Visualize results
% Creates a 2×2 figure with measurement, prediction, gamma, and RMSE maps
```

### Adapting for Your Data

To use your own solar irradiance data:

```matlab
% Replace synthetic generation with your data loading
pred_map = load_your_prediction_data();     % [ny, nx]
meas_maps = load_your_measurement_data();   % [ny, nx, nt]

% Adjust tolerance parameters for your use case
dta = 50;    % Adjust based on forecast horizon
tta = 30;    % Adjust based on temporal resolution
idt = 50;    % Adjust based on acceptable error

% Rest of the workflow remains the same
[gammaMap, ...] = computeGammaMetrics(pred_map, meas_maps, ...);
```

## Parameter Selection Guidelines

### Tolerance Parameters

The choice of DTA, TTA, and IDT depends on your forecast type and acceptable error levels (for exemple):

| Forecast Horizon | DTA (km) | TTA (min) | IDT (W/m²) | Application |
|------------------|----------|-----------|------------|-------------|
| 0-3 hours (nowcast) | 20-40 | 15-30 | 30-60 | Real-time control |
| 3-6 hours (short-term) | 40-60 | 30-60 | 50-80 | Grid management |
| 6-24 hours (medium-term) | 60-100 | 60-120 | 80-120 | Day-ahead planning |
| 24-72 hours (long-term) | 100-200 | 120-240 | 120-200 | Energy trading |

**Rule of thumb:** A good forecast should achieve GPR > 90% with appropriate tolerances.

### Resolution Parameters

- Set `res_ini` and `res_end` to match your data's spatial resolution
- Set `timestep_ini` and `timestep_end` to match temporal resolution
- Use rescaling (res_ini ≠ res_end) when comparing datasets with different resolutions
- Temporal rescaling uses linear interpolation to match frequencies

## Output Interpretation

### Gamma Map (gammaMap)
- **γ < 1**: Pixel passes criteria (good agreement)
- **γ = 1**: Pixel exactly at tolerance boundary
- **γ > 1**: Pixel fails criteria (poor agreement)
- Lower gamma values indicate better agreement

### RMSE Map (rmseMap)
- Shows temporal variability for each pixel
- Complements gamma analysis with absolute error metric
- Units are W/m² for solar irradiance

### Global Metrics
- **GPR (Gamma Pass Rate)**: Percentage of pixels with γ ≤ 1
  - Target: > 90% for good forecasts
  - 80-90%: Acceptable for longer horizons
  - < 80%: Indicates systematic forecast issues
  
- **Mean γ**: Average gamma across all pixels
  - Target: < 0.5 for excellent forecasts
  - 0.5-0.8: Good forecasts
  - > 1.0: Poor forecast quality

## Visualization

The example generates a 4-panel figure:
1. **Measurement Map**: Shows actual observed irradiance
2. **Prediction Map**: Shows forecast irradiance
3. **Gamma Map**: Color-coded gamma values with γ=1 contour
4. **RMSE Map**: Temporal error distribution with γ=1 contour overlay

Black contour lines indicate regions where γ > 1 (failed criteria).

## Tips for Real Data

1. **Temporal Alignment**: Ensure prediction and measurement times are synchronized
2. **Spatial Registration**: Verify both datasets use the same coordinate system
3. **Data Quality**: Remove invalid measurements (clouds, sensor errors) before analysis
4. **Temporal Window**: Include sufficient images before/after (±12 recommended)
5. **Resolution Matching**: Use rescaling features if datasets have different resolutions
6. **Parameter Tuning**: Start with conservative tolerances and adjust based on results

## Common Issues and Solutions

**Issue**: GPR very low (< 50%)
- Solution: Check if tolerances are too strict; verify data alignment

**Issue**: All gamma values are Inf or NaN
- Solution: Check data ranges; ensure measurements exist for all time steps

**Issue**: Edge effects visible in gamma map
- Solution: This is normal due to limited neighborhood at boundaries

**Issue**: Computation very slow for large datasets
- Solution: Consider spatial downsampling or process regions separately

## Citation

If you use this code, please cite the software itself (see [`CITATION.cff`](CITATION.cff)) and, in academic work, the companion manuscript:

> [Author names anonymized for peer review], "A Tolerance-Based Framework for Spatio-Temporal Forecast Validation Using the Γ-Index," preprint submitted to an Elsevier journal, 2026 (see [`Gamma_Index_Cyril_Voyant.pdf`](Gamma_Index_Cyril_Voyant.pdf)).

This manuscript is currently under peer review; the reference above will be updated with the final author list, journal, and DOI once published. The Γ-Index method itself originates from radiotherapy dose verification (Low, Harms, Mutic, Purdy, 1998, *Medical Physics*).

## License

The **code** in this repository (`computeGammaMetrics.m`, `main.m`) is MIT-licensed — see [`LICENSE`](LICENSE). This applies to the software only: the satellite datasets used in the companion paper (HelioClim-3, CAMS Radiation Service) are third-party and not redistributed here; `main.m` uses synthetic data for exactly this reason (see [Understanding the Example](#understanding-the-example-mainm)).

## Support

For questions about:
- **Code usage**: Refer to the example in `main.m`
- **Parameter selection**: See guidelines above
- **Adapting for your data**: Follow the template in "Adapting for Your Data" section

## Version History

- **v1.0** (2026-01-28): Initial release with 3D Gamma Index implementation
  - Core gamma calculation with DTA/TTA/IDT
  - Spatial and temporal rescaling
  - Comprehensive visualization tools
  - Example with synthetic data

---

**Note**: The synthetic data in `main.m` is for demonstration only. Real solar irradiance data should be obtained from appropriate meteorological sources or satellite observations.
