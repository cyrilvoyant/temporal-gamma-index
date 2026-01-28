% **********************************************************************
% *       Initialization and Gamma/RMSE map calculation                *
% **********************************************************************
%
% This script:
% 1 Generates measurement maps (meas_maps) with realistic dynamics
% 2 Generates a predicted map (pred_map) with some uncertainty
% 3 Applies spatial and temporal rescaling if necessary
% 4 Computes Gamma and RMSE considering 12 images before/after
% 5 Displays results as maps
%
% -----------------  Main parameters  -----------------

%  Tolerances for Gamma calculation
dta = 50;    % Spatial tolerance (km) 
idt = 50;   % Intensity deviation tolerance (W/m²)
tta = 30;   % Temporal tolerance (minutes)

%  Spatial and temporal resolution parameters
res_ini = 1;     % Initial pixel resolution (km)
res_end = 1;     % Target resolution after rescaling (km)
timestep_ini = 1;  % Initial interval between images (h)
timestep_end = 1; % Target interval after rescaling (h)

%  Image dimensions
ny = 50; % Number of rows
nx = 50; % Number of columns
nt = 24;  % Number of temporal maps (currently 12 before/after)

% ----------------------------------------------------------------------
% 1 Generation of measurement maps and prediction map
% ----------------------------------------------------------------------

% Coordinate grid to simulate spatial gradient
[x, y] = meshgrid(linspace(-1, 1, nx), linspace(-1, 1, ny));

% Base map (ideal measurements without noise)
meas_base = 800 + 200 * peaks(ny) .* exp(-(x.^2 + y.^2)/2);

% Initialize maps
meas_maps = zeros(ny, nx, nt); % Measured maps (one per time t)
pred_map = zeros(ny, nx);      % Predicted map

% ----------------------------------------------------------------------
% 2 Generate measurement maps with realistic noise
% ----------------------------------------------------------------------

noise_std = 10;  % Noise standard deviation (W/m²), adjusts uncertainty

for t = 1:nt
    noise = randn(ny, nx) * noise_std/2; % Generate Gaussian noise
    meas_maps(:,:,t) = meas_base + noise; % Add noise

    % Manage out-of-range values
    meas_maps(:,:,t) = min(meas_maps(:,:,t), 1000); 
    meas_maps(:,:,t) = max(meas_maps(:,:,t), 0);
end

% ----------------------------------------------------------------------
% 3 Generate predicted map
% ----------------------------------------------------------------------

%  Create smoothing filter to simulate prediction
filter = fspecial('average', 5);  % 5x5 averaging filter

%  Apply filter to base map (spatial smoothing)
pred_map = conv2(meas_base, filter, 'same');

%  Add realistic noise to prediction
pred_map = pred_map + randn(ny, nx) * noise_std/3;  % Add reduced noise

%  Manage out-of-range values
pred_map = min(pred_map, 1000);
pred_map = max(pred_map, 0);

% ----------------------------------------------------------------------
% 4 Select temporal images (up to 12 before and 12 after)
% ----------------------------------------------------------------------
if nt == 1
    warning('Only 1 image available in meas_maps. Temporal component will not be used.');
    meas_selected = meas_maps; % No temporal selection possible
else
    t_center = round(nt / 2);
    t_min = max(1, t_center - 12);
    t_max = min(nt, t_center + 12);
    meas_selected = meas_maps(:, :, t_min:t_max);
end

% ----------------------------------------------------------------------
% 5 Call computeGammaMetrics to calculate Gamma and RMSE
% ----------------------------------------------------------------------
[gammaMap, gammaMean, gammaMax, rmse, gpr, rmseMap] = ...
    computeGammaMetrics(pred_map, meas_selected, dta, tta, idt, res_ini, res_end, timestep_ini, timestep_end);

% ----------------------------------------------------------------------
% 6 Display Gamma and RMSE results
% ----------------------------------------------------------------------

figure;

%  Measurement map (first selected time step)
subplot(2,2,1);
imagesc(meas_selected(:,:,1));
colorbar;
title('Measurement Map', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$X$ Position (pixels)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Y$ Position (pixels)', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 14, 'FontName', 'Helvetica');
caxis([700 1000]);

%  Prediction map
subplot(2,2,2);
imagesc(pred_map);
colorbar;
title('Prediction Map', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$X$ Position (pixels)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Y$ Position (pixels)', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'FontSize', 14, 'FontName', 'Helvetica');
caxis([700 1000]);

%  Gamma map with contour for gamma > 1
subplot(2,2,3);
imagesc(gammaMap);
colorbar;
caxis([0 2]); 
title(['Gamma Map $(DTA=', num2str(dta), '$ km, $TTA=', num2str(tta), '$ min, $IDT=', num2str(idt), '$ W/m$^2$)'], ...
    'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$X$ Position (pixels)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Y$ Position (pixels)', 'Interpreter', 'latex', 'FontSize', 14);
hold on;
xlim([0, 50]);
ylim([0, 50]);
contour(gammaMap, [1 1], 'k', 'LineWidth', 2); % Contour for gamma > 1
hold off;
set(gca, 'FontSize', 14, 'FontName', 'Helvetica');

% RMSE map with contour for gamma > 1
subplot(2,2,4);
imagesc(rmseMap);
colorbar;
title('RMSE Map', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('$X$ Position (pixels)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Y$ Position (pixels)', 'Interpreter', 'latex', 'FontSize', 14);
hold on;
xlim([0, 50]);
ylim([0, 50]);
contour(gammaMap, [1 1], 'k', 'LineWidth', 2); % Reuse contour for gamma > 1
caxis([min(rmseMap(:)), 50]);
hold off;
set(gca, 'FontSize', 14, 'FontName', 'Helvetica');

% Calculate and display the percentage of area where gamma > 1
gammaArea = sum(gammaMap(:) > 1) / numel(gammaMap) * 100;

% Add global statistics with actual calculated values
annotation('textbox', [0.13, 0.02, 0.75, 0.06], 'String', ...
    sprintf('Mean \\gamma: %.3f, Max \\gamma: %.2f, GPR: %.1f%%, Area (\\gamma > 1): %.2f%%, TTA: %d min', ...
    gammaMean, gammaMax, gpr, gammaArea, tta), ...
    'Interpreter', 'tex', 'FontSize', 14, 'FontName', 'Helvetica', ...
    'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'BackgroundColor', 'white');
