function [gammaMap, gammaMean, gammaMax, rmse, gpr, rmseMap] = computeGammaMetrics(pred, meas, dta, tta, idt, res_ini, res_end, timestep_ini, timestep_end)
% ************************************************************************
% *  Gamma and RMSE map calculation with DTA and TTA consideration       *
% ************************************************************************
%
% -----------------  Objective  -----------------
% 1 Adjust spatial resolution (res_ini -> res_end)
% 2 Adjust temporal frequency (timestep_ini -> timestep_end)
% 3 Search for pixels in a neighborhood defined by DTA and TTA
% 4 Calculate Gamma according to the true 3D formula
% 5 Calculate RMSE on all available images

% -----------------  Inputs  -----------------
%   pred  : Predicted SSI map ([ny, nx])
%   meas  : Temporal sequence of measured SSI maps ([ny, nx, nt])
%   dta   : Spatial tolerance (km)
%   tta   : Temporal tolerance (minutes)
%   idt   : Intensity deviation tolerance (W/m²)
%   res_ini : Initial pixel resolution (km)
%   res_end : Target resolution after rescaling (km)
%   timestep_ini : Initial temporal interval between images (h)
%   timestep_end : Target temporal interval after rescaling (h)


% -----------------  Spatial rescaling  -----------------
if res_ini ~= res_end
    scaleFactorSpatial = res_end / res_ini; 
    warning('Spatial rescaling applied (Factor: %.2f)', scaleFactorSpatial);
    pred = imresize(pred, scaleFactorSpatial);
    meas = imresize(meas, scaleFactorSpatial);
end

[ny, nx, nt] = size(meas); % Update dimensions after rescaling

% -----------------  Temporal rescaling  -----------------
if timestep_ini ~= timestep_end
    scaleFactorTemporal = timestep_ini / timestep_end;
    nt_new = round(nt * scaleFactorTemporal);
    warning('Temporal rescaling applied (Factor: %.2f, Final images: %d)', scaleFactorTemporal, nt_new);

    % Temporal interpolation
    meas_resized = zeros(ny, nx, nt_new);
    time_original = linspace(0, (nt-1) * timestep_ini * 60, nt); % Time in minutes
    time_interpolated = linspace(0, (nt-1) * timestep_ini * 60, nt_new); % New time values

    for i = 1:ny
        for j = 1:nx
            meas_resized(i, j, :) = interp1(time_original, squeeze(meas(i, j, :)), ...
                                            time_interpolated, 'linear', 'extrap');
        end
    end
    meas = meas_resized;
    nt = nt_new; % Update number of images
end

% ----------------- Initialize result matrices -----------------
gammaMap = inf(ny, nx); % Gamma values map
rmseMap = zeros(ny, nx); % RMSE error map

% ----------------- Neighborhood search (DTA, TTA) and Gamma calculation -----------------
[X, Y] = meshgrid(1:nx, 1:ny);
dmax = 3 * dta / (2 * res_end); % Maximum spatial distance in pixels
tmax = 3 * tta / (2 * timestep_end * 60); % Maximum temporal distance in indices

for i = 1:ny
    for j = 1:nx
        gamma_min = inf; % Initialize minimum gamma

        for t = 1:nt
            %  Define spatial and temporal neighborhood
            x_min = max(1, j - round(dmax));
            x_max = min(nx, j + round(dmax));
            y_min = max(1, i - round(dmax));
            y_max = min(ny, i + round(dmax));
            t_min = max(1, t - round(tmax));
            t_max = min(nt, t + round(tmax));

            %  Extract neighborhood
            [Xt, Yt] = meshgrid(x_min:x_max, y_min:y_max);
            It = meas(y_min:y_max, x_min:x_max, t_min:t_max);

            %  Adjust dimensions for compatibility
            [nb_y, nb_x, nb_t] = size(It);  % Extracted neighborhood size

            % Calculate spatial distance
            d_space = ((j - Xt).^2 + (i - Yt).^2) / dta^2;
            d_space = repmat(d_space, [1, 1, nb_t]); % Extend to 3D

            % Calculate temporal distance with correct dimensions
            t_indices = t_min:t_max;
            d_time = ((t - t_indices).^2) / tta^2;
            % Create 3D array with correct dimensions [nb_y, nb_x, nb_t]
            d_time = reshape(d_time, [1, 1, nb_t]); 
            d_time = repmat(d_time, [nb_y, nb_x, 1]);

            % Calculate intensity distance
            d_intensity = ((pred(i, j) - It).^2) / idt^2;

            %  Avoid Inf values
            d_space(isinf(d_space)) = NaN;
            d_time(isinf(d_time)) = NaN;
            d_intensity(isinf(d_intensity)) = NaN;

            %  Calculate Gamma with compatible dimensions
            gamma_local = sqrt(d_space + d_time + d_intensity);

            %  Correction before taking min()
            gamma_local(isinf(gamma_local)) = NaN;
            gamma_min = min(gamma_min, min(gamma_local(:), [], 'omitnan'));
        end

        %  Update maps
        if isinf(gamma_min)
            gamma_min = NaN; % Final correction to avoid Inf propagation
        end
        gammaMap(i, j) = gamma_min;
        rmseMap(i, j) = sqrt(mean((pred(i, j) - squeeze(meas(i, j, :))).^2, 'omitnan'));
    end
end

% ----------------- Calculate global metrics -----------------
gammaMean = mean(gammaMap(~isinf(gammaMap(:))), 'omitnan');
gammaMax = max(gammaMap(~isinf(gammaMap(:))), [], 'omitnan');
rmse = mean(rmseMap(:), 'omitnan');
gpr = mean(gammaMap(:) <= 1, 'omitnan') * 100;

end
