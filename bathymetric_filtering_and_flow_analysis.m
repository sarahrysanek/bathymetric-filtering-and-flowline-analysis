%% Load DEM (FULL resolution retained for filtering + flow)
[Z_full, R] = readgeoraster('combined_filled_dem.tif');

Z_full = double(Z_full);
Z_full(isnan(Z_full)) = 0;

dx_full = 40; % meters

%% IMPORTANT NOTE
% The DEM is downsampled ONLY for spectral (FFT) analysis.
% All filtering and flow routing use the original 40 m DEM.

%% Create downsampled copy for FFT

downsample_factor = 2;

Z_fft = Z_full;
dx_fft = dx_full;

if downsample_factor > 1
    Z_fft = Z_fft(1:downsample_factor:end, 1:downsample_factor:end);
    dx_fft = dx_fft * downsample_factor;
end

fprintf('FFT resolution: %.1f m | Full resolution: %.1f m\n', dx_fft, dx_full);

%% Compute power spectrum (FFT uses Z_fft ONLY)

[nrows, ncols] = size(Z_fft);

FFTbathy = fftshift(fft2(Z_fft));
power = abs(FFTbathy).^2;

fx = (-floor(ncols/2):ceil(ncols/2)-1) / (dx_fft * ncols);
fy = (-floor(nrows/2):ceil(nrows/2)-1) / (dx_fft * nrows);
[Fx, Fy] = meshgrid(fx, fy);
Fr = sqrt(Fx.^2 + Fy.^2);

Fr_flat = Fr(:);
P_flat  = power(:);

%% Raw power spectrum

figure
subplot(2,1,1)
loglog(Fr_flat, P_flat, '.')
xlabel('Radial Frequency (1/m)')
ylabel('Power')
title('Raw Power Spectrum')
grid on

subplot(2,1,2)
loglog(1./Fr_flat, P_flat, '.')
xlabel('Wavelength (m)')
ylabel('Power')
title('Raw Power Spectrum (Wavelength)')
grid on

%% Normalize spectrum (remove power-law trend)

valid = (Fr_flat > 0) & (P_flat > 0);
Fr_valid = Fr_flat(valid);
P_valid  = P_flat(valid);

logF = log(Fr_valid);
logP = log(P_valid);

G = [ones(size(logF)), logF];
model = G \ logP;

fit = G * model;
residual = logP - fit;

normalized_power = exp(residual);

fprintf('Power-law slope: %.2f\n', model(2));

figure
plot(logF, logP, '.')
hold on
[logF_sorted, idx] = sort(logF);
plot(logF_sorted, fit(idx), 'r-', 'LineWidth', 1.5)
xlabel('log(Frequency)')
ylabel('log(Power)')
title('Power Spectrum with Fit')
grid on

%% Radially averaged spectrum (with uncertainty shading)

num_bins = 500;
edges = logspace(log10(min(Fr_valid)), log10(max(Fr_valid)), num_bins+1);
centers = sqrt(edges(1:end-1).*edges(2:end));

mean_power = nan(num_bins,1);
std_power  = nan(num_bins,1);

for i = 1:num_bins
    idx = Fr_valid >= edges(i) & Fr_valid < edges(i+1);
    if any(idx)
        mean_power(i) = mean(normalized_power(idx));
        std_power(i)  = std(normalized_power(idx));
    end
end

% Remove empty bins
valid_bins = ~isnan(mean_power);

centers     = centers(valid_bins);
mean_power  = mean_power(valid_bins);
std_power   = std_power(valid_bins);

%% Plot radially averaged spectrum

wavelength_bins = 1 ./ centers;

% ensure column vectors
wavelength_bins = wavelength_bins(:);
mean_power      = mean_power(:);
std_power       = std_power(:);

upper = mean_power + std_power;
lower = max(mean_power - std_power, 1e-10);

figure
semilogx(wavelength_bins, mean_power, 'b-', 'LineWidth', 1.5)
hold on

patch([wavelength_bins; flipud(wavelength_bins)], ...
      [upper; flipud(lower)], ...
      'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none')

xlabel('Radial Wavelength (m)')
ylabel('Normalized Power')
title('Radially Averaged Detrended Power Spectrum')

xline(1000, 'r--', 'LineWidth', 1.5, 'Label', '1 km')
xline(5000, 'b--', 'LineWidth', 1.5, 'Label', '5 km')

xlim([0 10000])
grid on
legend('Mean power', '±1 std dev')

%% NOTE:
% Large Gaussian kernels introduce edge artifacts near DEM boundaries.
% To minimize boundary contamination:
%   1. NaNs are temporarily filled during filtering
%   2. Original NaN regions are restored afterward
%   3. Edge regions within ~3 sigma of the boundary are masked
% Flow routing is performed only on the preserved interior region.

%% High-pass filtering (FULL resolution)

cutoffs = [7000, 3000, 500]; % meters

filtered_DEMs = struct();

for i = 1:length(cutoffs)

    cutoff = cutoffs(i);

    sigma = cutoff / (2 * dx_full);

    % Temporary filled DEM for filtering
    Z_filter = Z_full;
    Z_filter(isnan(Z_filter)) = 0;

    % Gaussian filter
    lowpass = imgaussfilt(Z_filter, sigma, ...
        'Padding', 'replicate');

    % High-pass
    highpass = Z_filter - lowpass;

    % Restore original NaNs
    highpass(isnan(Z_full)) = NaN;

    % Remove edge-affected zone
    buffer = ceil(3 * sigma);

    highpass(1:buffer,:) = NaN;
    highpass(end-buffer+1:end,:) = NaN;
    highpass(:,1:buffer) = NaN;
    highpass(:,end-buffer+1:end) = NaN;

    % Store
    fieldname = sprintf('hp_%dm', cutoff);
    filtered_DEMs.(fieldname) = highpass;

    fprintf('Created high-pass DEM: %s\n', fieldname);

end

%% Export filtered DEMs

fields = fieldnames(filtered_DEMs);

for i = 1:length(fields)

    name = fields{i};
    Z_out = filtered_DEMs.(name);

    filename = sprintf('filtered_dem_%s.tif', name);

    geotiffwrite(filename, Z_out, R, 'CoordRefSysCode', 32616);

    fprintf('Saved: %s\n', filename);
end

%% Flow routing for each filtered DEM

flow_results = struct();

for i = 1:length(fields)
    
    name = fields{i};
    filename = sprintf('filtered_dem_%s.tif', name);
    
    fprintf('Running flow analysis for %s...\n', name);
    
    DEM = GRIDobj(filename);
    
    % Invert bathymetry
    DEM.Z = -DEM.Z;
    DEM.Z(abs(DEM.Z) < 1e-10) = NaN;
    
    % Optional smoothing
    %DEM = filter(DEM, 'mean', [3 3]);
    
    FD = FLOWobj(DEM);
    A = flowacc(FD);
    
    % Store results
    flow_results.(name).DEM = DEM;
    flow_results.(name).A   = A;

    A_plot = A;
    A_plot.Z(A_plot.Z <= 0) = NaN;
    
    % Plot
    figure
    imageschs(DEM, log(A_plot))
    title(sprintf('Flow accumulation (%s)', name))
end

%% Preview stream networks before export

% Minimum contributing area threshold
area_threshold = 1e3;

fields = fieldnames(flow_results);

for i = 1:length(fields)

    name = fields{i};

    fprintf('Previewing flow network for %s...\n', name);

    DEM = flow_results.(name).DEM;
    A   = flow_results.(name).A;

    % Extract stream network
    S = STREAMobj(FLOWobj(DEM), A > area_threshold);

    % Optional: keep largest connected component
    %S = klargestconncomps(S, 1);

    % Plot
    figure

    % Background DEM
    imageschs(DEM)
    hold on

    % Overlay stream network
    plot(S, 'b', 'LineWidth', 1)

    title(sprintf('Flow Network Preview (%s)', name))

end
%% Export flow line networks (for QGIS)

% Minimum contributing area threshold (tune this!)
area_threshold = 1e5;  % in number of cells (adjust based on your DEM size)

fields = fieldnames(flow_results);

for i = 1:length(fields)
    
    name = fields{i};
    
    fprintf('Extracting flow network for %s...\n', name);
    
    DEM = flow_results.(name).DEM;
    A   = flow_results.(name).A;
    
    % Extract stream network
    S = STREAMobj(FLOWobj(DEM), A > area_threshold);
    % S = klargestconncomps(S, 1); % keep largest connected network (optional)
    
    % Convert to map coordinates
    GS = STREAMobj2shape(S, 'type', 'map');
    
    % Export to shapefile
    shapename = sprintf('flow_lines_1e5_%s.shp', name);
    shapewrite(GS, shapename);
    
    fprintf('Saved flow lines: %s\n', shapename);
end
