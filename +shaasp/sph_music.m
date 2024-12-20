function [M, T, P] = sph_music(alpha, num_src, r, k, options)
% SPH_MUSIC - Source localisation from SH coefficients (alphas).
% Spherical Harmonic MUlitple SIgnal Classification source localisation.
%
% Syntax:  [M, T, P] = shmusic(alpha, num_src, r, k, options)
%
% Inputs:
%   alpha - Interior sound field coefficients [N,1], [N,K] or [N,K,T].
%   num_src - Number of sources to localise.
%   r - Steering radius for near-field localisation.
%   k - Steering wave number for near-field localisation.
%
%   Options:
%   ngrid - [360] number of steering directions per (theta, phi) axis.
%   logscale - [true] convert MUSIC spectra to dB log scale.
%   f_plot - [false] plot the MUSIC spectra. 
%   h_fig - figure handle to plot the MUSIC spectra on.
%
% Outputs:
%   M - [T,P] MUSIC spectra for (theta, phi) meshgrid positions.
%   T - Theta meshgrid values for M.
%   P - Phi meshgrid values for M.
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 19-Dec-2024
% Last revision: 19-Dec-2024

arguments
    alpha (:,:,:)
    num_src = 1
    r = []
    k = []
    options.ngrid = 360
    options.logscale = true;
    options.f_plot = false;
    options.h_fig = [];
end

% Truncation order of sound field.
N = sqrt(size(alpha, 1)) - 1;

% Get spatial correlation matrix Ra from alpha.
if ndims(alpha) <=2
    % Alpha [N,1] or [N,K].
    alpha = permute(alpha, [1,3,2]);  % [N,1,K]
    Ra = pagemtimes(alpha, 'none', alpha, 'ctranspose');  % [N,N,K]
else
    % Alpha [N,K,T], average over time dimension.
    alpha = permute(alpha, [1,4,2,3]);  % [N,1,K,T]
    Ra = pagemtimes(alpha, 'none', alpha, 'ctranspose');  % [N,N,K,T]
    Ra = mean(Ra, 4, "omitmissing");  % [N,N,K]
end

% Create localisation / plotting grid.
t = linspace(0, pi, options.ngrid).';
p = linspace(-pi, pi, options.ngrid).';
[T,P] = meshgrid(t,p);
M = zeros(size(T));

% Far-field MUSIC.
if isempty(r) || isinf(r)
    % Frequency smoothing.
    Rsmooth = mean(Ra, 3);
    [U,~,~] = svd(Rsmooth);
    Un = U(:, num_src+1:end);

    % Far-field steering vector.
    y = (-1i).^shaasp.SPHMacros.n_set(N) .* (4*pi) .* conj(shaasp.sph_ynm(N, T(:), P(:)));

    % Far-field MUSIC spectra.
    M(:) = 1 ./ sum(abs( y * Un ).^2, 2);
    M = reshape(M, [length(t), length(p)]);

% Near-field MUSIC.
else
    % Near-field steering vector.
    y = 1i .* permute(k(:),[2,3,1]) .* shaasp.sph_hn(N,k,r) .* conj(shaasp.sph_ynm(N, T(:), P(:)));  % [Q,N,K]

    nbin = size(alpha, 3);
    Mband = zeros(length(T(:)), nbin);
    for ibin = 1:nbin
        % Narrow-band noise sub-space.
        Rband = Ra(:,:,ibin);
        [U,~,~] = svd(Rband);
        Un = U(:, num_src+1:end);
        
        % Near-field MUSIC spectra.
        Mband(:,ibin) = 1 ./ sum(abs( y(:,:,ibin) * Un ).^2, 2);  % v2. correct.
    end

    % Average MUSIC spectra over frequencies.
    M = mean(Mband, 2, 'omitmissing');
    M = reshape(M, [length(t), length(p)]);
    M = M ./ max(abs(M(:)));  % Normalise M.
    if options.logscale
        M = 10 .* log10(M);  % Log scale.
    end
end

% Plot MUSIC spectra.
if ~nargout || options.f_plot
    
    if options.logscale
        cmin = -40;
        M(M < cmin) = cmin;
    else
        cmin = 10^(-40/20);
        M(M < cmin) = cmin;
    end

    if isempty(options.h_fig)
        options.h_fig = figure('Color', [1,1,1]);
    end

    surf(rad2deg(T), rad2deg(P), M, 'EdgeColor', 'interp');
    
    view(-90, -90);
    xlabel('theta');
    xlim(rad2deg([min(T(:)), max(T(:))]));
    ylabel('phi');
    ylim(rad2deg([min(P(:)), max(P(:))]));
    axis('square');
end
end