function [p] = sfd_greens_sphcoord(source_type, k, source_rtp, mic_rtp, s, w)
% [p] = sfd_greens_sphcoord(source_type, k, source_rtp, mic_rtp, s, w)
% Wrapper function for shaasp.sfd_greens().
% Take input as spherical coordinates instead of cartesian coordinates.
% (r,t,p) == (radius, theta, phi) == (radius, elevation, azimuth).
% Theta = elevation angle downwards from +z-axis, 0 to pi.
% Phi = azimuth angle counterclockwise from +x-axis >> +y-axis, 0 to 2pi.

    % Default inputs.
    if (nargin < 6), w = 1; end
    if (nargin < 5), s = 1; end

    % Convert to cartesian.
    if any(strcmp(source_type, {'plane-wave','planewave','pw', ...
                                'planewave-','pw-','planewave+','pw+'}))
        source_xyz = shaasp.rtp2xyz(k, source_rtp(:,2), source_rtp(:,3));
    else
        source_xyz = shaasp.rtp2xyz(source_rtp);
    end
    mic_xyz = shaasp.rtp2xyz(mic_rtp);
    
    % Solve greens.
    p = shaasp.sfd_greens(source_type, k, source_xyz, mic_xyz, s, w);
    
end