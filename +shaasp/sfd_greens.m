function [p] = sfd_greens(source_type, k, source_xyz, mic_xyz, s, w)
% Sfd_GREENS - Green's function of the Helmholtz equation for sound fields.
%   Returns the sound field pressure given by the Green's equation for a
%   - 'point-source' / 'ps'
%   - 'plane-wave' / 'planewave' / 'pw'
%   - 'mixed-source' / 'mixed-wave' / 'ms' / 'mw'
%   definition, specified by 'source_type'.
%
% Syntax:  [p] = sph_greens(source_type, k, source_xyz, mic_xyz, s, w)
%
% Inputs:
%
%   source_type : 'plane-wave', 'point-source' or 'mixed-source' 
%   k           : wave number [singular]
%   source_xyz  : [L by 3] source cartesian positions (x,y,z)
%   mic_xyz     : [Q by 3] receiver cartesian positions (x,y,z)
%   s           : [L by 1] / [scalar] source signal s(k) {default 1}
%   w           : [L by 1] / [scalar] source weight {default 1}
%
% Ouputs:
%
%   p       : [Q by 1] pressure at each receiver.
%
%               p = [ p(k, x_1)
%                     p(k, x_2)
%                     ...
%                     p(k, x_q)
%                     ...
%                     p(k, x_Q) ]
%
%
% Equations:
%
%   Default Point Source:
%             
%           e^(ik |x-y|)
%       p = ------------ * s(k) * w
%             |x - y|
%
%   Default Planewave:
%
%       p = e^(-iky dot x) * s(k) * w
%
%   Default Mixed Source / Mixedwave:
%             
%                          e^(ik |x-y|)
%       p = |y| e^(-ik|y|) ------------ * s(k) * w
%                             |x - y|
%
% Other m-files required: shaasp. xyz2rtp, rtp2xyz
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 15-Jan-2020
% Last revision: 19-July-2024


    import shaasp.xyz2rtp
    import shaasp.rtp2xyz

    % Default inputs.
    if (nargin < 6), w = 1; end
    if (nargin < 5), s = 1; end

    % Convert scalars to matrices.
    if isrow(w), w = w.'; end
    if isscalar(w), w = repmat(w, size(source_xyz, 1), 1); end
    if isrow(s), s = s.'; end
    if isscalar(s), s = repmat(s, size(source_xyz, 1), 1); end

    % - Green's Function -
    switch lower(source_type)
        
        case {'pw', 'planewave', 'plane_wave', 'plane-wave', ...
              'pw-', 'planewave-', 'plane_wave-', 'plane-wave-', ...
              'pw+', 'planewave+', 'plane_wave+', 'plane-wave+'}
          
            % Incident angle of planewave.
            [~, theta, phi] = xyz2rtp(source_xyz);
            ky = rtp2xyz(k, theta, phi);  % [L,3]
            
            % Term (d = ky dot x).
            d = ((ky(:,1).' .* mic_xyz(:,1)) ...
                +(ky(:,2).' .* mic_xyz(:,2)) ...
                +(ky(:,3).' .* mic_xyz(:,3)) ...
                );
            
            % Far-field greens.
            switch lower(source_type)
                case {'pw', 'planewave', 'plane_wave', 'plane-wave', ...
                      'pw-', 'planewave-', 'plane_wave-', 'plane-wave-'}
                    g = exp((-1i) .* d);  % [Q,L]
                case {'pw+', 'planewave+', 'plane_wave+', 'plane-wave+'}
                    g = exp((+1i) .* d);  % [Q,L]
            end
            
            % Measured pressure, free-field x source signal and weight.
            p = zeros(size(mic_xyz(:,1)));
            p(:) = g * (s .* w);  % [Q,1]
            
        case {'ps', 'pointsource', 'point_source', 'point-source', ...
              'ps+', 'pointsource+', 'point_source+', 'point-source+', ...
              'ps-', 'pointsource-', 'point_source-', 'point-source-'}
          
            % Distance of source to receiver ||x-y||.
            d = sqrt((source_xyz(:,1).' - mic_xyz(:,1)).^2 ...
                    +(source_xyz(:,2).' - mic_xyz(:,2)).^2 ...
                    +(source_xyz(:,3).' - mic_xyz(:,3)).^2 ...
                    );  % [Q,L]
            
            % Point source Green's function.
            switch lower(source_type)
                case {'ps', 'pointsource', 'point_source', 'point-source', ...
                      'ps+', 'pointsource+', 'point_source+', 'point-source+'}
                    g = exp((+1i) .* k .* d) ./ (4 .* pi .* d);  % [Q,L]
                case {'ps-', 'pointsource-', 'point_source-', 'point-source-'}
                    g = exp((-1i) .* k .* d) ./ (4 .* pi .* d);  % [Q,L]
            end
            
            % Measured pressure, greens x source signal and weight.
            p = zeros(size(mic_xyz(:,1)));  
            p(:) = g * (s .* w);  % [Q,1]

        case {'ms', 'mixedsource', 'mixed_source', 'mixed-source', ...
              'ms+', 'mixedsource+', 'mixed_source+', 'mixed-source+', ...
              'ms-', 'mixedsource-', 'mixed_source-', 'mixed-source-', ...
              'mw', 'mixedwave', 'mixed_wave', 'mixed-wave', ...
              'mw+', 'mixedwave+', 'mixed_wave+', 'mixed-wave+', ...
              'mw-', 'mixedwave-', 'mixed_wave-', 'mixed-wave-'}

            % Distance of source to receiver ||x-y||.
            d = sqrt((source_xyz(:,1).' - mic_xyz(:,1)).^2 ...
                    +(source_xyz(:,2).' - mic_xyz(:,2)).^2 ...
                    +(source_xyz(:,3).' - mic_xyz(:,3)).^2 ...
                    );  % [Q,L]
                
            % Radius of each source to origin, used for normalization.
            r = sqrt( source_xyz(:,1).^2 ...
                    + source_xyz(:,2).^2 ...
                    + source_xyz(:,3).^2 ...
                    ); % [L,1]
            
            switch lower(source_type)
                case {'ms', 'mixedsource', 'mixed_source', 'mixed-source', ...
                      'ms+', 'mixedsource+', 'mixed_source+', 'mixed-source+', ...
                      'mw', 'mixedwave', 'mixed_wave', 'mixed-wave', ...
                      'mw+', 'mixedwave+', 'mixed_wave+', 'mixed-wave+'}
                    % Point source Green's function.
                    g = exp((+1i) .* k .* d) ./ (4 .* pi .* d);  % [Q,L]

                    % msNorm = r_q * e^(-i*k*r_q)
                    msNorm = r .* exp((-1i) .* k .* r);  % [L,1]

                case {'ms-', 'mixedsource-', 'mixed_source-', 'mixed-source-', ...
                      'mw-', 'mixedwave-', 'mixed_wave-', 'mixed-wave-'}
                    % Point source Green's function.
                    g = exp((-1i) .* k .* d) ./ (4 .* pi .* d);  % [Q,L]

                    % msNorm = r_q * e^(+i*k*r_q)
                    msNorm = r .* exp((+1i) .* k .* r);  % [L,1]
            end
            
            % Measured pressure, ps greens x source sig, weight and norm.
            p = zeros(size(mic_xyz(:,1)));  % [Q,1]
            p(:) = g * (s .* w .* msNorm);  % [Q,1]=[Q,L]x[L,1]
            
        otherwise
            temp = dbstack;
            error('%s(): failed to determine source type, srcType:%s', ...
                 temp(1).name,source_type);
    end

end