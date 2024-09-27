function [p] = sfd_greens(source_type, k, source_xyz, mic_xyz, s, w)
% Sfd_GREENS - Green's function of the Helmholtz equation for sound fields.
%   Returns the sound field pressure given by the Green's equation for a
%   - 'pointsource' / 'ps'
%   - 'planewave' / 'pw'
%   - 'mixedsource' / 'mixedwave' / 'ms' / 'mw'
%   definition, specified by 'source_type'.
%
% Syntax:  [p] = sph_greens(source_type, k, source_xyz, mic_xyz, s, w)
%
% Inputs:
%
%   source_type : 'planewave', 'pointsource' or 'mixedsource' 'mixedwave' 
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
% Last revision: 27-Sept-2024

    arguments
        source_type {mustBeMember(source_type, ...
            ["pw", "planewave", "pw-", "planewave-", "pw+", "planewave+", ...
             "ps", "pointsource", "ps-", "pointsource-", "ps+", "pointsource+", ...
             "ms", "mixedsource", "ms-", "mixedsource-", "ms+", "mixedsource+", ...
             "mw", "mixedwave", "mw-", "mixedwave-", "mw+", "mixedwave+"])}
        k {mustBeScalarOrEmpty}
        source_xyz (:,3)
        mic_xyz (:,3)
        s (:,1) = 1
        w (:,1) = 1
    end

    % Convert scalars to matrices.
    if isscalar(w)
        w = repmat(w, size(source_xyz, 1), 1); 
    end
    if isscalar(s)
        s = repmat(s, size(source_xyz, 1), 1);
    end

    % - Green's Function -
    switch lower(source_type)
        
        % Planewave source.
        case {'pw', 'planewave', 'pw-', 'planewave-', 'pw+', 'planewave+'}         
            
            % Incident angle of planewave.
            [~, src_theta, src_phi] = shaasp.xyz2rtp(source_xyz);
            ky_hat = shaasp.rtp2xyz(k, src_theta, src_phi);  % [L by 3]
            
            % Distance between each src and mic (d = ky dot x) [Q by L].
            d = mic_xyz * ky_hat.';
            
            if contains(source_type, '+')
                pw_sign = +1i;
            else
                pw_sign = -1i;
            end

            g = exp(pw_sign .* d);  % Free-field Green's function [Q by L]
            p = g * (s .* w);  % Measured pressure [Q by 1]
            return;
            
        % Point source.
        case {'ps', 'pointsource', 'ps+', 'pointsource+', 'ps-', 'pointsource-'}
            
            % Distance of source to receiver ||x-y|| [Q by L].
            d = sqrt( (source_xyz(:,1).' - mic_xyz(:,1)).^2 ...
                    + (source_xyz(:,2).' - mic_xyz(:,2)).^2 ...
                    + (source_xyz(:,3).' - mic_xyz(:,3)).^2 );
            
            if contains(source_type, '-')
                ps_sign = -1i;
            else
                ps_sign = +1i;
            end

            g = exp(ps_sign .* k .* d) ./ (4 .* pi .* d);  % Free-field.
            p = g * (s .* w);  % Weighted measurement [Q by 1]
            return;

        case {'ms', 'mixedsource', 'ms+', 'mixedsource+', 'ms-', 'mixedsource-', ...
              'mw', 'mixedwave', 'mw+', 'mixedwave+', 'mw-', 'mixedwave-'}

            % Distance of source to receiver ||x-y|| [Q by L].
            d = sqrt( (source_xyz(:,1).' - mic_xyz(:,1)).^2 ...
                    + (source_xyz(:,2).' - mic_xyz(:,2)).^2 ...
                    + (source_xyz(:,3).' - mic_xyz(:,3)).^2 );
                
            % Radius of each source to origin, used for normalization.
            src_r = rssq(source_xyz, 2);  % [L by 1].

            if contains(source_type, '-')
                mw_sign = -1i;
            else
                mw_sign = +1i;
            end

            % Point source Green's function.
            g = exp(mw_sign .* k .* d) ./ (4 .* pi .* d);  % [Q,L]

            % mw_norm = r_q * e^(-i*k*r_q)
            mw_norm = src_r .* exp(-mw_sign .* k .* src_r);  % [L,1]
           
            p = g * (s .* w .* mw_norm);  % [Q,1] = [Q,L] [L,1]
            return;
            
        otherwise
            temp = dbstack;
            error('%s(): failed to determine source type, srcType:%s', ...
                 temp(1).name, source_type);
    end

end