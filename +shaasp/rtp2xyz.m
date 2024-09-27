function [x, y, z] = rtp2xyz(r, t, p)
% RTP2XYZ - Spherical coordinates to cartesian coordinates.
%   (radius, theta / elevation, phi / azimuth) -> (x, y, z)
%
% Syntax:  [x, y, z] = rtp2xyz(r, t, p)
%
% Inputs:
%
%   rtp     : [Q by 3] (r,t,p) position of each point (Q points total).
%             (radius, elevation 0 to pi, azimuth 0 to 2pi).
%   - or -    (output type depends on number of outputs).
%   r,t,p   : [Q by 1] vecotr of each points r, t, and p position.
%
% Outputs:
%
%   xyz     : [Q by 3] (x,y,z) position of each point (Q points total).
%   - or -    (input type depends on number of inputs).
%   x,y,z   : [Q by 1] vector of each points x,y,z position.
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
% Creation: 27-Jan-2021
% Last revision: 19-July-2024


    if (nargin == 1)
        % Split matrix input into vectors.
        temp = r;
        r = temp(:,1);
        t = temp(:,2);
        p = temp(:,3);
    end
    
    % Make sure positions are columns.
    if size(r,2) > size(r,1), r = r.'; end
    if size(t,2) > size(t,1), t = t.'; end
    if size(p,2) > size(p,1), p = p.'; end
    
    x = r .* sin(t) .* cos(p);
    y = r .* sin(t) .* sin(p);
    z = r .* cos(t);
    
    if (nargout == 1)
        % Combine outputs into matrix.
        x = [x, y, z];
    end
    
end