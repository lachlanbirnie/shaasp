function [bn] = sph_bn_cardioid(N,k,r,options)
% SPH_BN_CARDIOID - Array weighting function for open cardioid array.
%   Function assumes a open cardioid sphere, 
%   where the sphere radius is equal to the receiver radii. 
%
% Syntax:  [bn] = sph_bn_cardioid(N,k,r)
%
% Inputs:
%
%   N   Order of the returned baffle term matrix.
%   k   [1 by K] vector of frequency (wave number) arguments.
%   r   [Q by 1] vector of radius arguments (m).
%
% Outputs:
%
%   bn  [(N+1)^2 by Q by K] matrix of rigid baffle equation.
%
%       bn(:,:,k) = [ b_0(k * r1) ... b_0(k * rQ) ]
%                   [    ...             ...      ]
%                   [ b_N(k * r1) ... b_N(k * rQ) ]
%
% Equation:
%
%   Array Baffle Equation:
%
%       b_n(kr) = j_n(kr) - i * d/dx j_n(kr)
%
%                           /  d/dx j_n(kr)           \
%       b_n(kr) = j_n(kr) - |  ------------ X h_n(kr) |     [rigid sphere]
%                           \  d/dx h_n(kr)           /
%
%       b_n(kr) = j_n(kr)                                   [open sphere]
%
%
% Other m-files required: +shaasp. sph_jn, sph_djndx
% Subfunctions: none
% MAT-files required: none
%
% See also: sph_bn(),  sph_bn_rigid()
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 13-Dec-2024
% Last revision: 10-Jan-2025

    arguments
        N (1,1) {mustBeNonnegative, mustBeInteger}
        k (1,1,:) {mustBeNonnegative}
        r (1,:) {mustBeNonnegative}
        options.orientation {mustBeMember(options.orientation, ["[N,Q]", "[Q,N]", "[N,Q,K]", "[Q,N,K]"])} = '[N,Q]'
    end

    import shaasp.sph_jn
    import shaasp.sph_djndx
    
    % Cardioid baffle equation.
    bn = sph_jn(N,k,r) - 1i .* sph_djndx(N,k,r);

    % Options orientation.
    switch options.orientation
        case {'[Q,N]', '[Q,N,K]'}
            bn = permute(bn, [2,1,3]);
        otherwise
    end

end