function [bn] = sph_bn_rigid(N,k,r,r_baff,options)
% SPH_BN - Spherical Rigid Baffle Equation, bn(kr).
%   Function assumes a rigid sphere, where the sphere radius is r_baff.
%   Use sph_jn() for a open sphere case. 
%
% Syntax:  [bn] = sph_bn_rigid(N,k,r,R)
%
% Inputs:
%
%   N       [1,1] Order of the returned baffle term matrix.
%   k       [K,1] vector of frequency (wave number) arguments.
%   r       [Q,1] vector of radius arguments (m).
%   r_baff  Radius of rigid baffle, default is same as r.
%
% Outputs:
%
%   bn  [(N+1)^2 by Q by K] matrix of rigid baffle equation.
%
%       bn(:,:,k) = [ b_0(k * r1) ... b_0(k * rQ) ]
%                   [    ...             ...      ]
%                   [ b_N(k * r1) ... b_N(k * rQ) ]
%
%   Note that the n'th order of each column increments as:
%       [0,1,1,1,2,2,2,2,2,3,3, ... N].',
%   to match with the order-mode index's of (n,m) = (0,0) (1,-1) (1,0) ...
%
% Equation:
%
%   Rigid Baffle Equation:
%
%                           /  d/dx j_n(kR)           \
%       b_n(kr) = j_n(kr) - |  ------------ X h_n(kr) |     [rigid sphere]
%                           \  d/dx h_n(kR)           /
%
%
% Other m-files required: +shaasp. sph_jn, sph_hn, sph_djndx, sph_dhndx
% Subfunctions: none
% MAT-files required: none
%
% See also: sph_jn,  sph_bn_cardioid.
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 28-Jan-2021
% Last revision: 10-Jan-2025

    arguments
        N (1,1) {mustBeNonnegative, mustBeInteger}
        k (1,1,:) {mustBeNonnegative}
        r (1,:) {mustBeNonnegative}
        r_baff (1,:) {mustBeNonnegative} = r
        options.orientation {mustBeMember(options.orientation, ["[N,Q]", "[Q,N]", "[N,Q,K]", "[Q,N,K]"])} = '[N,Q]'
    end
    
    import shaasp.sph_jn
    import shaasp.sph_djndx
    import shaasp.sph_dhndx
    import shaasp.sph_hn
    
    % Rigid baffle equation.
    bn = sph_jn(N,k,r) - (sph_djndx(N,k,r_baff) ./ sph_dhndx(N,k,r_baff)) .* sph_hn(N,k,r);

    % Options orientation.
    switch options.orientation
        case {'[Q,N]', '[Q,N,K]'}
            bn = permute(bn, [2,1,3]);
        otherwise
    end

end