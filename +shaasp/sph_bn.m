function [bn] = sph_bn(N,k,r,r_baff)
% SPH_BN - Spherical Rigid Baffle Equation, bn(kr).
%   Function assumes a rigid sphere, where the sphere radius is equal to
%   the receiver radii. 
%   Use sph_jn() for a open sphere case. 
%
% Syntax:  [bn] = sph_bn(N,k,r,R)
%
% Inputs:
%
%   N       Order of the returned baffle term matrix.
%   k       [1 by K] vector of frequency (wave number) arguments.
%   r       [Q by 1] vector of radius arguments (m).
%   r_baff  Radius of rigid baffle, default is same as r.
%
% Outputs:
%
%   b  [Q by (N+1)^2 by K] matrix of rigid baffle terms,
%      where b(a,b,c) = b_[n(b)](k(c)*r(a));
%
%       b(:,:,k_1) = [ b_0(k_1*r_1) b_1(k_1*r_1) ... b_N(k_1*r_1)
%                      b_0(k_1*r_2) b_1(k_1*r_2) ... b_N(k_1*r_2)
%                      ...
%                      b_0(k_1*r_Q) b_1(k_1*r_Q) ... b_N(k_1*r_Q) ]
%
%       b(:,:,k_2) = [ b_0(k_2*r_1) b_1(k_2*r_1) ... b_N(k_2*r_1)
%                      b_0(k_2*r_2) b_1(k_2*r_2) ... b_N(k_2*r_2)
%                      ...
%                      b_0(k_2*r_Q) b_1(k_2*r_Q) ... b_N(k_2*r_Q) ]
%
%   Note that the n'th order of each column increments as:
%       [0,1,1,1,2,2,2,2,2,3,3, ... N],
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
%       b_n(kr) = j_n(kr)                                   [open sphere]
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
% Last revision: 07-Jan-2025

    arguments
        N (1,1) {mustBeInteger}
        k
        r
        r_baff = r;
    end
    
    import shaasp.sph_jn
    import shaasp.sph_djndx
    import shaasp.sph_dhndx
    import shaasp.sph_hn
    
    bn = sph_jn(N,k,r) - (sph_djndx(N,k,r_baff) ./ sph_dhndx(N,k,r_baff)) .* sph_hn(N,k,r);

end