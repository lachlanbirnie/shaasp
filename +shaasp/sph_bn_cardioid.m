function [bn] = sph_bn_cardioid(N,k,r)
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
% See also: sph_jn(),  sph_bn()
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 13-Dec-2024
% Last revision: 13-Dec-2024

    import shaasp.sph_jn
    import shaasp.sph_djndx
    
    bn = sph_jn(N,k,r) - 1i .* sph_djndx(N,k,r);

end