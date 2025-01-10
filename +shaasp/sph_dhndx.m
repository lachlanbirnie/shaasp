function [dhdx] = sph_dhndx(N,k,r,options)
% SPH_DHNDX - Differential of spherical Hankel function: d/dx h_n(x)
%   Differential of first kind spherical Hankel function.
%   d/dx h_n(x) where x = k * r.
%
% Syntax:  [dhdx] = sph_dhndx(N,k,r)
%
% Inputs:
%
%   N   | Order of the differential spherical Hankel function matrix.
%   k   | [K,1] vector of wave numbers (frequency) arguments.
%   r   | [Q,1] vector of radius arguments (m).
%
% Outputs:
%
%   dhdx    : [(N+1)^2 by Q by K] matrix of d/dx h_n(x) terms.
%
% Equations:
%
%   First kind spherical Hankel function:
%
%       h_n(x) = sqrt(pi/2) * 1/sqrt(x) * H_(n+1/2)(x)
%
%       where H_n(x) is the first kind Hankel Function,
%       given by MATLAB inbuilt function 'besselh()'
%
%   Differential first kind spherical Hankel function:
%
%       d/dx h_n(x) = ( (n/x) * h_(n)(x) ) - h_(n+1)(x)
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
% Creation: 28-Jan-2021
% Last revision: 10-Jan-2025

    arguments
        N (1,1) {mustBeNonnegative, mustBeInteger}
        k (1,1,:) {mustBeNonnegative}
        r (1,:) {mustBeNonnegative}
        options.orientation {mustBeMember(options.orientation, ["[N,Q]", "[Q,N]", "[N,Q,K]", "[Q,N,K]"])} = '[N,Q]'
    end
    
    % Implicit inputs.
    K = length(k);
    Q = length(r);

    % Value of n for all nm pairs [0 : (N+1)^2].
    vec_n = repelem((0:N), 2.*(0:N)+1).';  %[N,1]

    % Expand function arguments into matrices [N,Q,K].
    arg_n = repmat(vec_n, [1, Q, K]);
    arg_n_plus1 = arg_n + 1;
    arg_k = repmat(k, [(N+1)^2, Q, 1]);
    arg_r = repmat(r, [(N+1)^2, 1, K]);
    
    % Spherical coefficient.
    sph_coe = sqrt(pi/2) .* (1 ./ sqrt(arg_k .* arg_r));
        
    % h_(n)(kr)
    hn = sph_coe .* besselh(arg_n + 0.5, arg_k .* arg_r);
    
    % h_(n+1)(kr)
    hn_plus1 = sph_coe .* besselh(arg_n_plus1 + 0.5, arg_k .* arg_r);
    
    % Differential first kind spherical Hankel function.
    dhdx = ( (arg_n ./ (arg_k .* arg_r)) .* hn) - hn_plus1;

    % Options orientation.
    switch options.orientation
        case {'[Q,N]', '[Q,N,K]'}
            dhdx = permute(dhdx, [2,1,3]);
        otherwise
    end

end