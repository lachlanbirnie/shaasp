function [djndx] = sph_djndx(N,k,r,options)
% SPH_DJNDX - Differential of spherical Bessel function: d/dx j_n(x)
%   Differential of first kind spherical Bessel function
%   d/dx j_n(x) where x = k*r.
%
% Syntax:  [djndx] = sph_djndx(N,k,r)
%
% Inputs:
%
%   N   | Order of the differential spherical Bessel function matrix.
%   k   | [K,1] vector of wave numbers (frequency) arguments.
%   r   | [Q,1] vector of radius arguments (m).
%
% Outputs:
%
%   djdx    : [(N+1)^2 by Q by K] matrix of d/dx j_n(x) terms.
%
% Equations:
%
%   First kind spherical Bessel function:
%
%       j_n(x) = sqrt(pi/2) * 1/sqrt(x) * J_(n+1/2)(x)
%
%       where J_n(x) is the first kind Bessel Function,
%       given by MATLAB inbuilt function 'besselj()'
%
%   Differential first kind spherical Bessel function:
%
%       d/dx j_n(x) = ( (n/x) * j_(n)(x) ) - j_(n+1)(x)
%
%   if x = 0:
%
%   d/dx jn(x) = 0      [for n != 1]
%   d/dx jn(x) = 1/3    [for n == 1]
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
    vec_n = repelem((0:N), 2.*(0:N)+1).';  % [N,1]
    
    % Expand function arguments into matrices [N,Q,K].
    arg_n = repmat(vec_n, [1, Q, K]);
    arg_n_plus1 = arg_n + 1;
    arg_k = repmat(k, [(N+1)^2, Q, 1]);
    arg_r = repmat(r, [(N+1)^2, 1, K]);
    
    % Spherical coefficient.
    sph_coe = sqrt(pi/2) .* (1 ./ sqrt(arg_k .* arg_r));
        
    % j_(n)(kr)
    jn = sph_coe .* besselj(arg_n + 0.5, arg_k .* arg_r);
    
    % j_(n+1)(kr)
    jn_plus1 = sph_coe .* besselj(arg_n_plus1 + 0.5, arg_k .* arg_r);
    
    % Address Bessel zeros.
    ind_x_is_zero = ((arg_k .* arg_r) == 0);
    if ~isempty(nonzeros(ind_x_is_zero))
        % j_(n)(kr)
        ind_n_is_zero = (arg_n == 0);
        jn(ind_x_is_zero & ind_n_is_zero) = 1;   % j_(n=0)(x=0) = 1
        jn(ind_x_is_zero & ~ind_n_is_zero) = 0;  % j_(n!=0)(x=0) = 0
        
        % j_(n+1)(kr)
        ind_np1_is_zero = (arg_n_plus1 == 0);
        jn_plus1(ind_x_is_zero & ind_np1_is_zero) = 1;   % j_(n=0)(x=0) = 1
        jn_plus1(ind_x_is_zero & ~ind_np1_is_zero) = 0;  % j_(n!=0)(x=0) = 0
    end
    
    % Differential first kind spherical Bessel function.
    djndx = ( (arg_n ./ (arg_k .* arg_r)) .* jn) - jn_plus1;

    % Address differential bessel zeros.
    if ~isempty(ind_x_is_zero)
        ind_n_is_one = (arg_n == 1);
        djndx(ind_x_is_zero & ind_n_is_one) = 1/3;  % d/dx jn(x) = 1/3    [for n == 1]
        djndx(ind_x_is_zero & ~ind_n_is_one) = 0;  % d/dx jn(x) = 0      [for n != 1]
    end

    % Options orientation.
    switch options.orientation
        case {'[Q,N]', '[Q,N,K]'}
            djndx = permute(djndx, [2,1,3]);
        otherwise
    end

end