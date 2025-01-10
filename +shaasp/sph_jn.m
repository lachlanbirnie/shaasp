function [jn] = sph_jn(N,k,r,options)
% SPH_JN - Spherical Bessel fuction of the first kind.
%   Returns the first kind spherial Bessel function,
%   for input argument vectors (k,r) in matrix form, 
%   for all harmonics up to N ([N+1]^2 total).
%
% Syntax:  [jn] = sph_jn(N,k,r)
%
% Inputs:
%
%   N   Order of the returned spherial Bessel matrix.
%   k   [K,1] vector of frequency (wave number) arguments.
%   r   [Q,1] vector of radius arguments (meters).
%
%   options
%       orientation - '[N,Q,K]' (default), '[Q,N,K]'
%
% Outputs:
%
%   jn  [(N+1)^2 by Q by K] matrix of the spherial Bessel function,
%       where J(a,b,c) = j_[n(a)](k(c)*r(b));
%
%       jn(:,:,k) = [ j_0(k * r1) ... j_0(k * rQ) ]
%                   [    ...             ...      ]
%                   [ j_N(k * r1) ... j_N(k * rQ) ]
%
%   Note that the n'th order of each row increments as:
%       [0,1,1,1,2,2,2,2,2,3,3, ... N].',
%   to match with the order-mode index's of (n,m) = (0,0) (1,-1) (1,0) ...
%
% Equation:
%            ____
%           | pi    1
%   jn(x) = | --   ___  J_(n+1 / 2)(x)
%           / 2   / x
%
%   where J_n(x) is the first kind Bessel Function,
%   given by MATLAB inbuilt function 'besselj()'
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
% Creation: 13-Feb-2019
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

    % Value of n for all nm pairs along the jn rows [(N+1)^2, 1].
    row_n = repelem((0:N).', 2 .* (0 : N) + 1);
    
    % Expand arguments to matrix of same size and dimension [N,Q,K].
    arg_n = repmat(row_n, [1, Q, K]);
    arg_k = repmat(k, [(N+1)^2, Q, 1]);
    arg_r = repmat(r, [(N+1)^2, 1, K]);
    
    % Solve Spherical Bessel Function.
    jn = sqrt(pi/2) .* (1 ./ sqrt(arg_k .* arg_r)) .* besselj(arg_n + 0.5, arg_k .* arg_r);
    
    % Address Bessel zeros.
    ind_x_is_zero = ((arg_k .* arg_r) == 0);
    if ~isempty(nonzeros(ind_x_is_zero))
        ind_n_is_zero = (arg_n == 0);
        jn(ind_x_is_zero & ind_n_is_zero) = 1;   % j_(n=0)(x=0) = 1
        jn(ind_x_is_zero & ~ind_n_is_zero) = 0;  % j_(n!=0)(x=0) = 0
    end

    % Options orientation for [jn].
    switch options.orientation
        case {'[Q,N]', '[Q,N,K]'}
            jn = permute(jn, [2,1,3]);
        otherwise
            % jn = permute(jn, [1,2,3]);
    end

end