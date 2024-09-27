function [jn] = sph_jn(N,k,r)
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
%   k   [1 by K] vector of frequency (wave number) arguments.
%   r   [Q by 1] vector of radius arguments (meters).
%
% Outputs:
%
%   jn  [Q by (N+1)^2 by K] matrix of the spherial Bessel function,
%       where J(a,b,c) = j_[n(b)](k(c)*r(a));
%
%       jn(:,:,k_1) = [j_0(k_1*r_1) j_1(k_1*r_1) ... j_N(k_1*r_1)
%                      j_0(k_1*r_2) j_1(k_1*r_2) ... j_N(k_1*r_2)
%                      ...
%                      j_0(k_1*r_Q) j_1(k_1*r_Q) ... j_N(k_1*r_Q)]
%
%       jn(:,:,k_2) = [j_0(k_2*r_1) j_1(k_2*r_1) ... j_N(k_2*r_1)
%                      j_0(k_2*r_2) j_1(k_2*r_2) ... j_N(k_2*r_2)
%                      ...
%                      j_0(k_2*r_Q) j_1(k_2*r_Q) ... j_N(k_2*r_Q)]
%
%   Note that the n'th order of each column increments as:
%       [0,1,1,1,2,2,2,2,2,3,3, ... N],
%   to match with the order-mode index's of (n,m) = (0,0) (1,-1) (1,0) ...
%
% Equation:
%
%   j_n(x) = sqrt(pi/2) * 1/sqrt(x) * J_(n+1/2)(x)
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
% Last revision: 19-July-2024


    % Input checking.
    validateattributes(N, {'double'},{'integer','>=',0,'scalar'});
    if iscolumn(k), k = k.'; end
    validateattributes(k, {'double'},{'vector','nrows',1,'>=',0});
    if isrow(r), r = r.'; end
    validateattributes(r, {'double'},{'vector','ncols',1,'>=',0});
                                
    % Implicit inputs.
    K = length(k);
    Q = length(r);

    % Value of n for all nm pairs [0 : (N+1)^2].
    v_n = @(N) repelem((0:N), 2.*(0:N)+1);
    
    % Expand arguments to matrix.
    arg_n = repmat(v_n(N), [Q, 1, K]);               % [Q,N,K]
    arg_k = repmat(reshape(k, [1,1,K]), [Q, (N+1)^2, 1]);   % [Q,N,K]   
    arg_r = repmat(r, [1, (N+1)^2, K]);                     % [Q,N,K]
    
    % Solve Spherical Bessel Function.
    jn = sqrt(pi/2) .* (1 ./ sqrt(arg_k .* arg_r)) .* besselj(arg_n + 0.5, arg_k .* arg_r);
    
    % Address Bessel zeros.
    ind_x_zero = ((arg_k .* arg_r) == 0);
    if ~isempty(nonzeros(ind_x_zero))
        ind_n_zero = (arg_n == 0);
        jn(ind_x_zero & ind_n_zero) = 1;   % j_(n=0)(x=0) = 1
        jn(ind_x_zero & ~ind_n_zero) = 0;  % j_(n!=0)(x=0) = 0
    end

end