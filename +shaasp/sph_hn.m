function [hn] = sph_hn(N,k,r)
% SPH_HN - Spherical Hankel function of the first kind.
%   Returns the first kind spherial Hankel function,
%   for input argument vectors (k,r) in matrix from, 
%   for all harmonics up to N ([N+1]^2 total).
%
% Syntax:  [hn] = sph_hn(N,k,r)
%
% Inputs:
%
%   N   Order of the returned spherial Hankel matrix.
%   k   [1 by K] vector of frequency (wave number) arguments.
%   r   [Q by 1] vector of radius arguments (m).
%
% Outputs:
%
%   hn  [Q by (N+1)^2 by K] matrix of the spherial Hankel function,
%       where h(a,b,c) = h_[n(b)](k(c)*r(a));
%
%       hn(:,:,k_1) = [h_0(k_1*r_1) h_1(k_1*r_1) ... h_N(k_1*r_1)
%                      h_0(k_1*r_2) h_1(k_1*r_2) ... h_N(k_1*r_2)
%                      ...
%                      h_0(k_1*r_Q) h_1(k_1*r_Q) ... h_N(k_1*r_Q)]
%
%       hn(:,:,k_2) = [h_0(k_2*r_1) h_1(k_2*r_1) ... h_N(k_2*r_1)
%                      h_0(k_2*r_2) h_1(k_2*r_2) ... h_N(k_2*r_2)
%                      ...
%                      h_0(k_2*r_Q) h_1(k_2*r_Q) ... h_N(k_2*r_Q)]
%
%   Note that the n'th order of each column increments as:
%       [0,1,1,1,2,2,2,2,2,3,3, ... N],
%   to match with the order-mode index's of (n,m) = (0,0) (1,-1) (1,0) ...
%
% Equaition:
%
%   h_n(x) = sqrt(pi/2) * 1/sqrt(x) * H_(n+1/2)(x)
%
%   where H_n(x) is the first kind Hankel Function,
%   given by MATLAB inbuilt function 'besselh()'
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
    validateattributes(r, {'double'},{'vector','ncols',1,'>',0});
                                
    % Implicit inputs.
    K = length(k);
    Q = length(r);
    
    % Value of n for all nm pairs [0 : (N+1)^2].
    v_n = @(N) repelem((0:N), 2.*(0:N)+1);

    % Expand arguments to matrix.
    arg_n = repmat(v_n(N), [Q, 1, K]);               % [Q,N,K]
    arg_k = repmat(reshape(k, [1,1,K]), [Q, (N+1)^2, 1]);   % [Q,N,K]
    arg_r = repmat(r, [1, (N+1)^2, K]);                     % [Q,N,K]
    
    % Solve first kind spherical Hankel function.
    hn = sqrt(pi/2) ...
     .* (1 ./ sqrt(arg_k .* arg_r)) ...
     .* besselh(arg_n + 0.5, arg_k .* arg_r);

end