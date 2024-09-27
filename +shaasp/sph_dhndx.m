function [dhdx] = sph_dhndx(N,k,r)
% SPH_DHNDX - Differential of spherical Hankel function: d/dx h_n(x)
%   Differential of first kind spherical Hankel function.
%   d/dx h_n(x) where x = k * r.
%
% Syntax:  [dhdx] = sph_dhndx(N,k,r)
%
% Inputs:
%
%   N   | Order of the differential spherical Hankel function matrix.
%   k   | [1 by K] vector of wave numbers (frequency) arguments.
%   r   | [Q by 1] vector of radius arguments (m).
%
% Outputs:
%
%   dhdx    : [Q by (N+1)^2 by K] matrix of d/dx h_n(x) terms.
%
% Equations:
%
%   First kind spherical Hankel function:
%
%       j_n(x) = sqrt(pi/2) * 1/sqrt(x) * H_(n+1/2)(x)
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

    % Expand function arguments into matrices [Q,N,K].
    arg_n = repmat(v_n(N), [Q, 1, K]);
    arg_n_plus1 = arg_n + 1;
    arg_k = repmat(reshape(k, [1,1,K]), [Q, (N+1)^2, 1]);
    arg_r = repmat(r, [1, (N+1)^2, K]);
    
    % Spherical coefficient.
    sph_coe = sqrt(pi/2) .* (1 ./ sqrt(arg_k .* arg_r));
        
    % h_(n)(kr)
    hn = sph_coe .* besselh(arg_n + 0.5, arg_k .* arg_r);
    
    % h_(n+1)(kr)
    hn_plus1 = sph_coe .* besselh(arg_n_plus1 + 0.5, arg_k .* arg_r);
    
    % Differential first kind spherical Hankel function.
    dhdx = ( (arg_n ./ (arg_k .* arg_r)) .* hn) - hn_plus1;

end