function [djndx] = sph_djndx(N,k,r)
% SPH_DJNDX - Differential of spherical Bessel function: d/dx j_n(x)
%   Differential of first kind spherical Bessel function
%   d/dx j_n(x) where x = k*r.
%
% Syntax:  [djndx] = sph_djndx(N,k,r)
%
% Inputs:
%
%   N   | Order of the differential spherical Bessel function matrix.
%   k   | [1 by K] vector of wave numbers (frequency) arguments.
%   r   | [Q by 1] vector of radius arguments (m).
%
% Outputs:
%
%   djdx    : [Q by (N+1)^2 by K] matrix of d/dx j_n(x) terms.
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
        
    % j_(n)(kr)
    jn = sph_coe .* besselj(arg_n + 0.5, arg_k .* arg_r);
    
    % j_(n+1)(kr)
    jn_plus1 = sph_coe .* besselj(arg_n_plus1 + 0.5, arg_k .* arg_r);
    
    % Address Bessel zeros.
    ind_x_zero = ((arg_k .* arg_r) == 0);
    if ~isempty(nonzeros(ind_x_zero))
        % j_(n)(kr)
        ind_n_zero = (arg_n == 0);
        jn(ind_x_zero & ind_n_zero) = 1;   % j_(n=0)(x=0) = 1
        jn(ind_x_zero & ~ind_n_zero) = 0;  % j_(n!=0)(x=0) = 0
        
        % j_(n+1)(kr)
        ind_np1_zero = (arg_n_plus1 == 0);
        jn_plus1(ind_x_zero & ind_np1_zero) = 1;   % j_(n=0)(x=0) = 1
        jn_plus1(ind_x_zero & ~ind_np1_zero) = 0;  % j_(n!=0)(x=0) = 0
    end
    
    % Differential first kind spherical Bessel function.
    djndx = ( (arg_n ./ (arg_k .* arg_r)) .* jn) - jn_plus1;

end