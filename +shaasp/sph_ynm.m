function [Y] = sph_ynm(N,theta,phi)
% SPH_YNM - Spherical Harmonic Basis Functions Y_(n,m)(theta, phi)
%   Calculate the spherical harmonic basis function for each element of the
%   input vectors (theta, phi) and each harmonic pair up to order N.
%   Uing Fourier Acoustics book definiton.
%
% Syntax:  [Y] = sph_ynm(N,theta,phi)
%
% Inputs:
%
%   N       Order of the returned spherical harmonic matrix.
%   t       [Q by 1] vector of theta positions (rad).
%   p       [Q by 1] vector of phi positions (rad).
%
% Outputs:
%
%   Y       [Q by (N+1)^2] matrix of spherical harmonic equations,
%           where each element of Y(i,j) = Y(jth nm, t(i), p(i)).
%
%           Y = [ Y_00(t_1,p_1) Y_1-1(t_1, p_1) ... Y_N+N(t_1,p_1)
%                 Y_00(t_2,p_2) Y_1-1(t_2, p_2) ... Y_N+N(t_2,p_2)
%                 ...
%                 Y_00(t_q,p_q) Y_1-1(t_q, p_q) ... Y_N+N(t_q,p_q)
%                 ...
%                 Y_00(t_Q,p_Q) Y_1-1(t_Q, p_Q) ... Y_N+N(t_Q,p_Q) ]
%
% Equations:
%
%   Spherical Harmonic Equation:
%                    __________________
%                   | (2n + 1) (n - m)!
%       Ynm(t,p) =  | -------- -------- Pnm(cos t) e^(i m p)
%                  /     4pi   (n + m)!
%
%   Associate Legendre Functions:
%
%       Pnm(x) = (-1)^m (1 - x^2)^(m/2) (dm/dx^2) Pn(x)  [for m >= 0]
%                ______
%                  ^ note the Condon-Shortley phase.               
%
%                            (n - |m|)!
%       Pn(-m)(x) = (-1)^|m| ---------- Pn|m|(x)         [for m < = 0]
%                            (n + |m|)! 
%
%   Legendre Polynomials:
%   
%                 1      d^n
%       Pn(x) = ------ ------ (x^2 - 1)^n
%               2^n n!  dx^n
%
%   Note that the above Associate Legendre Function and the Legendre
%   Polynomials are the same definition as the MATLAB inbuilt
%   function: 'legendre()'.
%   Note that MATLAB includes the Condon-Shortley phase.
%
% Example: 
%   Line 1 of example
%   Line 2 of example
%   Line 3 of example
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
% Creation: 13-Feb-2024
% Last revision: 19-July-2024

    % Input checking.
    validateattributes(N, {'double'},{'integer','>=',0,'size',[1,1]});
    if isrow(theta), theta = theta.'; end
    validateattributes(theta, {'double'},{'vector','ncols',1});
    if isrow(phi), phi = phi.'; end
    validateattributes(phi, {'double'},{'vector','ncols',1,'size',size(theta)});

    % Implicit inputs.
    Q = length(theta);

    % Value of n for all nm pairs [0 : (N+1)^2].
    v_n = @(N) repelem((0:N), 2.*(0:N)+1);

    % Value of m for all nm pairs [0 : (N+1)^2].
    v_m = @(N) (1:(N+1)^2) - v_n(N).^2 - v_n(N) - 1;

    % Get set of harmonic (n,m)'s in vectors and matrices.
    vec_n = v_n(N);                   % [1,(N+1)^2]
    vec_m = v_m(N);                   % [1,(N+1)^2]
    arg_phi = repmat(phi, [1, (N+1)^2]);         % [Q,(N+1)^2]
    arg_m = repmat(v_m(N), [Q, 1]);   % [Q,(N+1)^2]
    
    % Solve Associate Legendre Function. 
    L = zeros(Q, (N+1)^2);

    for n = (0 : N)
        
        Ln = legendre(n, cos(theta));  % [(N+1), Q] for terms (n, m = 0>n).
        
        for m = (0 : n)
            
            % Positive modes.
            L(:, n^2+n+m+1) = Ln(m+1,:);
            
            % Negative modes.
            L(:, n^2+n-m+1) = (-1)^(m) .* (factorial(n-m)/factorial(n+m)) .* Ln(m+1,:);

        end % m
    end % n
    
    % Solve normalisation term.
    norm = sqrt((2.*vec_n + 1) ./ (4*pi)) ...  % [1,(N+1)^2]
        .* sqrt(factorial(vec_n - vec_m) ./ factorial(vec_n + vec_m));
    
    % Solve exponential term.
    E = exp(1i .* arg_m .* arg_phi);  % [Q,(N+1)^2]
    
    % Solve Spherical Harmonic Function.
    Y = norm .* L .* E;  % [Q,(N+1)^2]

end