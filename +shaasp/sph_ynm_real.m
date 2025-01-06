function [Y] = sph_ynm_real(N,theta,phi)
% SPH_YNM_REAL - Real Spherical Harmonic Basis Functions Y_(n,m)(theta,phi)
%   Calculate real spherical harmonic basis function for each element of
%   input vectors (theta, phi) and each harmonic pair up to order N.
%
% Syntax:  [Y] = sph_ynm_real(N,theta,phi)
%
% Inputs:
%
%   N       Order of the returned spherical harmonic matrix.
%   t       [Q by 1] vector of theta positions (rad).
%   p       [Q by 1] vector of phi positions (rad).
%
% Outputs:
%
%   Y       [Q by (N+1)^2] matrix of real spherical harmonic equations,
%           where each element of Y(i,j) = Y(nm, t(i), p(i)).
%
%           Y = [ Y_00(t1,p1) Y_1-1(t1, p1) ... Y_N+N(t1,p1)
%                 Y_00(t2,p2) Y_1-1(t2, p2) ... Y_N+N(t2,p2)
%                 ...
%                 Y_00(tq,pq) Y_1-1(tq, pq) ... Y_N+N(tq,pq)
%                 ...
%                 Y_00(tQ,pQ) Y_1-1(tQ, pQ) ... Y_N+N(tQ,pQ) ]
%
% Equations:
%
%   Real Spherical Harmonic Equation ~Ynm(t,p):
%
%   [for m = 0]
%                     ________
%                    |(2n + 1)
%        ~Ynm(t,p) = | ------  Pnm(cos t)
%                    /  4pi
%
%   [for m > 0]
%                        ________  ________
%                     _ |(2n + 1) |(n - m)!
%        ~Ynm(t,p) = /2 | ------  | ------  Pnm(cos t) cos(m p)
%                       /  4pi    /(n + m)!
%   [for m < 0]
%                        ________  __________
%                     _ |(2n + 1) |(n - |m|)!
%        ~Ynm(t,p) = /2 | ------  | --------  Pn|m|(cos t) sin(|m| p)
%                       /  4pi    /(n + |m|)!
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: sph_ynm,  OTHER_FUNCTION_NAME2
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 06-Jan-2025
% Last revision: 06-Jan-2025
%
% Extra Equations:
%
%   Complex Spherical Harmonic Equation:
%                    __________________
%                   | (2n + 1) (n - m)!
%       Ynm(t,p) =  | -------- -------- Pnm(cos t) e^(i m p)
%                   /   4pi    (n + m)!
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
%   Note that MATLAB includes legendre the Condon-Shortley phase.

    % Input checking.
    arguments
        N (1,1) {mustBeInteger}
        theta (:,1) {mustBeNumeric}
        phi (:,1) {mustBeNumeric}
    end

    if length(theta) ~= length(phi)
        error('Theta and Phi must be equal lengths');
    end

    % Implicit inputs.
    Q = length(theta);  % Number of positions.
    
    % Associate Legendre Function and Exponential term. 
    L = zeros(Q, (N+1)^2);
    E = zeros(Q, (N+1)^2);
    for n = (0 : N)
        Ln = legendre(n, cos(theta));  % [(N+1),Q] for terms (n, m = 0->n).
        for m = (0 : n)
            % Positive and zero modes.
            L(:, n^2+n+m+1) = Ln(m+1,:);
            E(:, n^2+n+m+1) = cos(m * phi);  % == 1 for m is zero.
            if m ~= 0
                % Negative modes.
                L(:, n^2+n-m+1) = Ln(abs(m)+1,:);
                E(:, n^2+n-m+1) = sin(abs(m) * phi);
            end
        end % m
    end % n

    % (n,m) values for each column in ~Ynm == [Q by (N+1)^2].
    col_n = repelem((0:N), 2 * (0:N) + 1);
    col_m = (1:(N+1)^2) - col_n.^2 - col_n - 1;
    
    % Normalisation term.
    norm = sqrt((2 .* col_n + 1) ...
                ./ (4 .* pi)) ...
        .* sqrt(factorial(col_n - abs(col_m)) ...
             ./ factorial(col_n + abs(col_m)));

    % sqrt(2) term for non-zero m.
    r2 = ones(size(norm));
    r2(col_m ~= 0) = sqrt(2);
    
    % Real Spherical Harmonic Function.
    Y = r2 .* norm .* L .* E;  % [Q,(N+1)^2]
end