function [Y] = sph_ynm(N,theta,phi,options)
% SPH_YNM - Spherical Harmonic Basis Functions Y_(n,m)(theta, phi, Options)
%   Calculate the spherical harmonic basis function for each element of the
%   input vectors (theta, phi) and each harmonic (n,m) pair up to order N.
%   Uing Fourier Acoustics book definiton.
%
% Syntax:  [Y] = sph_ynm(N,theta,phi)
%
% Inputs:
%
%   N       Order of the returned spherical harmonic matrix.
%   theta   [Q,1] elevation positions [0,pi] downwards from z-axis.
%   phi     [Q,1] azimuth positions [0, 2pi) counterclockwise from x-axis.
%
%   options
%       orientation = '[N,Q]' (default), '[Q,N]'.
%
% Outputs:
%
%   Y       [(N+1)^2 by Q] matrix of spherical harmonic basis functions,
%           where each element of Y(i,q) = Y(nm(i), theta(q), phi(q)).
%
%           Y = [ Y_00(theta[1], phi[1]) ... Y_00(theta[Q], phi[Q]) ]
%               [             ...                   ...             ]
%               [ Y_NN(theta[1], phi[1]) ... Y_NN(theta[Q], phi[Q]) ]
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
%   Y = shaasp.sph_ynm(4, pi/2, pi/3);
%   Y = shaasp.sph_ynm(4, pi/2, pi/3, 'orientation', '[Q,N]');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: sph_jn,  sph_ynm_real
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 13-Feb-2024
% Last revision: 10-Jan-2025

    arguments
        N (1,1) {mustBeNonnegative, mustBeInteger}
        theta (1,:) {mustBeNumeric}
        phi (1,:) {mustBeNumeric}
        options.orientation {mustBeMember(options.orientation, ["[N,Q]", "[Q,N]"])} = '[N,Q]'
    end

    if length(theta) ~= length(phi)
        error('Arguments theta and phi must be same length');
    end

    % Solve Associate Legendre Function. 
    L = zeros((N+1)^2, length(theta));  % [(N+1)^2, Q]
    for n = (0 : N)
        Ln = legendre(n, cos(theta));  % [(N+1), Q] for terms (n, m = 0>n).
        for m = (0 : n)
            % Positive and zero modes.
            L(n^2+n+m+1, :) = Ln(m+1, :);
            if (m ~= 0)
                % Negative modes.
                L(n^2+n-m+1, :) = (-1)^(m) .* (factorial(n-m)/factorial(n+m)) .* Ln(m+1, :);
            end
        end % m
    end % n

    % Create pairs of (n,m) along column vector.
    col_n = repelem((0:N).', 2.*(0:N)+1);  % [(N+1)^2, 1]
    col_m = (1:(N+1)^2).' - col_n.^2 - col_n - 1;  % [(N+1)^2, 1]
    
    % Solve normalisation term.
    norm = sqrt((2.*col_n + 1) ./ (4*pi)) ...  % [(N+1)^2, 1]
        .* sqrt(factorial(col_n - col_m) ./ factorial(col_n + col_m));
    
    % Solve exponential term.
    E = exp(1i .* col_m .* phi);  % [(N+1)^2, Q]
    
    % Solve Spherical Harmonic Function.
    Y = norm .* L .* E;  % [(N+1)^2, Q]

    % Options orientation for [Y].
    switch options.orientation
        case '[Q,N]'
            Y = Y.';
        otherwise
            % Y = Y;
    end

end