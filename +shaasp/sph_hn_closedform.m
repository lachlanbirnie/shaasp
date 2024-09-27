function [hn, hnm] = sph_hn_closedform(N,k,r)
% SPH_HN_CLOSEDFORM - Get the closed form expression for hn(kr).
%
% Syntax:  [hn, hnm] = sph_hn_closedform(N,k,r)
%
% Outputs:
%
%   hn = [N+1, k, r]
%   hnm = [(N+1)^2, k, r] (padded for nm pairs)
%
% See also: sph_hn, sph_hn2_closedform
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 03-02-2023
% Last revision: 03-02-2023


    kvec = reshape(k, [1, length(k)]);  % [1,k]
    rvec = reshape(r, [1, 1, length(r)]);  % [1,1,r]
    
    krvec = reshape(kvec .* rvec, [1, length(k)*length(r)]);  % [1,kr]
    
    ekr = exp(1i .* krvec) ./ krvec;  % [1,kr]
    
    hnkr = zeros(N+1, length(k)*length(r));  % [N+1, kr]
    
    for n = (0:N)
        m = (0:n).';
        c1 = (factorial(n+m) ./ factorial(n-m)).';  % [1,m]
        c2 = (1i).^m ./ (factorial(m) .* (2.*krvec).^m);  % [m,kr]
        hnkr(n+1,:) = (-1i).^(n+1) .* ekr .* (c1 * c2);  % [(n),kr]
    end
    
    hn = reshape(hnkr, [N+1, length(k), length(r)]);  % [N+1, k, r]
    
    if (nargout > 1)
        ind = repelem((0:N)+1, 2.*(0:N)+1);
        hnm = hn(ind, :, :);  % [(N+1)^2, k, r]
    end

end