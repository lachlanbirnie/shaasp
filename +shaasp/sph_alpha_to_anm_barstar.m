function [anm] = sph_alpha_to_anm_barstar(alpha)
% SPH_ALPHA_TO_ANM - Convert spherical coefficients to PWD coefficients
% Convert the standard spherical harmonic coefficients of alpha_nm(k)
% to the planewave decomposition coefficients of ((a_nm(k))^*)^*.
%
% Syntax:  [anm] = sph_alpha_to_anm_barstar(alpha)
%
% Inputs:
%    alpha - [(N+1)^2, K] set of sph coefficients
%
% Outputs:
%    anm - [(N+1)^2, K] set of pwd coefficients
%
% Example: 
%    alpha = rand(25,1) + 1i .* rand(25,);
%    [anm] = alpha_to_anm(alpha);
%
% Other m-files required: none
% Subfunctions: [a_n_negm] = nm_to_nminusm(a_n_m)
% MAT-files required: none
%
% See also: none
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 2021
% Last revision: 27-June-2023

    N = sqrt(size(alpha, 1)) - 1;
    n = repelem((0:N).', 2.*(0:N)+1);
    m = (1:(N+1)^2).' - n.^2 - n - 1;
    anm = (1./(4*pi)) .* (1./(+1i).^n) .* (-1).^(-m) .* nm_to_nminusm(alpha);

end


function [a_n_negm] = nm_to_nminusm(a_n_m) 
% NM_TO_NMINUSM - Reorder sph coefficients from a_nm to a_n(-m)
% Convert the standard spherical harmonic coefficients of alpha_nm(k)
%
% Syntax:  [a_n_negm] = nm_to_nminusm(a_n_m) 
%
% Inputs:
%    a_n_m - [(N+1)^2, K] set of sph coefficients
%
% Outputs:
%    a_n_negm - [(N+1)^2, K] set of sph coefficients
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: 2021
% Last revision: 27-June-2023

    N = sqrt(size(a_n_m, 1)) - 1;
    a_n_negm = zeros(size(a_n_m));
    for n = (0 : N)
        ind = ((n^2+n-n+1) : (n^2+n+n+1));
        a_n_negm(ind, :, :) = flipud(a_n_m(ind, :, :));
    end

end