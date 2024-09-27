function [invb_soft] = sph_bninv_softknee(N,k,r,amp)
% Soft knee radial filter (without 1/4pi i^n) from B. Bernschutz.
% https://d-nb.info/1156013852/34
%
% amp = gain amplitude in dB (default 62 dB).

    import shaasp.sph_bn

    if (nargin < 4), amp = 62; end  % Default amp of 62 dB.
    a = 10^(amp/20);  % Gain constant.
    b = sph_bn(N,k,r);  % bn(kr) function.
    invb = 1 ./ b;
    invb_soft = (2*a)./pi .* invb./abs(invb) .* atan(pi./(2*a) .* abs(invb));
    
end