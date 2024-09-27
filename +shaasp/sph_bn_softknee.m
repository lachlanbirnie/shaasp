function [b_soft] = sph_bn_softknee(N,k,r,amp)
% Soft knee bn(kr) (without 1/4pi i^n) from B. Bernschutz.
% https://d-nb.info/1156013852/34
%
% amp = gain amplitude in dB (default 62 dB).

    import shaasp.sph_bn

    if (nargin < 4), amp = 62; end  % Default amp of 62 dB.
    a = 10^(amp/20);  % Gain constant.
    b = sph_bn(N,k,r);  % bn(kr) function.
    b_soft = pi./(2*a) .* b./abs(b) .* 1./acot((2*a)./pi .* abs(b));
    
end