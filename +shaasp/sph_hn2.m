function [hn2] = sph_hn2(N,k,r)
% Wrapper function for second kind spherical Hankel function.
    import shaasp.sph_hn
    hn2 = conj(sph_hn(N, k, r));
end