function [betas, betas_ell] = sph_betas_point_source(N, k, source_rtp)
% Lachlan Birnie
% 10-Jan-2025

arguments
    N {mustBeScalarOrEmpty}
    k (1,1,:)
    source_rtp (:,3)
end

jn_mat = shaasp.sph_jn(N, k, source_rtp(:,1));  % [N,L,K]
ynm_mat = shaasp.sph_ynm(N, source_rtp(:,2), source_rtp(:,3));  % [N,L]

% betas = [1] .* [1,1,K] .* [N,Q,K] .* [N,Q]
betas = +1i .* k .* jn_mat .* conj(ynm_mat);

% Coefficients of each source.
if nargout > 1
    betas_ell = betas;  % [N,L,K]
end

% Total sound field coefficients.
betas = sum(betas, 2);  % [N, 1, K]

end