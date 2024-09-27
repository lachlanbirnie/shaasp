function [betas] = sph_betas_point_source(N, k, source_rtp)
% Lachlan Birnie
% 27-Sept-2024

arguments
    N {mustBeScalarOrEmpty}
    k (:,1)
    source_rtp (:,3)
end

k = permute(k, [2,3,1]);  % [1,1,K]
jn_mat = shaasp.sph_jn(N, k, source_rtp(:,1));  % [Q,N,K]
ynm_mat = shaasp.sph_ynm(N, source_rtp(:,2), source_rtp(:,3));  % [Q,N]

% betas = [1] .* [1,1,K]    .*    [N,Q,K]    .* [N,Q]
betas = +1i .* k .* permute(jn_mat, [2,1,3]) .* ynm_mat';
betas = sum(betas, 2);  % [N, 1, K]

end