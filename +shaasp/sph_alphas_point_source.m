function [alphas, alphas_ell] = sph_alphas_point_source(N, k, source_rtp, kind)
% kind = 'first' / '+' / 'outgoing' (default) 
%     or 'second' / '-' / 'incoming'
% Lachlan Birnie
% 10-Jan-2025

arguments
    N {mustBeScalarOrEmpty}
    k (1,1,:) 
    source_rtp (:,3)
    kind {mustBeMember(kind, ["first", "outgoing", "+", "second", "incoming", "-"])} = '+'
end

if contains(kind, {'first', 'outgoing', '+'})
    coe_sign = +1i;
    hn_mat = shaasp.sph_hn(N, k, source_rtp(:,1));  % [N,Q,K]
elseif contains(kind, {'second', 'incoming', '-'})
    coe_sign = -1i;
    hn_mat = shaasp.sph_hn2(N, k, source_rtp(:,1));
end

ynm_mat = shaasp.sph_ynm(N, source_rtp(:,2), source_rtp(:,3));  % [N,Q]

% alphas = [1] .* [1,1,K] .* [N,Q,K] .* [N,Q]
alphas = coe_sign .* k .* hn_mat .* conj(ynm_mat);

% Coefficients of each source.
if nargout > 1
    alphas_ell = alphas;  % [N,L,K]
end

% Total sound field coefficients.
alphas = sum(alphas, 2);  % [N, 1, K]

end