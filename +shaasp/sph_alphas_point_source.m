function [alphas] = sph_alphas_point_source(N, k, source_rtp, kind)
% kind = 'first' / '+' / 'outgoing' (default) 
%     or 'second' / '-' / 'incoming'
% Lachlan Birnie
% 27-Sept-2024

arguments
    N {mustBeScalarOrEmpty}
    k 
    source_rtp (:,3)
    kind {mustBeMember(kind, ["first", "outgoing", "+", "second", "incoming", "-"])} = '+'
end

if contains(kind, {'first', 'outgoing', '+'})
    coe_sign = +1i;
    hn_mat = shaasp.sph_hn(N, k, source_rtp(:,1));  % [Q,N,K]
elseif contains(kind, {'second', 'incoming', '-'})
    coe_sign = -1i;
    hn_mat = shaasp.sph_hn2(N, k, source_rtp(:,1));
end

k = permute(k, [2,3,1]);  % [1,1,K]
ynm_mat = shaasp.sph_ynm(N, source_rtp(:,2), source_rtp(:,3));  % [Q,N]

% alphas = [1] .* [1,1,K]    .*    [N,Q,K]         .* [N,Q]
alphas = coe_sign .* k .* permute(hn_mat, [2,1,3]) .* ynm_mat';
alphas = sum(alphas, 2);  % [N, 1, K]

end