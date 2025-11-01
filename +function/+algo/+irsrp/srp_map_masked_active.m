function srp_sum = srp_map_masked_active(X, f, tauTbl, l_idx, q_idx, W)
% P_all : nBin × nFrm × nPair   
% f     : nBin × 1              (Hz)
% tauTbl: nGrid × nPair
% W     : nBin × nFrm × nPair   (TF mask)
%
% Return：
% srp_sum: nGrid × 1

[nBin,~,~] = size(X);
nGrid = size(tauTbl,1);

srp_sum = zeros(nGrid,1,'single');

% Precompute exponentials for all pairs/grids (same as before)
E = exp( 2j*pi*f(:) .* permute(tauTbl,[3 2 1]) ); % nBin x nPair x nGrid

% for p = 1:numel(l_idx)
%     % skip if this pair has no active TF bins
%     Wp = W(:,:,p);                         % nBin x nFrm logical
%     if ~any(Wp(:)), continue; end
% 
%     % PHAT cross-spectrum for this pair
%     P = X(:,:,l_idx(p)) .* conj( X(:,:,q_idx(p)) );
%     P = P ./ (abs(P) + eps);
% 
%     % === key trick ===
%     % only accumulate bins that survive the mask, removing the frame dimension now
%     % S_freq(k) = sum_l P(k,l) for l where Wp(k,l)=1
%     P = P .* Wp;                           % zero-out pruned TF bins
%     S_freq = sum(P, 2);                    % nBin x 1 complex
% 
%     % keep only frequency rows that still have energy after masking
%     active_k = (S_freq ~= 0);
%     if ~any(active_k), continue; end
% 
%     % restrict exponentials to active frequencies only
%     Ep = reshape(E(:,p,:), nBin, []);      % nBin x nGrid
%     Ep = Ep(active_k, :);                  % K_active x nGrid
%     s  = S_freq(active_k).';               % 1 x K_active
% 
%     % complex GEMM is fast; real/imag split is optional
%     srp_sum = srp_sum + real(s * Ep).';    % nGrid x 1
% end

for p = 1:numel(l_idx)
    Wp = W(:,:,p);
    if ~any(Wp(:)), continue; end

    % PHAT cross-spectrum
    P = X(:,:,l_idx(p)) .* conj( X(:,:,q_idx(p)) );
    P = P ./ (abs(P) + eps);
    P = P .* Wp;                 % 掛遮罩
    S_freq = sum(P, 2);          % 沿時間加總 → nBin x 1

    active_k = (S_freq ~= 0);
    if ~any(active_k), continue; end

    % ←← Lazy exponentials：只對 active_k × N 現算 e^{j2π f τ}
    tau_p = tauTbl(:,p);         % N x 1
    Ep = exp( 2j*pi * f(active_k) .* tau_p.' );  % |K_active| x N
    s  = S_freq(active_k).';     % 1 x |K_active|

    srp_sum = srp_sum + real(s * Ep).';   % N x 1
end
end