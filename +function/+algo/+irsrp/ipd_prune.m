function W = ipd_prune(W, X, f, tau_row, l_idx, q_idx, dTH, delta_step)
% Function  : ipd_prune
% Creation  : 24-07-25
% Author    : Weiting Lai
% Version   : -
% Description:
%   Inter-Phase Difference (IPD) based time-frequency bin removal
%   Based on: Dang et al., ICCAIS 2023, eq.(18)‑(23)
% -------------------------------------------------------------
% Inputs:
%   W         : nBin × nFrm × nPair   logical mask (pair-wise mask)
%   X         : nBin × nFrm × M       complex STFT of microphone signals
%   f         : nBin × 1              frequency vector (Hz)
%   tau_row   : 1 × nPair             TDOA values of the localized source for each mic pair
%   l_idx     : 1 × nPair             indices of first microphone in each pair
%   q_idx     : 1 × nPair             indices of second microphone in each pair
%   dTH       : scalar (radians)      IPD distance threshold (e.g., 0.3 rad)
%   delta_step: scalar (optional)     step size for δ grid search (default: 0.05 rad)
%
% Output:
%   W(:,:,p) = W(:,:,p) & ~bad_p   % for each pair p, mask out inconsistent bins

if nargin < 8, delta_step = 0.05; end
deltas = -pi/3 : delta_step :  pi/3;   % grid of phase offsets to test

[nBin,nFrm,~] = size(X);
nPair         = numel(l_idx);

for p = 1:nPair
    l = l_idx(p);  q = q_idx(p);

    % ---------- 1) Predicted phase difference ----------
    phi_pred = 2*pi * f(:) * tau_row(p);  % nBin × 1, expected phase difference

    % ---------- 2) Observed inter-channel phase difference (IPD) ----------
    ipd = angle( X(:,:,q) .* conj( X(:,:,l) ) );   % nBin × nFrm

    best_cnt = -inf;                     % initialize best score
    best_bad = false(nBin,nFrm);        % initialize best "bad" bin mask

    % ---------- 3) Brute-force search for optimal phase offset δ ----------
    for d = deltas
        dphi1 = angle(exp(1j*(ipd - (phi_pred + d))));
        dphi2 = angle(exp(1j*(ipd - (phi_pred - d))));
        bad   = min(abs(dphi1),abs(dphi2)) < dTH;   % Eq. (20)
        cnt   = nnz(bad & W(:,:,p));                % count only currently‑active bins
        if cnt > best_cnt
            best_cnt = cnt;
            best_bad = bad;
        end
    end
    W(:,:,p) = W(:,:,p) & ~best_bad;                % Eq. (23)
end
end