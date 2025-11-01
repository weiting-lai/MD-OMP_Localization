function [alpha, pos, w_s] = md_omp(p, G, s)
% Function  : md_omp
% Creation  : 01-11-25
% Author    : Weiting Lai
% Version   : v1.0 (01-11-25)
%
% Description
%
%   Multi-Dictionary Orthogonal Matching Pursuit (MD-OMP) for near-field
%   source position estimation using multi-frame, multi-frequency data.
%   This implementation vectorises the joint-correlation across all time
%   frames and frequency bins, uses QR backslash (mldivide) for the
%   per-bin least-squares updates, and masks selected atoms to avoid
%   reselection. It returns the selected grid indices and the estimated
%   source weights in the same tensor layout as the input measurements.
%
% Inputs
%
%   p    : [M × T × J] complex measurements (microphone STFT slices),
%          with M microphones, T time frames, J frequency bins.
%   G    : [M × N × J] dictionary (per-bin Green’s function matrix);
%          G(:,n,f) is the atom from candidate grid index n at bin f.
%   s    : [scalar] sparsity level (number of sources to select).
%
% Outputs
%
%   alpha : [N × T × J] sparse coefficient tensor; rows at indices 'pos'
%           contain the per-bin LS solutions, other rows are zero.
%   pos   : [s × 1] selected candidate grid indices (1-based).
%   w_s   : [s × T × J] stacked per-iteration solutions (for convenience);
%           by construction, alpha(pos,:,:) == w_s and alpha(~pos,:,:) == 0.
%
% Algorithm (notation as in the paper)
%
%   Joint correlation (group-atom selection):
%       z_n = sum_{f=1..J} sum_{t=1..T} | ⟨ R_{:,t,f}, G_{:,n,f} ⟩ | .
%       λ = arg max_n z_n .
%
%   Support augmentation:
%       Append n = λ to the support; mask G(:,λ,:) to prevent reselection.
%
%   Group orthogonal projection (per frequency bin, all frames at once):
%       Let D_f = G(:,support,f). Solve D_f * W_f = P_f  via backslash,
%       where P_f = p(:,:,f). Update residual R = P − D*W.
%
%   Complexity (per iteration):
%       O(M·N·T·J) for correlations  +  O(M·k·T·J) for LS (k = iter).
%       Vectorization via pagemtimes and QR backslash accelerates both.
%
% Reference
%
%   W.-T. Lai et al., "Sound Source Localization using Multi-Dictionary
%   Orthogonal Matching Pursuit in Reverberant Environments", Sections II–IV,
%   esp. joint-correlation (eq. 13), selection (eq. 14), and projection (eq. 15).

M = size(p,1);
N = size(G,2);
T = size(p,2);
J = size(p,3);

pos   = zeros(s,1);
alpha = zeros(N,T,J, 'like', G);
w_s   = zeros(s,T,J, 'like', G);

dict  = zeros(M,s,J, 'like', G);   % selected atoms per bin
r_s   = p;                          % residual

% Precompute conjugate-transpose of G for all bins once:
% G_t: N x M x J (so pagemtimes(G_t, r_s) -> N x T x J)
G_t = conj(pagetranspose(G));

% selection mask to avoid re-choosing the same index
selected = false(N,1);

for k = 1:s
    % ---- 1) Correlation for all frames & bins in one shot (no for-t) ----
    % corr_full: N x T x J
    corr_full = abs(pagemtimes(G_t, r_s));       % G^H r_s
    % score: N x 1 (sum over frames and bins)
    score = sum(corr_full, [2 3]);               % requires R2018b+; else use sum(sum(corr_full,3),2)
    score = reshape(score, [N,1]);
    score(selected) = -inf;                      % mask out previously chosen atoms

    [~, idx] = max(score);                       % pick best atom (shared across bins)
    pos(k) = idx;
    selected(idx) = true;

    % ---- 2) Append selected atom to per-bin dictionaries ----
    dict(:,k,:) = G(:,idx,:);

    % ---- 3) LS update for ALL frames at once, per bin (avoid pinv) ----
    % Solve: dict(:,1:k,j) * w_s(1:k,:,j) = p(:,:,j)
    % Use QR via backslash; multiple RHS (T) solved in one call
    for j = 1:J
        Dj = dict(:,1:k,j);          % M x k
        w_s(1:k,:,j) = Dj \ p(:,:,j);% k x T
    end

    % ---- 4) Residual update (batched GEMM) ----
    r_s = p - pagemtimes(dict(:,1:k,:), w_s(1:k,:,:));
end

alpha(pos,:,:) = w_s;
end
