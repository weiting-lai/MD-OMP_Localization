function [alpha, pos, w_s] = wb_omp(p, G, s)
% Function  : wb_omp
% Creation  : 01-11-25
% Author    : Weiting Lai
% Version   : v1.0 (01-11-25)
%
% Description
%
%   Wideband Orthogonal Matching Pursuit (WB-OMP) for near-field source
%   localization. This variant follows the paper’s WB-OMP: it first
%   averages the STFT snapshots across time, reducing the problem to a
%   frequency-only (wideband) selection, then performs group selection
%   across all frequency bins and per-bin least-squares projection using
%   QR backslash. It is a simplified variant of MD-OMP obtained by
%   time-averaging.
%
% Inputs
%
%   p    : [M × T × J] complex measurements (microphone STFT slices).
%          If T > 1, time snapshots are averaged:  p_use = mean(p,2)
%          (resulting in M × 1 × J).
%   G    : [M × N × J] dictionary (per-bin Green’s function matrix);
%          G(:,n,f) is the atom for candidate n at bin f.
%   s    : [scalar] sparsity level (number of sources to select).
%
% Outputs
%
%   alpha : [N × 1 × J] sparse coefficient tensor; rows at indices 'pos'
%           contain the per-bin LS solutions, others are zero.
%   pos   : [s × 1] selected candidate grid indices (1-based).
%   w_s   : [s × 1 × J] stacked solutions on the active support
%           (so that alpha(pos,:,:) == w_s).
%
% Algorithm (WB-OMP as MD-OMP with time averaging)
%
%   1) Time-average:  p_use = mean_t p  →  M × 1 × J.
%   2) Group correlation across bins:
%         z_n = ∑_{f=1..J} | ⟨ r_f , g_{n,f} ⟩ | ,
%      pick λ = arg max_n z_n; mask λ to avoid reselection.
%   3) Per-bin projection (all bins, single RHS):
%         D_f = G(:,support,f),  solve  D_f w_f = p_use(:,1,f)  via “\”.
%   4) Residual update (all bins):  r = p_use − D w .
%   5) Repeat 2–4 for s iterations.
%
% Complexity (per iteration)
%
%   O(M·N·J) for correlations  +  O(M·k·J) for LS (k = iteration index),
%   accelerated via pagemtimes and QR backslash; no pinv used.
%
% Reference
%
%   W.-T. Lai et al., “Sound Source Localization using Multi-Dictionary
%   Orthogonal Matching Pursuit in Reverberant Environments,” wideband
%   OMP variant (time-averaged MD-OMP) as described in the simulation
%   section.

M = size(p,1);
N = size(G,2);
T = size(p,2);
J = size(G,3);

% --- Time-average to make it wideband (freq-only) ---
% Resulting p_use is M x 1 x J
if T > 1
    p_use = mean(p, 2);              % average snapshots across time
else
    p_use = p;                        % already M x 1 x J
end

pos   = zeros(s,1);
alpha = zeros(N,1,J, 'like', G);
w_s   = zeros(s,1,J, 'like', G);

dict  = zeros(M,s,J, 'like', G);     % selected atoms per bin
r_s   = p_use;                        % residual (freq-only path)

% Precompute conjugate-transpose of G for all bins once:
% G_t: N x M x J (so pagemtimes(G_t, r_s) -> N x 1 x J)
G_t = conj(pagetranspose(G));

% selection mask to avoid re-choosing the same index
selected = false(N,1);

for k = 1:s
    % ---- 1) Correlation across all bins in one shot ----
    % corr_full: N x 1 x J
    corr_full = abs(pagemtimes(G_t, r_s));        % G^H r_s
    % score: N x 1 (sum over time dim (size 1) and bins)
    score = sum(corr_full, [2 3]);                % R2018b+; else sum(sum(corr_full,3),2)
    score = reshape(score, [N,1]);
    score(selected) = -inf;                       % mask out previously chosen atoms

    [~, idx] = max(score);                        % pick best atom (shared across bins)
    pos(k) = idx;
    selected(idx) = true;

    % ---- 2) Append selected atom to per-bin dictionaries ----
    dict(:,k,:) = G(:,idx,:);

    % ---- 3) LS update for ALL bins (single RHS per bin since T=1 after avg) ----
    % Solve: dict(:,1:k,j) * w_s(1:k,1,j) = p_use(:,1,j)
    for j = 1:J
        Dj = dict(:,1:k,j);                        % M x k
        w_s(1:k,1,j) = Dj \ p_use(:,1,j);          % k x 1
    end

    % ---- 4) Residual update (batched GEMM) ----
    r_s = p_use - pagemtimes(dict(:,1:k,:), w_s(1:k,:,:));   % M x 1 x J
end

alpha(pos,:,:) = w_s;
end
