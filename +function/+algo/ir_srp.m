function [est_xyz, sum_bin] = ir_srp(X, mic_xyz, vs_xyz, f, s_amount, params)
% Function  : ir_srp
% Creation  : 24-07-25
% Author    : Weiting Lai
% Version   : -
% Description
%
%   Implement Iteratively Reweighted Steered Response Power with Phase 
%   Transform (IR-SRP-PHAT) (G-IRLS) from the reference paper.
%
% Reference  
%   Dang et al., An iteratively reweighted steered response power approach 
%   to multisource localization using a distributed microphone network,
%   J. Acoust. Soc. Amer., 2024.
% -------------------------------------------------------------
% IR-SRP
%   • screens bins with coherence & SNR masks
%   • loops: peak → IPD‑based bin removal → recompute map
%   • stops automatically when too few bins survive
%
% INPUTS
%   X         : nBin × nFrm × M   – STFT of all microphones
%   mic_xyz   : M   × 3           – microphone coordinates [m]
%   vs_xyz    : N   × 3           – search‑grid coordinates [m]
%   f         : nBin × 1          – frequency vector [Hz]
%   s_amount  : scalar            – # of source positions to pick
%   params  : struct with fields (defaults shown)
%       .muTH    = 0.8    coherence threshold
%       .lamTH   = 0      local‑SNR threshold (dB)
%       .B       = 2      IPD band‑width in #radians  (≈ ±B·π / f_max)
%       .min_d   = 0.20   minimum spacing between outputs [m]
%       .max_src = 20     hard cap on #iterations
%
% OUTPUT
%   est_xyz  : K×3  found source positions, one per row
%
% CALL
%   est_xyz = ir_srp_iter(X, mic_xyz, vs_xyz, f, params);
% ─────────────────────────────────────────────────────────────────────
% ─── defaults ────────────────────────────────────────────────────────
import function.algo.irsrp.*
p  = @(n,v) ~isfield(params,n) || isempty(params.(n));
if p('muTH'),     params.muTH   = 0.8;   end        % coherence
if p('lamTH'),    params.lamTH  = 0;     end        % SNR (dB)
if p('B'),        params.B      = 2;     end        % IPD half‑width
if p('Lv'),       params.Lv     = 2;     end        % noise-only frame
if p('dTH'),      params.dTH  = 0.50;    end        % IPD threshold
if p('delta_step'), params.delta_step  = 0.05; end  % δ step size
if p('min_d'),    params.min_d  = 0.20;  end        % [m]
if p('max_src'),  params.max_src = 20;   end

% ─── constants & mic‑pair bookkeeping ───────────────────────────────
c = 340;
[l_idx,q_idx] = find(triu(ones(size(mic_xyz,1)),1));   % upper‑tri pairs
pair_amount   = numel(l_idx);

% TDOA LUT τ (Ngrid × Npair)
tauTbl = tdoa_table(mic_xyz, vs_xyz, c, l_idx, q_idx);

% ─── two‑stage T‑F mask (Γ ∧ Λ) ─────────────────────────────────────
Wmask = tf_mask(X, params.muTH, params.lamTH, params.B, params.Lv);         % nBin × nFrm

% ─── iterative loop ────────────────────────────────────────────────
nBin = size(X,1);
est_xyz = zeros(params.max_src,3);   % pre‑allocate
src_cnt = 0;
sum_bin = 0;

while src_cnt < s_amount
    sum_bin = sum_bin + nnz(Wmask);
    % 1) SRP map on the *current* T‑F support ------------------------
    % srp_grid = srp_map_masked(X, f, tauTbl, l_idx, q_idx, Wmask);
    srp_grid = srp_map_masked_active(X, f, tauTbl, l_idx, q_idx, Wmask);
    [~,best] = max(srp_grid);        % single frame? we summed inside
    
    % nothing left? --------------------------------------------------
    if srp_grid(best)==0  
        break;  
    end
    
    % 2) store peak if spacing allows -------------------------------
    cand = vs_xyz(best,:);
    if src_cnt==0 || all(vecnorm(est_xyz(1:src_cnt,:) - cand,2,2) ...
                          > params.min_d)
        src_cnt = src_cnt + 1;
        est_xyz(src_cnt,:) = cand;
    end
    
    % 3) prune bins whose IPD matches this source -------------------
    Wmask = ipd_prune(Wmask, X, f, tauTbl(best,:), ...
                      l_idx, q_idx, params.dTH, params.delta_step);
    
    % stop conditions: empty mask OR hard cap -----------------------
    if ~any(Wmask(:)) || src_cnt >= params.max_src,  break;  end
end
est_xyz = est_xyz(1:src_cnt,:);      % trim
end

function tauTbl = tdoa_table(mic_xyz, vs_xyz, c, l_idx, q_idx)
d = pdist2(vs_xyz, mic_xyz,'euclidean');         % Ngrid × M
tauTbl = (d(:,q_idx) - d(:,l_idx)) / c;          % Ngrid × Npair
end
