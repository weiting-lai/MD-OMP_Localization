function W_TF = tf_mask(X, muTH, lamTH, B, Lv)
% Function  : tf_mask
% Creation  : 24-07-25
% Author    : Weiting Lai
% Version   : -
% T‑F bin initialised selection  (Dang et al., ICCAIS 2023, eq.(7)‑(15))
% -------------------------------------------------------------
% 
% INPUTS
% X     : nBin × nFrm × M   complex STFT of all microphones
% muTH  : coherence threshold     (e.g. 0.7 ~ 0.9)
% lamTH : SNR threshold, dB       (e.g. -5  or  0)
% B     : half‑window size for time averaging  (typ. 1–3)
% Lv    : #noise‑only frames for P_vm(k) estimate (default = B)
%
% RETURN
% W_TF  : nBin × nFrm × nPair   logical mask
%         1 ⇒ bin passes both coherence & SNR tests

if nargin < 5, Lv = B; end                % default Lv
[nBin,nFrm,M] = size(X);
K       = 2*B + 1;                        % window length
nPair   = M*(M-1)/2;
W_TF    = false(nBin,nFrm,nPair,'logical');

pair = 0;
for m1 = 1:M-1
    for m2 = m1+1:M
        pair = pair + 1;

        % ---------- 1) time average (E{·}) ----------
        % cross‑spectrum
        C = movmean( X(:,:,m1) .* conj(X(:,:,m2)), K, 2);   % nBin×nFrm
        % auto‑power
        P1 = movmean( abs(X(:,:,m1)).^2 , K, 2);
        P2 = movmean( abs(X(:,:,m2)).^2 , K, 2);

        % ---------- 2) Magnitude‑squared coherence ----------
        mu = abs(C) ./ sqrt( P1 .* P2 + eps );              % eq. (7)

        % ---------- 3) Local SNR  λ -------------------------
        % noise power P_vm(k) : eq. (11)
        Pv1 = mean( abs(X(:,1:Lv,m1)).^2 , 2);              % nBin×1
        Pv2 = mean( abs(X(:,1:Lv,m2)).^2 , 2);

        % from broadcast to nFrm
        Pv1 = repmat(Pv1, 1, nFrm);
        Pv2 = repmat(Pv2, 1, nFrm);

        % signal power P_sm(k,l)
        Ps1 = abs(X(:,:,m1)).^2;
        Ps2 = abs(X(:,:,m2)).^2;

        % eq. (10): λ = min( Ps/Pv - 1 )
        lambda = min( Ps1./Pv1 - 1 , Ps2./Pv2 - 1 );
        lambda_dB = 10*log10( max(lambda,0) + eps );

        % ---------- 4) two-layer weighting function ----------
        CohOk = (mu      > muTH);
        SNROk = (lambda_dB > lamTH);

        W_TF(:,:,pair) = CohOk & SNROk;   % eq. (14)
    end
end
end
