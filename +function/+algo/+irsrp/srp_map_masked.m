function srp_sum = srp_map_masked(X, f, tauTbl, l_idx, q_idx, W)
% Function  : srp_map_masked
% Creation  : 24-07-25
% Author    : Weiting Lai
% Version   : -
% Computes SRP-PHAT map with time-frequency masking (e.g. coherence or SNR mask)
% using precomputed TDOA lookup table (tauTbl)
% -------------------------------------------------------------
% 
% INPUTS
% X       : nBin × nFrm × M         complex STFT of all microphones
% f       : nBin × 1                frequency vector [Hz]
% tauTbl  : nGrid × nPair           TDOA lookup table
% l_idx   : nPair × 1               first mic index for each pair
% q_idx   : nPair × 1               second mic index for each pair
% W       : nBin × nFrm × nPair     T-F mask (e.g. logical or soft mask)
%
% RETURN
% srp_sum : nGrid × 1               SRP response for each candidate location
[nBin, ~, ~] = size(X);
nGrid = size(tauTbl,1);
srp_sum = zeros(nGrid,1,'single');

% vectorised exponentials   E: nBin × nPair × nGrid
E = exp( 2j*pi*f(:) .* permute(tauTbl,[3 2 1]) );

for p = 1:numel(l_idx)
    P = X(:,:,l_idx(p)) .* conj( X(:,:,q_idx(p)) );
    P = P ./ (abs(P) + eps);          % PHAT
    P = P .* W(:,:,p);                       % apply mask
    
    Ep = reshape(E(:,p,:), nBin, []); % nBin × nGrid
    spec = real(P).' * real(Ep) - imag(P).' * imag(Ep);  % nFrm × nGrid
    srp_sum = srp_sum + sum(spec,1).';                   % sum over frames
end
end
