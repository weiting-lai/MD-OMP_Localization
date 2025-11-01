function estimated_pt = srp_phat(X, mic_xyz, vs_xyz, f, s_amount, min_d)
% Function    : srp_phat
% Creation    : 24-07-25
% Author      : Weiting Lai
% Version     : 2.0
% Description : Optimize GCC for-loop to improve execution speed
% -------------------------------------------------------------
%
% INPUTS
%   X         : nBin × nFrm × M   – STFT of all microphones
%   mic_xyz   : M   × 3           – microphone coordinates [m]
%   vs_xyz    : N   × 3           – search‑grid coordinates [m]
%   f         : nBin × 1          – frequency vector [Hz]
%   s_amount  : scalar            – # of source positions to pick
%   min_d     : scalar            – minimum distance between picks [m]
%
% RETURN
%   estimated_pt : s_amount × 3   – selected grid points (rows of vs_xyz)

%% Constants and microphone‑pair indices
c = 340;    % sound speed [m/s]
[l_idx,q_idx] = find(triu(ones(size(mic_xyz,1)),1));  % upper‑tri indices
pair_amount   = numel(l_idx);

%% (1) Time‑delay look‑up table  τ  (N × P) – fully vectorised
dists = pdist2(vs_xyz, mic_xyz,'euclidean');      % N × M
tau   = (dists(:,q_idx) - dists(:,l_idx)) / c;    % N × P      Δt/c

%% (2) SRP‑PHAT map – accumulate every pair
[nBin,nFrm,~] = size(X);
srp = zeros(size(vs_xyz,1), nFrm,'single');

% Exponentials for **all** grid points & pairs:  E = e^{j2π f τ}
E = exp(2j*pi*f(:) .* permute(tau,[3 2 1]));      % nBin × P × N

for p = 1:pair_amount
    % Normalised cross‑spectrum for this pair (ϕ‑hat)
    P = X(:,:,l_idx(p)) .* conj( X(:,:,q_idx(p)) );
    P = P ./ (abs(P) + eps);                      % avoid division by 0
    
    E_p  = reshape(E(:,p,:), nBin, []);      % nBin × nGrid   (works even if nGrid==1)
    spec = real(P).' * real(E_p) - imag(P).' * imag(E_p);   % nFrm × nGrid
    srp  = srp + spec.';                                    % accumulate (Ngrid × nFrm)
end

%% (3) Pick the strongest grid points with spacing > min_d
srp_sum      = sum(srp,2);                         % collapse over frames
[~,index]    = sort(srp_sum,'descend');

% initialisation min_d & est_pt
pti = zeros(s_amount,1);
estimated_pt = zeros(s_amount,3);
pti(1) = index(1);
estimated_pt(1,:) = vs_xyz(pti(1),:);

% remove points shorter than min_d
left = 1;
right = 2;
while left < s_amount
    if_allowed = 1;
    for i = 1:left
        % distance between found points and the next largest point
        dist = norm(estimated_pt(i,:)-vs_xyz(index(right),:));
        if(dist < min_d)
            if_allowed = 0;
            break;
        end
    end
    % store new point
    if(if_allowed)
        left = left + 1;
        pti(left) = index(right);
        estimated_pt(left,:) = vs_xyz(pti(left),:);
    end
    right = right + 1;
end
end
