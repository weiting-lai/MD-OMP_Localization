% Monte Carlo Tests — MD-OMP under Varying Reverberation Time (RT)
%
% This script reproduces the Monte Carlo simulation described in Subsection 
% V-C¹ of the referenced paper, evaluating localization performance under 
% varying reverberation times (RTs) using the proposed Multi-Dictionary 
% Orthogonal Matching Pursuit (MD-OMP) method and baseline algorithms.
%
% What this script does (key settings aligned with the paper):
%   • Microphones: perimeter configuration (via pos.perimetermic);
%   • Virtual candidate grid: 30×45×5 over x–y–z (via pos.generate_grid),
%     then shifted by room_mid to match the room coordinates.
%   • Frequency band and binning: 125–2000 Hz with 61 uniformly spaced
%     bins (derived from the current window/FFT settings).
%   • Signals: 4 normalised speech sources (2 s segment each) convolved
%     with RIRs; 10 dB AWGN added to the microphone signals.
%   • STFT: Hamming window length = 1024, 50% overlap, FFT length = 2048.
%   • Methods (toggle in `methods`): MD-OMP, WB-OMP, NF-SRP-PHAT, IR-SRP, 
%     G-IRLS. (For MD-SBL, please see https://github.com/gerstoft/SBL)
%   • Metric: success rate = fraction of sources with position error
%     < 0.2 m (‖ŷ − y‖₂ < 0.2 m). Using the final averaged success rates 
%     stored in LE, you should be able to reproduce the line chart shown in 
%     Fig. 3 of the paper¹.
%
% Reference:
% ¹ W.-T. Lai et al., “Sound Source Localization using Multi-Dictionary
%   Orthogonal Matching Pursuit in Reverberant Environments”
%
% Author:   Weiting Lai
% Date:     01-11-25
close all; clear
import function.*
addpath 'dataset/speech'
addpath 'dataset/rir'
%% Read room data
load('rir_rt02_1.mat');
c = 340; % m/s
fs = 16000;
room = [3.36, 4.48, 4.1];
room_mid = [1.68, 2.24, 1.5];
rt60_codes = 2:8;  % RT in 0.2–0.8 s
n_trials = 100;  % Monte Carlo trials
if_save = true;  % Enable/disable saving results to .mat file
if_plot = false; % Enable/disable plotting of localization results
plot_method = 'mdomp'; % when plotting: 'mdomp'|'wbomp'|'nfsrp'|'irsrp'|'girls'
%% Select the localisation methods
methods = struct( ...
    'mdomp', true, ...   % MD-OMP
    'wbomp', true, ...  % WB-OMP
    'nfsrp', true, ...  % NF-SRP-PHAT
    'irsrp', true, ...  % IR-SRP
    'girls', false ...   % G-IRLS
);
%% Create microphones (fixed random wall)
mic_xyz = pos.perimetermic();
mic_amount = size(mic_xyz,1);
%% Create virtual source points
[vs_xyz, vs_amount] = pos.generate_grid([-1.4 1.4 30], [-2.0 2.0 45],...
    [-0.2 0.2 5]);
vs_xyz = vs_xyz + room_mid;
%% Find the greens function matrix -> SW
freq_l = 125;
freq_h = 2000;
numbin = 8*(freq_h-freq_l)/250+1;
S_matrix = zeros(mic_amount,vs_amount*numbin);
c = 340;
window = 1024;
nfft = window*2;

k_l = 2*pi*freq_l/c;
k_h = 2*pi*freq_h/c;
k_n = linspace(k_l, k_h, numbin);
freq_n = linspace(freq_l, freq_h, numbin);

bin_l = freq_l*window*2/fs;
bin_h = freq_h*window*2/fs;
% k = linspace(k_l,k_h,(bin_h-bin_l+1));
S_row = 1:mic_amount;
S_col = 1:vs_amount;

for k = k_n
    temp = dict.greens_matrix('pointsource', k, vs_xyz, mic_xyz);
    S_matrix(S_row,S_col) = temp;
    S_col = S_col + vs_amount;    
end

% normalisation
A = [];
[m,n]=size(S_matrix);
for i=1:n
    A(1,i)=norm(S_matrix(:,i));
end
A = repmat(A,m,1);
S_matrix = S_matrix./A;
S_matrix = reshape(S_matrix,[mic_amount,vs_amount,numbin]);

%% Initialize localization errors
LE = struct();   % averaged success rates across Monte Carlo tests
le = struct();   % temporary storage of a single trial

if methods.mdomp, LE.mdomp = []; end
if methods.wbomp, LE.wbomp = []; end
if methods.nfsrp, LE.nfsrp = []; end
if methods.irsrp, LE.irsrp = []; end
if methods.girls, LE.girls = []; end
for rt60 = rt60_codes
    if methods.mdomp, le.mdomp = []; end
    if methods.wbomp, le.wbomp = []; end
    if methods.nfsrp, le.nfsrp = []; end
    if methods.irsrp, le.irsrp = []; end
    if methods.girls, le.girls = []; end

    for h = n_trials
        filename = sprintf('rir_rt0%d_%d.mat', rt60, h);
        data = load(filename);
        rir_s1 = data.rir_s1;
        rir_s2 = data.rir_s2;
        rir_s3 = data.rir_s3;
        rir_s4 = data.rir_s4;
        source = data.source;
        %% Generate Speech Signal
        fs = 16000;
        % source1
        s1 = audioread("speech1.wav");
        source1_p = s1(48001:128000);
        
        % source2
        s2 = audioread("speech2.wav");
        source2_p = s2(48001:128000);
        
        % source3
        s3 = audioread("speech3.wav");
        source3_p = s3(48001:128000);
        
        % source4
        s4 = audioread("speech4.wav");
        source4_p = s4(48001:128000);
        
        % normalisation
        source1_p = source1_p./norm(source1_p);
        source2_p = source2_p./norm(source2_p);
        source3_p = source3_p./norm(source3_p);
        source4_p = source4_p./norm(source4_p);
        
        %% Simulate microphone recordings
        rec_n = fs*6-1;
        rec_s1 = zeros(rec_n, mic_amount);
        rec_s2 = zeros(rec_n, mic_amount);
        rec_s3 = zeros(rec_n, mic_amount);
        rec_s4 = zeros(rec_n, mic_amount);
        % Convolution
        for i = 1:mic_amount
            % Convolve speech with RIRs to generate each source recording
            rec_s1(:,i) = conv(source1_p, rir_s1(i,1:16000));
            rec_s2(:,i) = conv(source2_p, rir_s2(i,1:16000));
            rec_s3(:,i) = conv(source3_p, rir_s3(i,1:16000));
            rec_s4(:,i) = conv(source4_p, rir_s4(i,1:16000));
        end
        rec = rec_s1 + rec_s2 + rec_s3 + rec_s4; % Sum over all sources
        rec = rec';
        
        rng('default')
        rec = awgn(rec,10, 'measured'); % Add 10 dB AWGN
        
        %% STFT
        w = 16001;
        timeframe = w:(w+32000-1);
        s = stft(rec(1,timeframe),fs,'Window',hamming(window,'periodic'),...
            'OverlapLength',window/2,'FFTLength',nfft);
        T = size(s,2);
        rec_stft = zeros(mic_amount,T,nfft);
        rec_stft(1,:,:) = s';
        for i = 2:mic_amount
            [s,freq,t] = stft(rec(i,timeframe),fs,'Window',hamming(window,'periodic'),...
                'OverlapLength',window/2,'FFTLength',nfft);
            rec_stft(i,:,:) = s';
        end
        
        
        %% Frequency bin selection
        rec_subband = zeros(mic_amount, T, numbin);
        bin_n = linspace(bin_l,bin_h,numbin);
        rec_k = 1;
        for bin = bin_n
            % Extract the RTF at the selected frequency bin
            R = rec_stft(:,:,bin+window);
            rec_subband(:,:,rec_k) = R;
            rec_k = rec_k + 1;
        end
        
        %% Methods
        % MD-OMP
        if methods.mdomp
            [~, pti] = algo.md_omp(rec_subband, S_matrix, 4);
            mdomp_pt = vs_xyz(pti,:);
            [~, dist] = algo.dsearchl10n(source(:,:), mdomp_pt(:,:));
            le.mdomp = [le.mdomp, sum(dist < 0.2)/4];
        end
        
        % WB-OMP
        if methods.wbomp
            rtf_avg = mean(rec_subband, 2);
            [~, pti] = algo.wb_omp(rtf_avg, S_matrix, 4);
            wbomp_pt = vs_xyz(pti,:);
            [~, dist] = algo.dsearchl10n(source(:,:), wbomp_pt(:,:));
            le.wbomp = [le.wbomp, sum(dist < 0.2)/4];
        end

        % NF-SRP-PHAT
        if methods.nfsrp
            X = permute(rec_subband, [3,2,1]);
            nfsrp_pt = algo.srp_phat(X, mic_xyz, vs_xyz, freq_n', 4, 0.2);
            [~, dist] = algo.dsearchl10n(source(:,:), nfsrp_pt(:,:));
            le.nfsrp = [le.nfsrp, sum(dist < 0.2)/4];
        end

        % IR-SRP
        if methods.irsrp
            params = struct('muTH',0.4, 'lamTH',0, 'B',2, 'Lv',2, ...
                            'min_d',0.2, 'max_src',10);
            X = permute(rec_subband, [3,2,1]);
            irsrp_pt = algo.ir_srp(X, mic_xyz, vs_xyz, freq_n', 4, params);
            [~, dist] = algo.dsearchl10n(source(:,:), irsrp_pt(:,:));
            le.irsrp = [le.irsrp, sum(dist < 0.2)/4];
        end

        % G-IRLS
        if methods.girls
            rtf_avg = mean(rec_subband, 2);
            girls_pt = algo.g_irls(rtf_avg, S_matrix, 0.0001, 10, 1, 0.01, 10, vs_xyz, 0.2);
            [~, dist] = algo.dsearchl10n(source(:,:), girls_pt(:,:));
            le.girls = [le.girls, sum(dist < 0.2)/4];
        end

    end % trial loop

    if methods.mdomp, LE.mdomp = [LE.mdomp, sum(le.mdomp)/...
            max(1, numel(le.mdomp))]; end
    if methods.wbomp, LE.wbomp = [LE.wbomp, sum(le.wbomp)/...
            max(1, numel(le.wbomp))]; end
    if methods.nfsrp, LE.nfsrp = [LE.nfsrp, sum(le.nfsrp)/...
            max(1, numel(le.nfsrp))]; end
    if methods.irsrp, LE.irsrp = [LE.irsrp, sum(le.irsrp)/...
            max(1, numel(le.irsrp))]; end
    if methods.girls, LE.girls = [LE.girls, sum(le.girls)/...
            max(1, numel(le.girls))]; end

end % rt60 loop
%% Optional: Save the success rate (LE) of the algorithms 
if if_save
    result_dir = 'result';
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    % Automatically find next available filename
    base_name = 'rt_test';
    i = 1;
    filename = fullfile(result_dir, sprintf('%s%02d.mat', base_name, i));
    
    while exist(filename, 'file')
        i = i + 1;
        filename = fullfile(result_dir, sprintf('%s%02d.mat', base_name, i));
    end
    
    % Save the LE struct
    save(filename, 'LE');
    fprintf('Results saved to %s\n', filename);
end

%% Optional: Plot the localisation result of a single trial
plot_cfg = struct( ...
    'enable', if_plot, ...              % toggle plotting on/off
    'selected_method', plot_method, ... % one of: 'mdomp'|'wbomp'|'nfsrp'|'irsrp'|'girls'
    'save_path', '', ...                % e.g., 'fig_room.png' ('' = no save)
    'view', [-35 20], ...               % azimuth, elevation
    'show_roi', true ...                % draw ROI box from vs_xyz
);
est_pt = eval([plot_cfg.selected_method '_pt']);
if plot_cfg.enable
    plt.plot_room(mic_xyz, source, est_pt, room, plot_cfg);
end

