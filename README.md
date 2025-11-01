# MD-OMP Monte Carlo (Varying Reverberation Time)

MATLAB implementation of the Monte Carlo simulation in **Subsection V-C** of the MD-OMP paper, evaluating localization accuracy under varying reverberation times (RT) for **MD-OMP** and several baselines (WB-OMP, NF-SRP-PHAT, IR-SRP, G-IRLS). Using the averaged success rates stored in `LE`, you can reproduce the line plot analogous to **Fig. 3** (success rate vs RT).

## Features
- Perimeter microphone layout & 3D candidate grid (30×45×5).
- Wideband processing over 125–2000 Hz with 61 bins.
- Multiple methods (toggle via `methods` struct).
- Optional 3D plotting of room, microphones, true and estimated sources.
- Automatic result saving with incremental filenames.

## Repo Structure (suggested)
```
.
├─ demo_monte_t60.m           % main script
├─ +function/
│  ├─ +pos/        % perimetermic(), generate_grid()
│  ├─ +dict/       % greens_matrix()
│  ├─ +algo/       % md_omp.m, wb_omp.m, srp_phat.m, ir_srp.m, g_irls.m, dsearchl10n.m
│  │  └─ +irsrp/   % tf_mask.m, ipd_prune.m, srp_map_masked_active.m, ...
│  └─ +plt/        % plot_room.m
├─ dataset/
│  ├─ rir/         % rir_rt0{2..8}_{1..100}.mat  (RIRs & source coords)
│  └─ speech/      % speech1.wav ... speech4.wav
└─ result/         % output *.mat (auto-created)
```

## Requirements
- MATLAB R2022b+ (uses `pagemtimes` and `pagetranspose`; R2024a used in the paper timing table).  
- Signal Processing Toolbox (for `stft`, `awgn`, etc.).

## Quick Start
1. Place the RIR files and four speech files as shown above.
2. Open `demo_monte_t60.m`, then set:
   ```matlab
   rt60_codes = 2:8;     % RT in 0.2–0.8 s
   n_trials   = 100;     % Monte Carlo trials
   methods = struct('mdomp',true,'wbomp',true,'nfsrp',true,'irsrp',true,'girls',false);
   if_save  = true;      % save results into /result/rt_testXX.mat
   if_plot  = false;     % per-trial 3D plot off by default
   plot_method = 'mdomp'; % when plotting: 'mdomp'|'wbomp'|'nfsrp'|'irsrp'|'girls'
   ```
3. Run:
   ```matlab
   >> demo_monte_t60
   ```
4. Results are saved to `result/rt_testXX.mat` with a struct `LE`:
   ```matlab
   LE.mdomp, LE.wbomp, LE.nfsrp, LE.irsrp, LE.girls
   ```

## Reproduce the Fig.-3 Style Plot (Success Rate vs RT)
```matlab
load('result/rt_test01.mat','LE');            % pick your file
rt = 0.2:0.1:0.8;                              % matches rt60_codes = 2:8
plot(rt, LE.mdomp, '-o'); hold on;
plot(rt, LE.wbomp, '-o'); plot(rt, LE.nfsrp, '-o'); plot(rt, LE.irsrp, '-o');
xlabel('RT (s)'); ylabel('Success rate'); ylim([0 1]); grid on; legend show;
```
> If multiple methods were enabled, each `LE.<method>` will have one success-rate per RT—sufficient to draw the multi-curve figure like Fig. 3.

## Optional 3D Plot of a Trial
Enable plotting in the script:
```matlab
if_plot = true; plot_method = 'mdomp';
```
The call `plt.plot_room(mic_xyz, source, est_pt, room, plot_cfg)` draws the room box, microphones, ground-truth sources, and the selected method’s estimates.

## Citation
If you use this code, please cite the MD-OMP paper:
> W.-T. Lai *et al.*, “Sound Source Localization using Multi-Dictionary Orthogonal Matching Pursuit in Reverberant Environments,” 2025.

## License
MIT (or your preferred license).