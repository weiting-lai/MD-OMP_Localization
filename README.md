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


## Citation
If you use this code, please cite the MD-OMP paper:
> W.-T. Lai *et al.*, “Sound Source Localization using Multi-Dictionary Orthogonal Matching Pursuit in Reverberant Environments,” 2025.

## License
MIT.