# Demo of Sound Source Localization using Multi-Dictionary Orthogonal Matching Pursuit

This MATLAB demo aims to identify 3-D sound source positions in a reverberant room using a perimeter microphone array (arranged around the room walls with slight height differences between each). It serves as the implementation of the Monte Carlo simulation in Subsection V-C of the MD-OMP paper, evaluating localization accuracy under varying reverberation times (RT) for MD-OMP and several baselines (WB-OMP, NF-SRP-PHAT, IR-SRP, G-IRLS). Using the averaged success rates stored in LE, you can reproduce the line plot analogous to Fig. 3.

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
│  ├─ Download RIR dataset: [Zenodo link](https://zenodo.org/records/17500576) (RIR dataset hosted on Zenodo; see link below)
└─ result/         % output *.mat
```

## Requirements
- MATLAB R2022b+ (uses `pagemtimes` and `pagetranspose`; R2024a used in the paper timing table).  
- Signal Processing Toolbox (for `stft`, `awgn`, etc.).


## Citation
If you use this code or dataset, please cite the paper below:
> W.-T. Lai, L. Birnie, T. Abhayapala, A. Bastine, and P. Samarasinghe,
Sound Source Localization using Multi-Dictionary Orthogonal Matching Pursuit in Reverberant Environments, IEEE Trans. Audio, Speech, Lang. Process., 2025.

## BibTeX
```
@ARTICLE{mdomp,
  author={Lai, Wei-Ting and Birnie, Lachlan and Abhayapala, Thushara and Bastine, Amy and Samarasinghe, Prasanga},
  journal={IEEE Transactions on Audio, Speech and Language Processing},
  title={Sound Source Localization using Multi-Dictionary Orthogonal Matching Pursuit in Reverberant Environments},
  year={2025},
  volume={},
  number={},
  pages={1-11},
  keywords={Estimation;Location awareness;Time-frequency analysis;Accuracy;Matching pursuit algorithms;Direction-of-arrival estimation;Microphone arrays;Dictionaries;Transfer functions;Noise;source localization;source position estimation;sparse representation;matching pursuit;spatio-temporal-spectral fusion},
  doi={10.1109/TASLPRO.2025.3624966}
}```

## License
MIT