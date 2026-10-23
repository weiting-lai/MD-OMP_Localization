# MD-OMP · 3D Source Localisation

This is a companion to **Sound Source Localization using Multi-Dictionary Orthogonal Matching Pursuit in Reverberant Environments**. The research contribution is MD-OMP localisation: estimating the 3D positions of multiple sound sources from microphone recordings.

The room view, measured/estimated coordinates and spatial errors show the localisation result. The source-extraction example and stereo player were added for this companion page to illustrate a downstream use of the estimated locations. They are separate from the paper's proposed method and localisation evaluation.

## Current listening example

Only one extraction pipeline is published:

```text
14-channel recording
  -> WPE dereverberation
  -> MD-OMP-position-guided spatially regularised ILRMA
  -> mask-based MVDR
  -> four saved source tracks
  -> centred two-channel stereo playback
```

The estimated coordinates guide initialisation, a prior during iterative extraction and geometric association of the resulting tracks. Ground-truth coordinates are used for source labels and localisation-error evaluation. Clean speech, reference activity masks and measured room-transfer filters are not used to generate the published source audio.

The four source outputs come from the **MD-OMP estimated-position A configuration without the optional postfilter**, with fixed source order and mean STOI **0.704360** on this recording. The per-source values are **0.647388, 0.767717, 0.785845, 0.616492**. This auxiliary extraction score is separate from the paper's localisation metrics. Distortion and interference remain, particularly for S1/S4. Results concern this clip and are not a new speech-enhancement benchmark.

The exact selected processing settings and output measurements are in `data/provenance.json` and `data/listening-evaluation.json`. Historical comparisons remain in the local research folders; they are excluded from the static site.

## Playback

- Click an estimated marker, a Speaker button or a table audition button to hear that source.
- **All Four Extracted Sources** plays their fixed-gain stereo sum.
- **Original Recording** provides the actual microphone-1 mixture for context.
- The six audio buffers share one playback clock. Switching preserves time.
- Keys **1–4** select sources; **Space** plays or pauses. Loop, mute, level and seek controls are available.
- The room camera has 3D/top views, drag orbiting and arrow-key controls. It affects only the visual view.

All published WAVs are 16 kHz, PCM16, two channels, 201728 frames (12.608 seconds). Each source's left and right channels are bit-exact copies of its verified mono output. The mix is `L = R = round((S1 + S2 + S3 + S4) / 4)` in PCM16. This centred stereo presentation keeps the room coordinates as a visual localisation result.

## Local preview

From the original experiment directory:

```powershell
python md_omp_demo/pipeline/serve.py
```

Open [the local preview](http://127.0.0.1:8765/). The public website is buildless HTML/CSS/JavaScript and needs no runtime processing server. Open it through HTTP so the JSON and WAV assets load correctly.

## Rebuild and verify

The local audio study must already exist at `gt_position_study/`. The selected `wpe-ilrma-euclidean0.05-n4096-k2-md-omp-A-mvdr-p0-gmatch` candidate and its measured WAV results are the authoritative inputs. The original `demo_journal.m`, recordings and localisation export remain unchanged.

```powershell
python md_omp_demo/pipeline/publish_localisation.py
python md_omp_demo/pipeline/validate.py
node --check md_omp_demo/app.js
python md_omp_demo/pipeline/package_site.py
```

`publish_localisation.py` republishes the measured outputs and the HTML/player templates without rerunning separation. It archives retired listening assets locally. `validate.py` checks the actual stereo channels, waveform identity, fixed-gain mix, STOI, estimated coordinates, relative links, original script hash and the exact MetaFOA base stylesheet. Numerical validation requires the local clean playback clips for evaluation only.

The static deployment files are in `dist/site/`; `dist/md-omp-demo.zip` contains the same allowlisted files. Copy that folder's contents into the intended GitHub Pages directory. Packaging excludes the research code, calibration models, alternate listening methods and local recordings. No remote publication is performed by these commands.

The base `styles.css` is copied verbatim from MetaFOA Git blob `de61f1aa7556575d903fd94456c1c484c737d987`. Page-specific layouts are in `demo.css`; the palette and Inter / IBM Plex Mono typography are retained.

## Paper

[Research repository](https://github.com/weiting-lai/MD-OMP_Localization) · [DOI](https://doi.org/10.1109/TASLPRO.2025.3624966)

```bibtex
@article{mdomp,
  author = {Lai, Wei-Ting and Birnie, Lachlan and Abhayapala, Thushara and Bastine, Amy and Samarasinghe, Prasanga},
  title = {Sound Source Localization using Multi-Dictionary Orthogonal Matching Pursuit in Reverberant Environments},
  journal = {IEEE Transactions on Audio, Speech and Language Processing},
  year = {2025},
  doi = {10.1109/TASLPRO.2025.3624966}
}
```

The paper, supplied setup photograph and speech recording are provided for this research demonstration. Existing repository licensing applies to the research code.
