# Bayesian Parameter Estimation of GW231002 and GW231123: Evidence for Intermediate-Mass Black Holes

Code accompanying my M.Sc. thesis (Physics, MANIT Bhopal, 2026), supervised by
Dr. Sumit Ghosh and Dr. Sadhana Singh.

## Summary

Bayesian parameter estimation on the LVK O4 events **GW231002** and
**GW231123** using [`bilby`](https://git.ligo.org/lscsoft/bilby) with the
`dynesty` dynamic nested sampler, run against three waveform approximants
(`IMRPhenomXPHM`, `IMRPhenomTPHM`, `IMRPhenomPv3`) per event to test the
robustness of the inferred source parameters to waveform systematics. Both
events show component masses consistent with the pair-instability mass gap
(60–130 M☉) and remnant masses in the intermediate-mass black hole (IMBH)
range, favoring a hierarchical-merger formation channel.

| | GW231002 | GW231123 |
|---|---|---|
| Primary mass (M☉) | 93 | 205 |
| Secondary mass (M☉) | 79 | 100 |
| Remnant mass (M☉) | 167 | 300 |
| χ_eff | −0.14 | 0.25 |
| Luminosity distance (Gpc) | 4.1 | 2.9 |

Full tables: [`results/GW231002/summary_table.csv`](results/GW231002/summary_table.csv),
[`results/GW231123/summary_table.csv`](results/GW231123/summary_table.csv),
[`results/comparison_table.csv`](results/comparison_table.csv).

## Repository layout

```
priors/           prior distributions (one file per event, shared across all 3 waveforms)
src/
  event_config.py per-event GPS time, data paths, likelihood/sampler settings
  data_utils.py   GWOSC strain loading
  remnant_fits.py remnant-mass/spin/ringdown fits + posterior summary helpers
  run_pe.py       single entry point for any (event, waveform) run
notebooks/        postprocessing: overlaid corner plots, regenerated summary tables
results/          summary tables + individual/combined plots per event
```

The original analysis was six separate near-duplicate scripts (one per
event × waveform). They've been consolidated here into one parameterized
pipeline (`src/run_pe.py`) plus two prior files, since the only real
differences between runs were the waveform approximant and the priors/
likelihood settings per event — everything else (data loading, PSD
estimation, sampler config, remnant fits, plotting) was identical.

## Setup

```bash
conda env create -f environment.yml
conda activate gw-pe-imbh
# or: pip install -r requirements.txt
```

Strain data is not included (GWOSC files are large and freely available).
Download the relevant 4096 s segments from [GWOSC](https://gwosc.org) for:

- GW231002 — GPS 1380292774.7 (H1 file starts at GPS 1380290560)
- GW231123 — GPS 1384782888.7 (H1 file starts at GPS 1384779776)

and place them under `data/`, or point to them directly with `--h1-file` /
`--l1-file`.

## Running an analysis

```bash
python src/run_pe.py \
    --event GW231002 \
    --waveform IMRPhenomXPHM \
    --h1-file data/H-H1_GWOSC_O4a_16KHZ_R1-1380290560-4096.hdf5 \
    --l1-file data/L-L1_GWOSC_O4a_16KHZ_R1-1380290560-4096.hdf5
```

Repeat with `--waveform IMRPhenomTPHM` / `IMRPhenomPv3` for the other two
runs used for the waveform-systematics comparison. Each run produces a
`outdir_<event>_<waveform>/` with the bilby result file, plus strain/PSD/
corner plots under `results/<event>/`.

Note: dynesty runs at `nlive=3000` take many hours on an HPC cluster (this
work used the MANIT Bhopal CCF HPC facility) — not intended to run on a
laptop in a reasonable time.

## Methodology (brief)

- Data conditioned with an 8 s analysis window (6 s pre-trigger, 2 s
  post-trigger), PSD estimated via Welch/median method over 256 s of
  off-source data (Tukey-windowed).
- Priors: `BBHPriorDict`/`PriorDict` with astrophysically motivated ranges
  for chirp mass, mass ratio, spins, sky location, and
  `UniformSourceFrame` for luminosity distance (see `priors/`).
- Likelihood: `GravitationalWaveTransient` with time (and, for GW231002,
  phase) marginalization.
- Sampler: `dynesty` via bilby, `nlive=3000`, `dlogz≈0.1`, `sample="rwalk"`.
- Remnant mass/spin/ringdown frequency from NR-calibrated phenomenological
  fits (`src/remnant_fits.py`), not full-waveform ringdown fitting.

Full derivations, priors table, and discussion are in the accompanying
M.Sc. thesis (not included in this repository).

## Author

Lakshita Jain — M.Sc. Physics, MANIT Bhopal
[LinkedIn](https://www.linkedin.com/in/lakshitajain) ·
[GitHub](https://github.com/Lakshita-Jain09)

## License

MIT — see [`LICENSE`](LICENSE).
