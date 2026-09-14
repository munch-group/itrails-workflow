
# Data files tracked with Git

- `simulated_alignment.maf` — 1 Mb dummy alignment of hg38/panTro5/gorGor5/ponAbe2
  simulated with msprime under a great-ape-like demography (splits at 240k/350k/700k
  generations, ancestral Ne 65k/100k/100k). Used by the default `example` analysis in
  `analyses.yml`; regenerate with `pixi run simulate-example`
  (see `scripts/simulate_data.py` for the exact parameters).

Not tracked (downloaded on demand):

- `example_alignment.maf` — the real 183 MB great-ape alignment from Zenodo;
  fetch with `pixi run download-example`.
