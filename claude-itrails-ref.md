# iTRAILS Reference Guide

Reference for **iTRAILS** ("Tree reconstruction of ancestry using incomplete lineage sorting"), a command-line tool for inferring ancestral population parameters and decoding gene-tree topologies along a genome alignment using a coalescent hidden Markov model (the TRAILS framework, [Rivas-González et al., 2024](https://doi.org/10.1371/journal.pgen.1010836)).

**Sources**: [itrails.readthedocs.io](https://itrails.readthedocs.io/en/docs-stable/) · [github.com/trails-phylogeny/itrails](https://github.com/trails-phylogeny/itrails)

---

## What it does

Given a four-genome alignment — three ingroup species **A**, **B**, **C** with topology ((A,B),C) plus an **outgroup** — iTRAILS:

1. **`itrails-optimize`** — fits the coalescent HMM by maximum likelihood: estimates ancestral effective population sizes and speciation times from the spatial pattern of incomplete lineage sorting (ILS) along the alignment.
2. **`itrails-viterbi`** — decodes the single most likely sequence of hidden states (gene-tree topology + discretized coalescence-time intervals) along the alignment.
3. **`itrails-posterior`** — computes, per position, the posterior probability of every hidden state.

Hidden states combine a topology (((A,B),C) with coalescences in the AB or ABC ancestor, or the ILS topologies ((A,C),B) / ((B,C),A)) with discretized time intervals for the first and second coalescent event. With `n_int_AB=3, n_int_ABC=3` there are 27 states (index 0–26).

An introgression-aware variant exists as separate entry points: `itrails-int-optimize`, `itrails-int-viterbi`, `itrails-int-posterior` (adds migration parameters `N_BC`, `t_m`, `m`; see `example_config_int.yaml` in the package).

---

## Installation

```bash
pip install itrails
# or
conda install conda-forge::itrails
```

- This document describes the **2.x** series (verified against 2.0.5), matching the docs-stable documentation. The old 0.1.x releases on PyPI have an incompatible CLI (`--output` is a directory, fixed output names `best_model.yaml`/`viterbi.csv`, no `--config-file`/`--n_cpu`/`--reference` on the decode commands) — pin `itrails >=2.0.5,<3`.

- Requires **Python ≥ 3.12**. Pinned deps: `numpy==1.26.4`, `scipy==1.13.0`, `numba==0.59.1`, `joblib`, `biopython`, `pyyaml`, `h5py`, `pandas`, `ray[default]`.
- Parallelism via ray/joblib. CPU use is capped at `min(n_cpu, SLURM_JOB_CPUS_PER_NODE or cpu_count())` — i.e. **Slurm-aware**: it never uses more cores than the job allocation, and clamps invalid/absent `n_cpu` to what's available.
- Sets `OMP/MKL/NUMEXPR/RAYON_NUM_THREADS=1` internally; parallelism is process-level.

---

## Input

A **MAF** (Multiple Alignment Format) file containing the aligned genomes. Sequence names in the MAF must match the entries of `species_list` (e.g. `["hg38", "panTro5", "gorGor5", "ponAbe2"]`). Order = A, B, C, outgroup. The parser (`itrails.read_data.maf_parser`) extracts alignment blocks for those species.

---

## Model parameters

Times are in **generations**; sizes are effective population sizes (individuals); rates are per site per generation. Internally everything is scaled by `mu` (values × `mu`, except `r` which becomes `r/mu`), but config files use natural units.

| Parameter | Meaning | Constraint |
|:---|:---|:---|
| `N_AB` | Ne of the AB ancestral population | fixed or optimized (required) |
| `N_ABC` | Ne of the ABC ancestral population | fixed or optimized (required) |
| `t_1` | Time from (ultrametric) sampling to the A–B split | see combination rules |
| `t_A`, `t_B`, `t_C` | Per-species time from its sampling to its join (A,B → AB split; C → ABC split). Use instead of `t_1` for non-ultrametric trees | see combination rules |
| `t_2` | Duration of the AB ancestral population (A–B split → ABC split) | fixed or optimized (required) |
| `t_3` | Duration from ABC split to join with the outgroup | at least one of `t_3` / `t_upper` |
| `t_upper` | Time from the start of the last discretized ABC interval to outgroup coalescence. If absent, computed from `t_3` and `N_ABC` via cutpoints | at least one of `t_3` / `t_upper` |
| `t_out` | Time from outgroup sampling to coalescence with ABC | optional; **must be fixed** if given |
| `r` | Recombination rate | fixed or optimized (required) |
| `mu` | Mutation rate | **must be fixed** (required) |

**Allowed combinations of time parameters** (anything else is an error; unspecified ones are derived): `{t_A,t_B,t_C}`, `{t_1,t_A}`, `{t_1,t_B}`, `{t_1,t_C}`, `{t_A,t_B}`, `{t_A,t_C}`, `{t_B,t_C}`, `{t_1}`.

**Discretization** (in `settings`): `n_int_AB`, `n_int_ABC` — number of time intervals in the AB and ABC ancestors (positive ints; 3 is the documented example). More intervals = more states = slower but finer-grained.

---

## `itrails-optimize`

```bash
itrails-optimize config.yaml --input aln.maf --output outdir/prefix
```

- `config_file` (positional, required): YAML config.
- `--input` (optional): MAF path — overrides `settings.input_maf` (warns if both).
- `--output` (optional): `directory/prefix` — overrides `settings.output_prefix`. The directory is created if missing.
- Input and output must each be given in the config or on the command line, else it errors.

### Config YAML

```yaml
fixed_parameters:
  mu: 2e-8                  # mu is always fixed

optimized_parameters:       # each: [starting, min, max]
  N_AB:    [50000, 5000, 500000]
  N_ABC:   [50000, 5000, 500000]
  t_1:     [240000, 24000, 2400000]
  t_2:     [40000, 4000, 400000]
  t_3:     [800000, 80000, 8000000]
  t_upper: [745069.3855, 74506.9385, 7450693.8556]
  r:       [1e-8, 1e-9, 1e-7]

settings:
  input_maf: path/to/alignment.maf     # optional if --input given
  output_prefix: path/to/outdir/prefix # optional if --output given
  n_cpu: 64                            # capped at Slurm allocation
  method: "Nelder-Mead"                # or "L-BFGS-B" (case-insensitive)
  species_list: ["hg38", "panTro5", "gorGor5", "ponAbe2"]
  n_int_AB: 3
  n_int_ABC: 3
  n_iter: 10000                        # optimizer iterations (default 10000)
```

A parameter may appear in `fixed_parameters` **or** `optimized_parameters`, never both. Starting values must lie within bounds; all values positive.

### Outputs (exact names, given `--output outdir/prefix`)

| File | Content |
|:---|:---|
| `outdir/prefix.starting_params.yaml` | Record of starting values/bounds and settings |
| `outdir/prefix.best_model.yaml` | Best parameter set so far: `fixed_parameters`, `optimized_parameters`, `results: {log_likelihood, iteration}`, `settings` (including resolved `input_maf` and `output_prefix`). **Written at start** (log-likelihood −inf) and updated during optimization — existence ≠ finished. |
| `outdir/prefix.optimization_history.csv` | Per-iteration parameter values, log-likelihood, elapsed time |

The `best_model.yaml` is the input config for the decoding steps.

---

## `itrails-viterbi` / `itrails-posterior`

```bash
itrails-viterbi   --config-file outdir/prefix.best_model.yaml \
                  --input aln.maf --output outdir/prefix
itrails-posterior --config-file outdir/prefix.best_model.yaml \
                  --input aln.maf --output outdir/prefix
```

- All arguments are flags (no positional). Config file optional — every model parameter and setting can instead be supplied/overridden on the command line: `--mu --t1 --t_A --t_B --t_C --t2 --t3 --t_upper --t_out --N_AB --N_ABC --r --n_cpu --species_list --reference --n_int_AB --n_int_ABC --cutpoints_AB --cutpoints_ABC`.
- `--reference SPECIES`: polarize output coordinates to one species' genome.
- Since `best_model.yaml` carries `input_maf`/`output_prefix` in its `settings`, the flags are optional — but passing `--input`/`--output` explicitly is more robust in pipelines.

### Outputs (exact names)

| Command | Files |
|:---|:---|
| `itrails-viterbi` | `prefix.viterbi.csv`, `prefix.hidden_states.csv` |
| `itrails-posterior` | `prefix.posterior.csv`, `prefix.hidden_states.csv` |

If `prefix.hidden_states.csv` already exists, the second command writes `prefix.hidden_states_2.csv` instead — **beware when declaring gwf outputs**; safest to use distinct prefixes per step (e.g. `prefix_vit`, `prefix_post`) so each target owns its files.

**`*.viterbi.csv`** — alignment segmented into blocks of constant most-likely state: `Block Index, Starting Position, End Position, Most Likely State`.

**`*.posterior.csv`** — per-position probability of each hidden state: block/position columns followed by one probability column per state.

**`*.hidden_states.csv`** — state index → interpretation: `State Index, Topology, Interval First Coalescence, Interval Second Coalescence, Shorthand Name` (e.g. `((A,B),C)` with bracket notation for the 2-species coalescence, shorthand like `1-1-2`).

---

## Pipeline / gwf integration notes

- **Chain**: `optimize` (MAF + config → `best_model.yaml`) → `viterbi` and `posterior` (each: `best_model.yaml` + MAF → CSVs). Declare `best_model.yaml` as the optimize target's output and the decode targets' input so gwf infers the dependency.
- **Resources**: optimize is the heavy step (many likelihood evaluations × `n_iter`); give it many cores and set `settings.n_cpu` = gwf `cores`. itrails clamps to `SLURM_JOB_CPUS_PER_NODE` anyway. Decoding is lighter but still parallel.
- **Progressive writes**: optimize writes/updates its outputs from iteration 0, so the write-to-tmp-then-move idiom doesn't apply; rely on gwf's job-status tracking for failure detection.
- **Same-name collision**: run viterbi and posterior with different `--output` prefixes to avoid the `hidden_states_2.csv` fallback making declared outputs wrong.
- **Test data**: example MAF on Zenodo: `https://zenodo.org/records/14930374/files/example_alignment.maf` (great apes; species `hg38 panTro5 gorGor5 ponAbe2`). Use a small `n_iter` (e.g. 5–50) for smoke tests.
- **Env**: name clash alert — the PyPI/conda package `itrails` vs. this repo also being named `itrails`. Never put the repo's parent dir on `sys.path` in a way that shadows the installed package.
