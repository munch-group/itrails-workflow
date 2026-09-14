# GWF (Grid WorkFlow) Reference Guide

Reference for **gwf**, a pragmatic Python workflow tool for building and running large scientific pipelines on HPC clusters. Workflows are plain Python; gwf infers the dependency graph from input/output filenames, submits only out-of-date work to a backend (Slurm, local, SGE, LSF, PBS), and is "fire-and-forget" (no `screen`/`tmux` needed).

**Source**: [gwf.app](https://gwf.app/) ([guide](https://gwf.app/guide/), [reference](https://gwf.app/reference/), [backends](https://gwf.app/backends/)) · [github.com/gwforg/gwf](https://github.com/gwforg/gwf)

Requires **Python 3.10+**. Originated at GenomeDK, Aarhus University.

---

## Mental Model

1. A **workflow** is a Python file (default `workflow.py`) that builds a `Workflow` object named `gwf`.
2. You add **targets** — named units of work, each with `inputs` (files it reads), `outputs` (files it writes), resource `options`, and a bash `spec`.
3. gwf **infers dependencies from files**: target B depends on target A if any of B's inputs is one of A's outputs. No manual edge declaration.
4. On `gwf run`, gwf computes the graph and **submits only targets that are out of date** to the active **backend**. A target is out of date when an output is missing, an input is newer than an output, or **the spec text changed** (gwf hashes each spec).
5. The same workflow runs unchanged on any backend — develop locally, submit to Slurm by switching one setting.

**Endpoint** = a target no other target depends on (the final products). **Partial execution** (gwf ≥ 2.2): targets whose inputs are missing *and* not produced by any other target are `skipped` rather than erroring.

---

## Installation & Project Setup

**Conda (recommended — required for the Conda executor and reproducible envs):**
```bash
conda config --add channels gwforg
conda install gwf
# or an isolated project env:
conda create -n myproject gwf dep1 dep2 ...
conda activate myproject
```

**pip:** `pip install gwf`

**Scaffold a project** (interactive — prompts for a backend, writes the files):
```bash
gwf init            # or: gwf init <dir>
```
This generates:

| File | Purpose |
|:-----|:--------|
| `workflow.py` | Defines `gwf = Workflow()` and your targets |
| `templates.py` | Reusable template functions (imported by `workflow.py`) |
| `.gwfconf.json` | Per-project config (e.g. chosen backend) — see [Configuration](#configuration) |

If you run any gwf command with no workflow file present, gwf offers to run `init` for you.

---

## Minimal Workflow

```python
from gwf import Workflow

gwf = Workflow()

gwf.target('MyTarget', inputs=[], outputs=['greeting.txt']) << """
echo hello world > greeting.txt
"""
```

The `<<` operator assigns the target's bash **spec**. Equivalent longhand:
```python
t = gwf.target('MyTarget', inputs=[], outputs=['greeting.txt'])
t.spec = "echo hello world > greeting.txt"
```

> The workflow object **must be named `gwf`** by default. To use a different name or file, pass `-f path.py:objname` (default is `workflow.py:gwf`).

### File-Inferred Dependencies

```python
gwf.target('TargetA', inputs=[], outputs=['x.txt']) << """
echo "this is x" > x.txt
"""

gwf.target('TargetC', inputs=['x.txt', 'y.txt'], outputs=['z.txt']) << """
cat x.txt y.txt > z.txt
"""
```
`TargetC` automatically depends on `TargetA` because it consumes `x.txt`. `TargetC` runs only after `TargetA` has produced its outputs.

### Named Inputs/Outputs (dict form)

Use dicts to give files labels and reference them downstream:
```python
foo = gwf.target(
    name='foo',
    inputs={'A': ['a1', 'a2'], 'B': 'b'},
    outputs={'C': ['a1b', 'a2b'], 'D': 'd'},
)

bar = gwf.target(
    name='bar',
    inputs=foo.outputs['C'],   # reference foo's labeled outputs
    outputs='result',
)
```

### Per-Target Resources & Defaults

```python
foo = gwf.target('foo', inputs=[...], outputs=[...], cores=8, memory='64gb')
print(foo.options)   # => {'cores': 8, 'memory': '64gb'}

# Defaults applied to every target in the workflow:
gwf = Workflow(defaults={'cores': 8, 'memory': '16gb'})
```
Resource keys are **backend-specific** — see [Backends](#backends). `protect=[...]` marks outputs that `gwf clean` must never delete.

---

## Templates (Reusable Targets)

A **template** is a function that returns an `AnonymousTarget` (a target with no name). Put templates in `templates.py`. Instantiate them with `target_from_template(name, template=...)`.

```python
from gwf import AnonymousTarget

def unzip(inputfile, outputfile):
    """A template for unzipping files."""
    inputs = [inputfile]
    outputs = [outputfile]
    options = {'cores': 1, 'memory': '2g'}
    spec = 'zcat {} > {}'.format(inputfile, outputfile)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
```

```python
gwf.target_from_template(
    name='UnzipGenome',
    template=unzip(inputfile='ponAbe2.fa.gz', outputfile='ponAbe2.fa'),
)
```

Real example with several outputs and named-arg formatting:
```python
def bwa_map(ref_genome, r1, r2, bamfile):
    """Map reads to a reference genome."""
    inputs = [r1, r2,
              '{}.amb'.format(ref_genome),
              '{}.ann'.format(ref_genome),
              '{}.pac'.format(ref_genome)]
    outputs = [bamfile]
    options = {'cores': 16, 'memory': '1g'}
    spec = '''
    bwa mem -t 16 {ref_genome} {r1} {r2} | \\
    samtools sort | \\
    samtools rmdup -s - {bamfile}
    '''.format(ref_genome=ref_genome, r1=r1, r2=r2, bamfile=bamfile)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

gwf.target_from_template(
    name='MapReads',
    template=bwa_map(ref_genome='ponAbe2', r1='R1.fastq.gz', r2='R2.fastq.gz', bamfile='Masala.bam'),
)
```

---

## Mapping Over Inputs

`Workflow.map(template_func, inputs, extra=None, name=None, **kwargs)` applies a template across many inputs and returns a `TargetList`.

```python
def transform_photo(path):
    inputs = {'path': path}
    outputs = {'path': path + '.new'}
    return AnonymousTarget(inputs=inputs, outputs=outputs, options={},
                           spec='./transform_photo {}'.format(path))

photos = gwf.glob('photos/*.jpg')   # glob relative to the workflow dir
gwf.map(transform_photo, photos)     # => targets transform_photo_0, _1, _2, ...
```

| Need | How |
|:-----|:----|
| Custom name prefix | `gwf.map(fn, photos, name='TransformPhoto')` → `TransformPhoto_0`, `_1`, … |
| Name from a function | pass `name=f(idx, target)` returning a string (see below) |
| Same extra arg for all | `gwf.map(fn, photos, extra={'width': 800})` |
| Different args per item | pass a list of dicts as `inputs` (keys = template params) |
| Chain stages | feed `.outputs` of one map into the next |

```python
def get_photo_name(idx, target):
    import os.path
    filename = os.path.splitext(os.path.basename(target.inputs['path']))[0]
    return 'transform_photo_{}'.format(filename)

gwf.map(transform_photo, photos, name=get_photo_name)   # transform_photo_dog, _cat, ...

# per-item args:
photos = [{'path': 'photos/dog.jpg', 'width': 600},
          {'path': 'photos/cat.jpg', 'width': 1000}]
gwf.map(transform_photo, photos, name=get_photo_name)

# chaining via .outputs:
transformed = gwf.map(transform_photo, photos, name=get_photo_name)
compressed  = gwf.map(compress_photo, transformed.outputs)
```

**Collecting** many outputs into one downstream target — `collect()` turns a list of dicts into one dict with **pluralized** keys:
```python
from gwf.workflow import collect

bundle = collect(compressed.outputs, ['path'])   # => {'paths': [...]}
gwf.target_from_template(
    name='zip_photos',
    template=zip_files(paths=bundle['paths'], output_path='photos.zip'),
)
```

---

## CLI Command Reference

**Global options come *before* the subcommand:**

| Option | Meaning | Default |
|:-------|:--------|:--------|
| `-f, --file PATH:OBJ` | Workflow file and object to load | `workflow.py:gwf` |
| `-b, --backend NAME` | Backend to use (`local`/`slurm`/`sge`/`lsf`/`pbs`) | from config / guessed |
| `-v, --verbose LEVEL` | `warning` / `info` / `debug` / `error` | `info` |
| `--no-color` / `--use-color` | Toggle colored output | auto |

Get help for any command with `gwf <command> --help`.

| Command | What it does | Key flags |
|:--------|:-------------|:----------|
| `gwf init [DIR]` | Scaffold `workflow.py`, `templates.py`, `.gwfconf.json` (interactive) | — |
| `gwf run [TARGETS...]` | Submit out-of-date targets (and their deps). No args = whole workflow. | `-d/--dry-run`, `-f/--force` (resubmit even if up to date), `-n/--no-deps` (named targets only, ignore deps), `-s/--status STATE` (repeatable), `-g/--group PAT` (repeatable) |
| `gwf status [TARGETS...]` | Show target states | `--endpoints`, `-f/--format {default,grouped,summary}`, `-s/--status STATE` (repeatable), `-g/--group PAT` |
| `gwf logs TARGET` | Show stdout (or stderr) of a target's last run | `-e/--stderr`, `--no-pager` |
| `gwf info [TARGETS...]` | Dump target metadata (inputs/outputs/spec/deps) | `-f/--format {json,pretty}` (default `json`) |
| `gwf cancel [TARGETS...]` | Cancel submitted/running targets. No args = **all** (prompts). | `-f/--force` (no confirm) |
| `gwf clean [TARGETS...]` | Delete output files of **non-endpoint** targets (intermediates) | `--all` (include endpoints), `-f/--force` (no confirm) |
| `gwf touch [TARGETS...]` | Bottom-up `touch` outputs + spec hashes so the workflow looks done | `-c/--create-missing` (also create missing files) |
| `gwf workers` | Start the **local backend** worker pool (see below) | `-n/--num-workers N`, `-p/--port P`, `-h/--host H` |
| `gwf config get/set/unset KEY [VALUE]` | Read/write `.gwfconf.json` | — |

`gwf clean` reports how much data it will delete and **never deletes `protect`ed files**. `gwf touch` is the rescue for accidentally deleted files — it makes gwf believe targets ran without re-running them.

### Filtering & glob patterns
`TARGETS` accept shell-glob name patterns; `--status`/`-s` is repeatable (OR semantics):
```bash
gwf status 'Align*'                          # name glob
gwf status -s running                         # one state
gwf status -s shouldrun -s completed          # either state
gwf status --endpoints --status running 'Map*'
gwf run 'Filter*' -n                          # run matching targets only, no deps
gwf run -f MapReads                           # force re-run a single target
```

---

## Target Status States

Shown by `gwf status` (symbol/color from the source):

| State | Symbol | Meaning |
|:------|:------:|:--------|
| `shouldrun` | ⨯ (magenta) | Out of date — will run next `gwf run` (output missing, input newer, or **spec changed**) |
| `submitted` | – (cyan) | Queued in the backend, not yet started |
| `running` | ↻ (blue) | Currently executing |
| `completed` | ✓ (green) | Outputs exist and are up to date — nothing to do |
| `failed` | ⨯ (red) | The job failed |
| `cancelled` | ⨯ (red) | The job was cancelled |
| `skipped` | ⨯ (red) | Required input missing and not produced by any target (partial execution) |

---

## Backends

Set the active backend per-invocation with `-b NAME`, or persistently with `gwf config set backend NAME`. Built-ins: `local`, `slurm`, `sge`, `lsf`, `pbs`.

### Slurm (the cluster target)

```bash
gwf config set backend slurm     # then plain `gwf run` / `gwf status` use Slurm
# or per command: gwf -b slurm run
```

**Target options (resource specs)** — set as keyword args to `target()` or keys in a template's `options` dict:

| Option | Maps to | Default | Notes |
|:-------|:--------|:--------|:------|
| `cores` (int) | `--cpus-per-task` | `1` | |
| `memory` (str) | `--mem` | `1` | Use explicit units, e.g. `'8g'`, `'64gb'` |
| `walltime` (str) | `--time` | `'01:00:00'` | `HH:MM:SS` (Slurm day form `D-HH:MM:SS` also valid) |
| `queue` (str) | `--partition` | — | Comma-separated list allowed |
| `account` (str) | `--account` | — | |
| `constraint` (str) | `--constraint` | — | |
| `qos` (str) | `--qos` | — | |
| `gres` (str) | `--gres` | — | Commonly GPUs, e.g. `gres='gpu:1'` |
| `mail_type` (str) | `--mail-type` | — | |
| `mail_user` (str) | `--mail-user` | — | |

**Backend config keys:**

| Key | Values | Default |
|:----|:-------|:--------|
| `backend.slurm.log_mode` | `full` (separate stdout/stderr) · `merged` · `none` | `full` |
| `backend.slurm.accounting_enabled` | `true` uses `sacct` for status | `true` |

Example target:
```python
gwf.target('Heavy', inputs=['in.bam'], outputs=['out.vcf'],
           cores=16, memory='32g', walltime='12:00:00',
           queue='normal', account='myproject', gres='gpu:1') << """
...
"""
```

### Local — and how to test on macOS

The **local backend** runs targets on a pool of worker processes you start yourself — ideal for developing/testing the *same* workflow on a Mac before submitting to Slurm. It needs two terminals: one running the worker daemon, one issuing gwf commands.

```bash
# Terminal 1 — start the worker pool (runs in the foreground; Ctrl-C to stop)
gwf -b local workers -n 4

# Terminal 2 — point gwf at the local backend and drive the workflow
gwf -b local run
gwf -b local status
```

Persist the backend so you can drop `-b local`:
```bash
gwf config set backend local
gwf workers -n 4        # terminal 1
gwf run                 # terminal 2
gwf status
```

**`gwf workers` flags / config:**

| Flag | Config key | Default |
|:-----|:-----------|:--------|
| `-n, --num-workers` | — | number of CPU cores |
| `-p, --port` | `local.port` | `12345` |
| `-h, --host` | `local.host` | `localhost` |

The local backend has **no target resource options** (`cores`/`memory`/etc. are accepted but ignored) — concurrency is just the number of workers.

**macOS gotchas when testing locally:**
- Workers execute each spec with **bash in the environment where `gwf workers` was started**. Activate your conda env / put tools on `PATH` *before* launching the worker pool. **Executors (Conda/Singularity/…) only take effect on the Slurm backend** — on local they are ignored, so the env must already be active.
- **Spec portability:** macOS ships BSD userland, the cluster ships GNU. Commands differ — e.g. macOS `gzcat` vs Linux `zcat`, `sed -i ''` vs `sed -i`, no `realpath`/`readlink -f` by default. Keep specs portable or install GNU coreutils (`brew install coreutils gnu-sed`).
- If port `12345` is busy (e.g. a stale worker), start the pool elsewhere and point the backend at it: `gwf workers -p 12346` plus `gwf config set local.port 12346`, so later `gwf run`/`gwf status` connect to the right port.
- Switching to the cluster is a one-liner — the workflow file does not change: `gwf config set backend slurm` (then `gwf run`).

### Other backends (brief)

| Backend | Requires | Target options | Notes |
|:--------|:---------|:---------------|:------|
| `sge` (Sun Grid Engine) | `smp` parallel env (`qconf -spl`) | `cores`, `memory`, `walltime`, `queue`, `account` | |
| `lsf` (IBM Spectrum LSF) | `bsub`, `bjobs` | `cores`, `memory` (def `4GB`), `queue` (def `normal`) | |
| `pbs` (Portable Batch System) | `qsub`, `qstat` | `cores`, `memory` (def `4GB`), `walltime` (`HH:MM:SS`), `queue` (def `normal`), `account` | |

---

## Executors (run targets inside an environment)

Executors wrap each target's spec to run inside a managed environment. **Currently Slurm-only** (other backends planned) — on the local backend they are ignored.

```python
from gwf import Workflow
from gwf.executors import Conda

# default for every target:
gwf = Workflow(executor=Conda("myenv"))

# or per-target:
gwf.target("Test", inputs=[], outputs=[], executor=Conda("myenv")) << """
echo runs inside the myenv conda environment
"""
```

| Executor | Constructor | Requirements |
|:---------|:------------|:-------------|
| `Bash` | `Bash()` (default) | none |
| `Conda` | `Conda(env, debug_mode=False)` — `env` = name or path | `conda` on `PATH` or `CONDA_EXE` set. Env must **already exist** (executor won't create it) |
| `Pixi` | `Pixi(project=None, env='default', debug_mode=False)` | `pixi` on `PATH` or `PIXI_EXE` set |
| `Singularity` | `Singularity(image, flags=None, debug_mode=False)` — `image` = `.sif` | `singularity` on `PATH` |
| `Apptainer` | `Apptainer(image, flags=None, debug_mode=False)` — `image` = `.sif` | `apptainer` on `PATH` |

---

## Configuration

All config is **per-project**, stored in `.gwfconf.json` in the workflow directory (aids reproducibility). Manage with `gwf config`:

```bash
gwf config get backend
gwf config set backend slurm
gwf config set verbose warning
gwf config unset backend
```

Dotted keys for backend-specific settings (e.g. `local.port`, `backend.slurm.log_mode`).

| Core key | Meaning | Default |
|:---------|:--------|:--------|
| `backend` | Active backend (= `--backend`) | `local` |
| `verbose` | Verbosity (= `--verbose`) | `info` |
| `no_color` | Disable colors if `true` | `false` |

Backends and plugins may register additional keys (see their tables above).

---

## Patterns & Tips

**Parameter sweep** — `itertools.product` over a grid:
```python
import itertools
xs, ys, zs = [0,1,2,4,5], ['cold','warm'], [0.1,0.2,0.3]
for x, y, z in itertools.product(xs, ys, zs):
    gwf.target(name='sim_{}_{}_{}'.format(x,y,z),
               inputs=['input.txt'],
               outputs=['output_{}_{}_{}.txt'.format(x,y,z)]) << """
    ./simulate {} {} {}
    """.format(x, y, z)
```

**Dynamic / reusable workflows** — a function that builds and returns a `Workflow` (put in `fancy.py`, import in `workflow.py`):
```python
import os.path
from gwf import Workflow

def my_fancy_workflow(output_dir='outputs/', summarize=True):
    w = Workflow()
    foo_output = os.path.join(output_dir, 'output1.txt')
    w.target('Foo', inputs=['input.txt'], outputs=[foo_output]) << """
    ./run_foo > {}
    """.format(foo_output)
    if summarize:                      # optional part of the pipeline
        bar_output = os.path.join(output_dir, 'output2.txt')
        w.target('Bar', inputs=[foo_output], outputs=[bar_output])
    return w
```
```python
# workflow.py
from fancy import my_fancy_workflow
gwf = my_fancy_workflow(output_dir='new_outputs/', summarize=False)
```

**External config (no code edits)** — read parameters from JSON:
```python
import json
from fancy import my_fancy_workflow
config = json.load(open('config.json'))
gwf = my_fancy_workflow(output_dir=config['output_dir'], summarize=config['summarize'])
```

**Very large workflows** (> ~50k targets) get slow because gwf hits the filesystem per file when scheduling. Split into **named sub-workflows** and run them individually:
```python
from gwf import Workflow
for sample in ['Sample1', 'Sample2', 'Sample3']:
    name = 'Analyse.{}'.format(sample)
    wf = Workflow(name=name)
    wf.target('{}.Filter'.format(sample), inputs=[sample], outputs=[...]) << "..."
    globals()[name] = wf            # expose for `-f workflow.py:Analyse.Sample1`
```
```bash
gwf -f workflow.py:Analyse.Sample1 run   # only schedules Sample1's targets
```

---

## Python API (quick reference)

```python
Workflow(working_dir=…, defaults=…, executor=…)
  .target(name, inputs, outputs, protect=None, group=None, executor=None, **options) -> Target
  .target_from_template(name, template, **options) -> Target   # template is an AnonymousTarget
  .map(template_func, inputs, extra=None, name=None, **kwargs)  -> TargetList
  .glob(pathname, ...) / .iglob(pathname, ...)   # glob relative to working_dir
  .shell(*args, **kwargs)                        # run a shell cmd while building the workflow
  classmethod from_path(path) / from_context(ctx)

AnonymousTarget(inputs, outputs, options, working_dir='.', protect=…, executor=None, spec='')
Target(name, inputs, outputs, options, ...)      # named; name must be a valid Python identifier
TargetList                                       # returned by .map(); has .inputs and .outputs
gwf.workflow.collect(list_of_dicts, keys)        # -> dict with pluralized keys (for fan-in)
```

- `inputs`/`outputs` may be a **str**, **list**, or **dict** (named). Dict form lets you reference `target.outputs['LABEL']` downstream.
- `options` is the backend-specific resource dict (`cores`, `memory`, `walltime`, …). Per-target options override `Workflow(defaults=…)`.

---

## Known Limitations & Sharp Edges

- **Executors are Slurm-only** right now — Conda/Singularity/Pixi/Apptainer have no effect on the local backend. For local testing, activate the environment before `gwf workers`.
- **Conda executor does not create envs** — the named env must already exist.
- **Spec changes trigger re-runs.** gwf hashes each spec; editing spec text (even whitespace) marks a target `shouldrun`. Use `gwf touch` to accept current outputs without re-running.
- **Dependencies are filename-based.** Two targets must agree on the *exact* output/input path string for an edge to form; mismatched relative paths silently break the link.
- **`FileProvidedByMultipleTargetsError`** if two targets declare the same output file; **`CircularDependencyError`** on cycles.
- **`gwf cancel`/`gwf clean` with no targets are sweeping** (all targets / all intermediates). They prompt unless `-f`; `clean` always spares `protect`ed files.
- **Reproducibility:** pin tool versions in a conda `environment.yml`, keep `.gwfconf.json` in the repo, and rely on the file-declared graph so the pipeline rebuilds deterministically across machines and backends.
