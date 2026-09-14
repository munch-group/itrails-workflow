# GWF (Grid Workflow) Reference

Local reference compiled from <https://gwf.app/> (gwf 2.1.1). GWF is a Python workflow manager for building and running scientific pipelines on Slurm, SGE, LSF, PBS, or locally. See `workflow.py` in this project for a worked example that uses most of the patterns below.

---

## Table of contents

1. [Core concepts](#core-concepts)
2. [Minimal workflow](#minimal-workflow)
3. [API: `Workflow`](#api-workflow)
4. [API: `AnonymousTarget` / `Target` / `TargetList`](#api-anonymoustarget--target--targetlist)
5. [API: `collect`](#api-collect)
6. [Dependency resolution](#dependency-resolution)
7. [Inputs / outputs (strings, lists, dicts)](#inputs--outputs-strings-lists-dicts)
8. [Resources & target options](#resources--target-options)
9. [Backends](#backends)
10. [Executors (Slurm only)](#executors-slurm-only)
11. [CLI reference](#cli-reference)
12. [Configuration (`.gwfconf.json`)](#configuration-gwfconfjson)
13. [Target states and status output](#target-states-and-status-output)
14. [Failure handling & the write-then-move idiom](#failure-handling--the-write-then-move-idiom)
15. [Partial workflow execution](#partial-workflow-execution)
16. [Cleaning & protecting outputs](#cleaning--protecting-outputs)
17. [Patterns](#patterns)
18. [Internal APIs (rarely needed)](#internal-apis-rarely-needed)
19. [Gotchas](#gotchas)

---

## Core concepts

- **Target** — a single unit of computation consuming zero or more input files and producing zero or more output files. Either named (via `gwf.target(...)`) or anonymous (via `AnonymousTarget(...)`, returned from a template function). Targets only (re)run when outputs are missing or older than inputs, or when the spec has changed.
- **Template function** — a plain Python function returning an `AnonymousTarget`. The reusable building block. Parameterize it, then instantiate many targets from it with `gwf.target_from_template` (one-off) or `gwf.map` (fan-out).
- **Dependency resolution is filename-based.** GWF matches a target's declared `inputs` against every other target's declared `outputs`; shared paths become edges. There is no explicit `depends_on=`.
- **Backend** — where targets run: `slurm`, `sge`, `lsf`, `pbs`, `local`. Chosen per project via `gwf config set backend <name>`; stored in `.gwfconf.json`.
- **Executor** — wraps the target spec in a runtime environment: `Bash` (default), `Conda`, `Pixi`, `Singularity`, `Apptainer`. **Slurm backend only.**
- **Endpoint** — a target that no other target depends on. `gwf clean` (without `--all`) preserves endpoint outputs; `gwf status --endpoints` shows only these.

---

## Minimal workflow

```python
# workflow.py
from gwf import Workflow

gwf = Workflow()

gwf.target('MyTarget', inputs=[], outputs=['greeting.txt']) << """
echo hello world > greeting.txt
"""
```

Run:

```bash
gwf workers &       # local backend only — start a worker pool
gwf run             # schedule & submit
gwf status          # progress
```

On re-run, `MyTarget` is skipped because `greeting.txt` exists and is newer than its (empty) inputs.

The `<<` operator assigns the trailing string as the target's Bash spec. This is pure Python — `gwf.target(...)` returns a `Target` object and `__lshift__` is implemented to set its spec.

---

## API: `Workflow`

```python
class gwf.Workflow(working_dir=NOTHING, defaults=NOTHING, executor=NOTHING)
```

- `working_dir` (str): directory containing the workflow file; every relative path used in targets is resolved against this.
- `defaults` (dict): per-target option defaults, overridden by options on individual targets. E.g. `{'account': 'myproj', 'cores': 8, 'memory': '16g', 'walltime': '04:00:00'}`.
- `executor` (Executor): default executor used for every target (Slurm backend only).

### `gwf.target(name, inputs, outputs, protect=None, group=None, executor=None, **options)`

Creates a named target and adds it to the workflow. Returns a `Target` instance.

- `name` (str): must be a valid Python identifier (letters, digits, underscores; no leading digit).
- `inputs`, `outputs`: str | list | dict (see [Inputs / outputs](#inputs--outputs-strings-lists-dicts)).
- `protect`: iterable of output paths `gwf clean` must not delete.
- `group` (str): optional group label (used by some plugins).
- `executor` (Executor): per-target override.
- `**options`: backend options — `cores`, `memory`, `walltime`, `account`, `queue`, `gres`, `constraint`, `qos`, `mail_type`, `mail_user`, ...

Use `<<` to assign a Bash spec:

```python
gwf.target('Unzip', inputs=['a.gz'], outputs=['a'], cores=1, memory='2g') << """
gzcat a.gz > a
"""
```

### `gwf.target_from_template(name, template, **options)`

Creates a target from an `AnonymousTarget` returned by a template function. Extra `**options` passed here override those baked into the template.

```python
gwf.target_from_template(
    'MergeMe',
    merge_names(collected['unique_me_paths'], 'results/merged.txt'),
    memory='32g',           # overrides the template's memory setting
)
```

### `gwf.map(template_func, inputs, extra=None, name=None, **kwargs)`

Fan-out helper. Calls `template_func(...)` once per item in `inputs`, registers each resulting target, returns a `TargetList` whose `.outputs` / `.inputs` properties expose all generated outputs/inputs as lists of dicts.

- `inputs`: iterable of
  - **scalars** (passed as the template's first positional arg), or
  - **dicts** (spread as `**kwargs` into the template — keys must match template parameter names), or
  - another map's `.outputs` (each dict in that list is spread as `**kwargs`).
- `extra` (mapping): constant kwargs added to every call, e.g. `extra=dict(me='Kasper')`.
- `name`:
  - `None` → `"<template_func_name>_0"`, `"_1"`, ... (default).
  - `str` → prefix: `"Foo_0"`, `"Foo_1"`, ...
  - `callable(idx, target) -> str` → full control; `target` is the `AnonymousTarget` built for that iteration, letting you name from its inputs.
- `**kwargs`: option overrides applied to every generated target.

Chaining:

```python
uppercased = gwf.map(uppercase_names, ['a.txt', 'b.txt'])
# uppercased.outputs == [{'uppercased_path': '...'}, {'uppercased_path': '...'}]
filtered   = gwf.map(divide_names, uppercased.outputs, extra=dict(me='Kasper'))
unique     = gwf.map(unique_names, filtered.outputs)
```

Custom naming from input paths:

```python
import os.path
def photo_name(idx, target):
    stem = os.path.splitext(os.path.basename(target.inputs['path']))[0]
    return f'transform_{stem}'

gwf.map(transform_photo, photos, name=photo_name)
```

### `gwf.glob(pathname, *args, **kwargs)` / `gwf.iglob(pathname, *args, **kwargs)`

Equivalent to `glob.glob` / `glob.iglob`, but relative patterns are resolved against the workflow's `working_dir`.

```python
photos = gwf.glob('photos/*.jpg')
```

### `gwf.shell(*args, **kwargs)`

Equivalent to `subprocess.check_output`, but run in the workflow's `working_dir`. Returns bytes — pass `universal_newlines=True` for `str`. Rarely needed in workflow files.

### Classmethods

- `Workflow.from_context(ctx)` — build from an internal gwf context object (used by the CLI).
- `Workflow.from_path(path)` — load a workflow from `path/workflow.py`, optionally selecting a named workflow with `path/workflow.py:NAME` (used for the multi-workflow pattern).

---

## API: `AnonymousTarget` / `Target` / `TargetList`

### `AnonymousTarget`

```python
class gwf.AnonymousTarget(
    inputs, outputs, options,
    group=None, working_dir='.', protect=NOTHING,
    executor=None, spec='',
)
```

Unnamed target — what a template function returns.

- `inputs`, `outputs`: str | list | dict. Dict keys become labels accessible via `target.outputs['key']` downstream.
- `options` (dict): backend options. Keys with value `None` are stripped before submission.
- `group` (str): optional group label.
- `working_dir` (str): defaults to `'.'`; set if the target must run in a different directory than the workflow root.
- `protect`: iterable of output paths excluded from `gwf clean`.
- `executor` (Executor | None): per-target runtime wrapper.
- `spec` (str): Bash script.

### `Target`

```python
class gwf.Target(
    name, inputs, outputs, options,
    group=None, working_dir='.', protect=NOTHING,
    executor=NOTHING, spec='',
)
```

Named target (subclass of `AnonymousTarget`). You rarely construct this directly — `gwf.target(...)` and `gwf.target_from_template(...)` do it for you. Name must be a valid Python identifier.

### `TargetList`

Returned by `gwf.map`. Has:

- `.inputs` — list of per-target input structures (strings/lists/dicts).
- `.outputs` — list of per-target output structures (strings/lists/dicts).

You feed `.outputs` into the next `gwf.map` call to chain stages.

---

## API: `collect`

```python
from gwf.workflow import collect
collect(outputs_list, keys)
```

Transforms a list-of-dicts into a dict-of-lists, pluralizing each key by appending `s`.

```python
# input
[
    {'path': 'a.jpg'},
    {'path': 'b.jpg'},
    {'path': 'c.jpg'},
]
collect(_, ['path'])
# -> {'paths': ['a.jpg', 'b.jpg', 'c.jpg']}
```

Typical use — fan-in after a map:

```python
collected = collect(unique_names_targets.outputs, ['unique_me_path'])
# -> {'unique_me_paths': ['...', '...', ...]}

gwf.target_from_template(
    'MergeMe',
    merge_names(collected['unique_me_paths'], 'results/merged.txt'),
)
```

Pluralization rule: `key` → `key + 's'`. Singular template labels (`path`, `sentinel`) work cleanly; plural labels (`paths`) yield awkward `pathss` — pick singular names.

---

## Dependency resolution

There is no `depends_on`. Edges are derived from paths:

1. Every target declares `inputs` (paths read) and `outputs` (paths written).
2. GWF indexes every output path → producing target.
3. For each target, each input path is looked up in that index; a match creates an edge.
4. Input paths not produced by any target are **unresolved** — they must already exist on disk. If they don't, that target and everything downstream is skipped (since gwf 2.2 — see [Partial execution](#partial-workflow-execution)).

**Errors**:

- Same path produced by two targets → `FileProvidedByMultipleTargetsError` at graph construction.
- Target listed as its own input → `CircularDependencyError` ("Target X depends on itself").

Verify what gwf thinks about your graph with `gwf -v debug status`.

---

## Inputs / outputs (strings, lists, dicts)

All three shapes are accepted:

```python
# string — single file
gwf.target('t', inputs='a.txt', outputs='b.txt') << "cp a.txt b.txt"

# list — multiple files, positional
gwf.target('t', inputs=['a', 'b'], outputs=['c']) << "cat a b > c"

# dict — labelled groups, accessible downstream
foo = gwf.target(
    name='foo',
    inputs={'A': ['a1', 'a2'], 'B': 'b'},
    outputs={'C': ['a1b', 'a2b'], 'D': 'd'},
)

bar = gwf.target(
    name='bar',
    inputs=foo.outputs['C'],   # reference by label
    outputs='result',
)
```

Dict form is what enables `collect` and labelled `target_from_template` output referencing. Use it when a target produces multiple logically distinct output files.

---

## Resources & target options

Per-target (inside `AnonymousTarget.options` or as kwargs on `gwf.target`):

| Option | Example | Backends | Notes |
|---|---|---|---|
| `cores` | `8` | all except local | int; default `1` |
| `memory` | `'16g'`, `'64gb'` | all except local | units: `k`, `m`, `g`, `gb`; default `1` (Slurm/SGE), `4GB` (LSF/PBS) |
| `walltime` | `'04:00:00'` | Slurm, SGE | `HH:MM:SS`; default `01:00:00` |
| `queue` | `'normal'` | Slurm, SGE, LSF, PBS | partition/queue; comma-separated for multiple (Slurm/SGE) |
| `account` | `'myproj'` | Slurm, SGE | Slurm account; SGE project |
| `gres` | `'gpu:1'` | Slurm | GPU / generic resources |
| `constraint` | `'haswell'` | Slurm | `--constraint` |
| `qos` | `'high'` | Slurm | `--qos` |
| `mail_type` | `'FAIL'` | Slurm | `--mail-type` |
| `mail_user` | `'me@x'` | Slurm | `--mail-user` |

Set workflow-wide defaults in the `Workflow` constructor:

```python
gwf = Workflow(defaults={'account': 'myproj', 'cores': 4, 'memory': '8g'})
```

A target's own options override workflow defaults. `target_from_template(..., **options)` overrides both.

---

## Backends

Switch with `gwf config set backend <name>` or the `-b <name>` CLI flag.

### `local`

Runs targets on the machine running `gwf`, via a worker pool. **You must start workers separately.**

```bash
gwf -b local workers -n 4     # start 4 workers
gwf -b local run              # in another terminal
```

Backend options:

- `local.host` (str, default `localhost`)
- `local.port` (int, default `12345`)

No target options.

### `slurm`

Submits via `sbatch`, polls via `squeue` / `sacct`.

Backend options:

- `backend.slurm.log_mode` — `full` (default, separate stdout/stderr), `merged`, or `none`.
- `backend.slurm.accounting_enabled` — use `sacct` for post-completion status (default `true`).

Supports the full resource option set listed above plus executors.

### `sge`

Sun Grid Engine. Requires a parallel environment named `smp` — check with `qconf -spl`. No backend-level options. Supports `cores`, `memory`, `walltime`, `queue`, `account` (SGE project).

### `lsf`

IBM Spectrum LSF. Requires `bsub`/`bjobs`. No backend-level options. Supports `cores` (default 1), `memory` (default `4GB`), `queue` (default `normal`).

### `pbs`

Portable Batch System. Requires `qsub`/`qstat`. No backend-level options. Same target options as LSF.

---

## Executors (Slurm only)

Wrap every target's spec in an activation step. As of 2.1.1, **Slurm backend only**. Pass via `Workflow(executor=...)` for a default, or per-target via `gwf.target(..., executor=...)`.

```python
from gwf import Workflow
from gwf.executors import Bash, Conda, Pixi, Singularity, Apptainer

gwf = Workflow(executor=Conda('myenv'))
```

### `Bash()`

The default. No activation, spec runs directly in Bash.

### `Conda(env, debug_mode=False)`

- `env` (str): environment name or absolute path to an env dir.
- `debug_mode` (bool): print activation details to the log.

Conda env must already exist — the executor will not create it.

### `Pixi(project=None, env='default', debug_mode=False)`

- `project` (str | None): path to Pixi project; defaults to workflow `working_dir`.
- `env` (str): environment name within the project.
- `debug_mode` (bool).

Env must exist (run `pixi install` first).

### `Singularity(image, flags=(), debug_mode=False)` / `Apptainer(image, flags=(), debug_mode=False)`

- `image` (str): path to `.sif` (or Apptainer image).
- `flags` (iterable[str]): extra args for `singularity exec` / `apptainer exec`, e.g. bind mounts.
- `debug_mode` (bool).

Per-target override example:

```python
gwf.target('Test', inputs=[], outputs=[], executor=Conda('myenv')) << """
echo hello
"""
```

---

## CLI reference

| Command | What it does |
|---|---|
| `gwf run` | Schedule & submit all out-of-date targets |
| `gwf run TARGET` | Run a specific target and its upstream dependencies |
| `gwf run -f path/wf.py:NAME` | Use a named workflow inside a file with multiple |
| `gwf status` | Show state of every target |
| `gwf status --endpoints` | Only leaf targets |
| `gwf status -s STATE` | Filter: `shouldrun`, `submitted`, `running`, `completed`, `failed`, `cancelled` |
| `gwf status 'Foo*'` | Glob pattern on target names |
| `gwf logs TARGET` | stdout from the most recent run |
| `gwf logs --stderr TARGET` (`-e`) | stderr |
| `gwf cancel TARGET` | Cancel a queued/running target (supports globs) |
| `gwf clean` | Delete non-endpoint (intermediate) outputs |
| `gwf clean --all` | Delete all outputs including endpoints |
| `gwf workers [-n N]` | Start N local workers (local backend only) |
| `gwf config set KEY VALUE` | Write setting to `.gwfconf.json` |
| `gwf config get KEY` | Read setting |
| `gwf info` | Diagnostics about the workflow and backend |
| `gwf --help` | Top-level help |
| `gwf <cmd> --help` | Per-command help |

Global flags:

- `-v debug|info|warning|error` — verbosity
- `-b <backend>` — one-off backend override
- `-f path/workflow.py[:NAME]` — use a non-default workflow file or named workflow

GWF does not archive logs — `gwf logs` always shows the most recent run.

---

## Configuration (`.gwfconf.json`)

Project-scoped JSON, stored next to `workflow.py`. Managed via `gwf config`:

```bash
gwf config set backend slurm
gwf config set verbose warning
gwf config set local.port 4321
gwf config set backend.slurm.log_mode merged
gwf config get backend
```

Core keys:

- `backend` (str, default `local`)
- `verbose` (str, default `info`) — one of `debug`, `info`, `warning`, `error`
- `no_color` (bool, default `false`)

Backend-specific keys use dotted paths, e.g. `backend.slurm.accounting_enabled`, `local.port`.

Commit `.gwfconf.json` if its values are shared across machines (e.g. Slurm defaults); keep it out of VCS if it encodes per-machine ports or paths.

---

## Target states and status output

`gwf status` output format:

```
<glyph> <name>   <percent>%    <reason>
```

| Glyph | State | Meaning |
|---|---|---|
| `✓` | completed | Outputs exist and are up-to-date |
| `↻` | running | Currently executing on the backend |
| `-` | submitted | In the backend's queue |
| `⨯` | shouldrun | Needs to (re)run — outputs missing, newer inputs, or spec changed |

Example reasons shown: `"not scheduled because it is a source"`, `"is up-to-date"`, `"spec has changed"`, `"a dependency was scheduled"`, `"has been submitted"`, `"is running"`.

For a live view of the whole graph during execution, keep `watch gwf status` in a second terminal.

---

## Failure handling & the write-then-move idiom

GWF determines success by checking whether outputs exist and are newer than inputs. It has no concept of exit codes or partial writes. This leads to two failure modes:

1. **No outputs produced.** Target exits non-zero before writing outputs. GWF correctly marks the target as `shouldrun` on the next `gwf run`.
2. **Partial/corrupt outputs produced.** Script writes outputs, then crashes (or outputs are truncated). GWF sees the files exist, marks the target `completed` — wrong. Downstream targets consume the corrupt output.

**Fix: the write-then-move idiom** (used throughout `workflow.py` in this project):

```python
spec = f"""
mkdir -p {output_dir}
<command> > {tmp_output_path} &&
    mv {tmp_output_path} {output_path}
"""
```

- Write to `/tmp` (or any scratch path).
- Move to the real output path **only after success** (`&&` short-circuits on failure).
- If the command crashes, the real output path never appears → GWF reruns the target.

For multi-output targets, write all tmp files, then move all of them under one `&&` chain so either all outputs appear or none do.

---

## Partial workflow execution

Since gwf 2.2, missing unresolved inputs are not a fatal error. Example:

```python
gwf.target("AnalyzeBlood",   inputs=["blood.txt"],              outputs=["levels.txt"])  << "..."
gwf.target("SummarizeLevels",inputs=["levels.txt"],             outputs=["summary.txt"]) << "..."
gwf.target("AnalyzeGenome",  inputs=["genome.fa","reference.fa"],outputs=["mutations.txt"])<< "..."
gwf.target("BuildReport",    inputs=["mutations.txt","levels.txt"],outputs=["report.pdf"]) << "..."
```

If `blood.txt` exists on disk but `genome.fa` does not, `gwf run` will schedule only `AnalyzeBlood` and `SummarizeLevels`. `AnalyzeGenome` and `BuildReport` are silently skipped until their inputs appear.

Use this deliberately when some data arrives later than the rest. Be careful: a typo in an input path also "disappears" as a partial-execution skip — use `gwf -v debug status` to see what was skipped and why.

---

## Cleaning & protecting outputs

```bash
gwf clean         # delete non-endpoint outputs
gwf clean --all   # delete everything declared as an output
```

Mark specific outputs as protected — they survive even `gwf clean --all`:

```python
gwf.target('TargetA', inputs=['a'], outputs=['b', 'c', 'd'], protect=['d']) << "..."
```

`protect` works the same way on `AnonymousTarget`.

Endpoint detection is automatic: any target whose outputs are not listed as another target's inputs is an endpoint.

---

## Patterns

### Parameter sweep

```python
import itertools
for x, y, z in itertools.product(xs, ys, zs):
    gwf.target(
        name=f'sim_{x}_{y}_{z}',
        inputs=['input.txt'],
        outputs=[f'output_{x}_{y}_{z}.txt'],
    ) << f"./simulate {x} {y} {z}"
```

### Map → collect → merge fan-in

From this project's `workflow.py`:

```python
uppercased = gwf.map(uppercase_names, input_files)
filtered   = gwf.map(divide_names, uppercased.outputs, extra=dict(me='Kasper'))
unique     = gwf.map(unique_names, filtered.outputs)

# Fan-in: gather all 'unique_me_path' outputs into one list
collected = collect(unique.outputs, ['unique_me_path'])

gwf.target_from_template(
    'MergeMe',
    merge_names(collected['unique_me_paths'], 'results/merged_me.txt'),
)
```

### Factory function that builds a workflow

```python
def build(output_dir='outputs/', summarize=True):
    w = Workflow()
    # ... add targets ...
    if summarize:
        w.target('Summary', ...)
    return w

gwf = build(output_dir='runs/2026/')
```

Drive from JSON config so you don't edit `workflow.py`:

```python
import json
cfg = json.load(open('config.json'))
gwf = build(**cfg)
```

### One workflow per sample (large workflows)

```python
for sample in ['S1', 'S2', 'S3']:
    name = f'Analyse.{sample}'
    wf = Workflow(name=name)
    wf.target(f'{sample}.Filter',    ...) << "..."
    wf.target(f'{sample}.Summaries', ...) << "..."
    globals()[name] = wf
```

Run just one: `gwf -f workflow.py:Analyse.S1 run`.

### Notebook as a target (sentinel pattern)

Executing a notebook in-place doesn't produce a new file, so use a sentinel:

```python
def run_notebook(path, dependencies, memory='8g', walltime='00:10:00', cores=1):
    sentinel = modify_path(path, base=f'.{Path(path).name}', suffix='.sentinel')
    return AnonymousTarget(
        inputs=[path] + dependencies,
        outputs={'sentinel': sentinel},
        options={'memory': memory, 'walltime': walltime, 'cores': cores},
        spec=f"jupyter nbconvert --to notebook --execute --inplace {path} && touch {sentinel}",
    )
```

Chain notebooks by appending each sentinel to the next notebook's `dependencies` list.

### Per-target Conda / container environment

```python
from gwf.executors import Conda, Singularity

gwf.target('A', inputs=[], outputs=['a'], executor=Conda('bio-env'))          << "run_bio.sh"
gwf.target('B', inputs=['a'], outputs=['b'], executor=Singularity('r.sif'))   << "Rscript x.R"
```

### Dynamic notebook dependency injection (from `workflow.py`)

Make a downstream notebook depend on every output produced by the workflow so far:

```python
notebook_dependencies = []
for t in gwf.targets.values():
    outs = t.outputs
    if isinstance(outs, dict):
        notebook_dependencies.extend(outs.values())
    elif isinstance(outs, list):
        notebook_dependencies.extend(outs)

for path in sorted(glob.glob('notebooks/*.ipynb')):
    target = gwf.target_from_template(
        os.path.basename(path),
        run_notebook(path, notebook_dependencies),
    )
    notebook_dependencies.append(target.outputs['sentinel'])
```

---

## Internal APIs (rarely needed)

These are documented for completeness — most workflows never touch them. Useful if you're writing a plugin or a custom `gwf` wrapper.

### `gwf.core.Graph`

```python
Graph(targets, provides, dependencies, dependents, unresolved)
```

Dependency graph. Attributes:

- `targets: dict[str, Target]`
- `provides: dict[path, Target]` — file → producing target
- `dependencies: dict[Target, set[Target]]`
- `dependents: dict[Target, set[Target]]`
- `unresolved: set[path]` — input paths no target produces

Methods:

- `Graph.from_targets(targets) -> Graph` — build from a target iterable. Raises `FileProvidedByMultipleTargetsError` or `CircularDependencyError`.
- `graph.dfs(root)` — depth-first traversal from a root target.
- `graph.endpoints() -> set[Target]` — targets with no dependents.

### `gwf.scheduling`

- `submit_workflow(endpoints, graph, fs, spec_hashes, backend, dry_run=False, force=False, no_deps=False)` — top-level submit.
- `submit_backend(target, dependencies, backend, spec_hashes)` — submit one target; preferred over `backend.submit()` directly because it applies defaults and strips `None` options.
- `get_status_map(graph, fs, spec_hashes, backend, endpoints=None)` — current state of every target.

### `gwf.backends`

- `BackendStatus` enum: `UNKNOWN=0`, `SUBMITTED=1`, `RUNNING=2`, `COMPLETED=3`, `FAILED=4`, `CANCELLED=5`.
- `create_backend(name, working_dir, config) -> type` — returns the backend **class** (not an instance).
- `discover_backends()` — find installed backend plugins.
- `guess_backend()` — pick one based on environment.
- `list_backends() -> list[str]`.

### `gwf.filtering`

- `filter_names(targets, patterns)` — glob-match target names. `patterns` is a list, e.g. `['Foo*', 'Bar*']`.
- `filter_generic(targets, filters)` — run targets through multiple filter instances; returns those passing all.
- `ApplyMixin` — base class for predicate filters. Subclass and implement `predicate(target) -> bool`; `apply(targets)` yields targets where predicate is true.

---

## Gotchas

- **Target names must be valid Python identifiers.** Dots are tolerated in some contexts (see large-workflows pattern, `Analyse.S1`) but safest to stick to `[A-Za-z_][A-Za-z0-9_]*`.
- **File-based dependencies**: two targets writing the same output path is a hard error (`FileProvidedByMultipleTargetsError`). Keep outputs unique per target.
- **`inputs=` must list everything the spec reads**; otherwise the producing target is not scheduled first, and your spec sees a stale or missing file.
- **Spec changes invalidate caches** — modifying the Bash script forces reruns. GWF stores a hash of each target's spec.
- **Write-then-move, always**, for any target that writes outputs (see [Failure handling](#failure-handling--the-write-then-move-idiom)).
- **Executors are Slurm-only** as of 2.1.1. Local/SGE/LSF/PBS ignore the `executor=` argument.
- **`gwf run` checks file mtimes.** If target B edits target A's output in place, GWF thinks A needs to run again. Don't mutate another target's outputs — write your own.
- **`collect()` pluralizes by appending `s`.** Singular labels (`path`, `sentinel`) give `paths`, `sentinels`. Plural labels (`paths`) give `pathss` — avoid.
- **`gwf logs` shows only the most recent run.** For persistent logs, redirect inside your spec or use a plugin.
- **Unresolved inputs silently skip downstream targets** (since 2.2). A typo in an input path won't error — run `gwf -v debug status` to see what was skipped.
- **Resources aren't enforced locally**: `cores` / `memory` / `walltime` are ignored by the `local` backend. Don't rely on them for concurrency limits in local runs.
- **The `gwf workers` process must keep running** on the local backend; kill it (Ctrl-C) when done, or orphaned targets will keep executing.

Upstream docs: <https://gwf.app/>. Source: <https://github.com/gwforg/gwf>.
