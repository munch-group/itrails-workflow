# %% [markdown]
# ---
# title: GWF workflow
# execute:
#   eval: false
# ---

# %% [markdown]
r"""

iTRAILS workflow: fits the TRAILS coalescent HMM (Rivas-González et al., 2024)
to a four-genome MAF alignment and decodes gene-tree topologies along the
genome. The heavy lifting lives in `itrails_workflow.py`, which a parent
project can import when this repository is used as a git submodule (see the
README). This file is the standalone entry point: `gwf run` in this directory
builds the workflow from the analyses listed in `analyses.yml`.

For each analysis, three targets are chained:

```plaintext
             alignment.maf      config.yaml
                    \              /
                  itrails_optimize_{name}
                          |
        optimize/{name}.best_model.yaml  (+ optimization_history.csv,
                   /           \            starting_params.yaml)
                  /             \
itrails_viterbi_{name}       itrails_posterior_{name}
         |                            |
viterbi/{name}.viterbi.csv   posterior/{name}.posterior.csv
(+ hidden_states.csv)        (+ hidden_states.csv)
```

All files for an analysis are written to `{output_dir}/{name}/`, with one
subfolder per step so outputs are easy to separate: `split/`, `optimize/`,
`viterbi/`, `posterior/`, `concat/`.

With a `window_size` key in the analysis, the alignment is instead split
into fixed-size windows at MAF-block boundaries (`itrails_split_{name}`,
writing to `split/`), each window is decoded as its own job, and the
per-window CSVs are concatenated into genome-wide `concat/{name}.viterbi.csv`
/ `concat/{name}.posterior.csv` (`itrails_concat_*_{name}`). The model is
fitted once for the whole alignment (`fit: genome`, default) or per window
(`fit: window`).

An analysis can also start from phased (g)VCFs instead of a MAF (`vcf` +
`samples` keys, optionally `fasta`): an `alignment/` step
(`itrails_vcf2maf_{name}`) reconstructs one haplotype per sample in the
shared reference coordinates — for closely related species all mapped to
the same outgroup genome. Without a `fasta`, all-sites genome VCFs are
converted from the records alone (uncalled positions become N). See the
README and the `itrails_workflow` docstring.

See `claude-itrails-ref.md` for the iTRAILS CLI, config format, and output
file details, and `claude-gwf-ref.md` for gwf itself.
"""

# %% [markdown]
"""
## Imports
"""

# %%
import glob
import os
from pathlib import Path

from gwf import AnonymousTarget, Workflow

from global_params import load_params
from itrails_workflow import itrails_workflow

# %% [markdown]
"""
## Notebook template

Executes a notebook in place once the workflow outputs it depends on exist.
"""


# %%
# task template function
def run_notebook(path, dependencies, memory='8g', walltime='00:10:00', cores=1):
    """
    Executes a notebook inplace and saves the output.
    """
    # path of output sentinel file
    sentinel = str(Path(path).parent / f'.{Path(path).name}.sentinel')

    # input specification
    inputs = [path] + dependencies
    # output specification mapping a label to each file
    outputs = {'sentinel': sentinel}
    # resource specification
    options = {'memory': memory, 'walltime': walltime, 'cores': cores}

    # commands to run in task (bash script)
    spec = f"""
    jupyter nbconvert --to notebook --execute --inplace {path} && touch {sentinel}
    """
    # return target
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# %% [markdown]
"""
## Workflow

The workflow is instantiated unconditionally so that gwf commands work even
before `analyses.yml` is filled in; targets are only added when it exists.
"""

# %%
gwf = Workflow(working_dir=os.getcwd())

params_file = Path(__file__).parent / 'analyses.yml'
if params_file.exists():
    params = load_params(params_file)

    gwf, targets = itrails_workflow(
        gwf=gwf,
        analyses=params.analyses,
        output_dir=params.output_dir,
        account=params.account,
    )

    # make notebooks depend on all output files from the workflow
    notebook_dependencies = []
    for target in gwf.targets.values():
        outputs = target.outputs
        if type(outputs) is dict:
            outputs = outputs.values()
        for output in outputs:
            # dict values may themselves be lists (e.g. the window MAFs
            # produced by the split target in windowed mode)
            if isinstance(output, list):
                notebook_dependencies.extend(output)
            else:
                notebook_dependencies.append(output)

    # run notebooks in sorted order nb01_, nb02_, ...
    for path in sorted(glob.glob('notebooks/*.ipynb')):
        target = gwf.target_from_template(
            os.path.basename(path).replace('.', '_'),
            run_notebook(path, notebook_dependencies))
        # make each notebook depend on all previous notebooks
        notebook_dependencies.append(target.outputs['sentinel'])

# %%
