"""
GWF workflow component for iTRAILS.

Fits the TRAILS coalescent HMM (Rivas-Gonzalez et al., 2024) to a four-genome
MAF alignment and decodes gene-tree topologies along the genome. Every step
writes to its own subfolder of `{output_dir}/{name}/`, so outputs are easy to
separate by step. For each analysis, three targets are created:

    itrails-optimize  (MAF + config.yaml -> optimize/{name}.best_model.yaml)
    itrails-viterbi   (best_model + MAF  -> viterbi/{name}.viterbi.csv)
    itrails-posterior (best_model + MAF  -> posterior/{name}.posterior.csv)

With a `window_size`, the alignment is instead split into fixed-size genomic
windows at MAF-block boundaries (lossless: iTRAILS treats blocks as
independent sequences) written to `split/`, each window is decoded separately
(`viterbi/w0000.viterbi.csv`, ...), and the per-window CSVs are concatenated
into genome-wide `concat/{name}.viterbi.csv` / `concat/{name}.posterior.csv`
annotated with window and block coordinates. The model is fitted once on the
whole alignment (`fit='genome'`, the default) or independently per window
(`fit='window'`, e.g. to study parameter variation along the genome).

Instead of a ready-made MAF, an analysis can start from phased (g)VCFs with
`vcf` + `samples` (4 sample IDs, first haplotype each; the reserved name REF
is the reference sequence itself), optionally plus `fasta`. This suits
closely related species all mapped to the same outgroup reference — e.g.
baboon samples against the macaque genome — where shared coordinates make
the genomes implicitly aligned. With a `fasta`, variant alleles are written
onto the reference sequence; without one, sequences come from the VCF
records alone (all-sites genome VCFs; uncalled positions become N). An
`alignment/` step materializes the MAF and everything downstream is
identical.

This module is written so the repository can be used as a git submodule of a
larger gwf project. The parent workflow imports `itrails_workflow` and grafts
the targets onto its own Workflow object:

    # parent workflow.py, with this repo as a submodule in ./itrails/
    import os, sys
    from gwf import Workflow

    gwf = Workflow(defaults={'account': 'my-project'})

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'itrails'))
    from itrails_workflow import itrails_workflow

    gwf, itrails_targets = itrails_workflow(
        gwf=gwf,
        analyses=[dict(name='primates',
                       maf='data/alignment.maf',
                       config='config/itrails_config.yaml',
                       window_size=2_000_000)],
        output_dir='steps/itrails',
    )

The module is deliberately *not* named `workflow.py` (that would collide with
the parent's own workflow module) nor importable as `itrails.workflow` (that
would shadow the installed `itrails` package).

See claude-itrails-ref.md for details on the iTRAILS CLI and file formats.
"""

import importlib.util
import json
import os
import re
from collections import defaultdict
from pathlib import Path

from gwf import AnonymousTarget, Workflow

# Absolute path to the directory holding this file, resolved at import time so
# that specs refer to this component's pixi environment no matter which
# directory the parent workflow is run from.
COMPONENT_DIR = Path(__file__).parent.resolve()

_SCRIPT_MODULES = {}


def _script_module(name):
    """
    A module from scripts/, loaded via importlib rather than sys.path so the
    parent project's namespace is untouched. Sharing the code between the
    targets (run time) and the workflow definition guarantees both derive
    the same windows/contigs from the same inputs.
    """
    if name not in _SCRIPT_MODULES:
        path = COMPONENT_DIR / 'scripts' / f'{name}.py'
        spec = importlib.util.spec_from_file_location(f'_itrails_{name}', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _SCRIPT_MODULES[name] = module
    return _SCRIPT_MODULES[name]


def _split_maf_module():
    return _script_module('split_maf')


def default_run_prefix():
    """
    Command prefix that runs the rest of a command line inside this
    repository's pixi `itrails` environment (works on both the local and
    slurm backends, unlike gwf executors which are slurm-only).
    """
    return f"pixi run --manifest-path {COMPONENT_DIR / 'pixi.toml'} -e itrails"


def _read_fai(fai_path):
    """[(contig, length)] in file order from a samtools faidx index."""
    contigs = []
    with open(fai_path) as f:
        for line in f:
            if line.strip():
                fields = line.split('\t')
                contigs.append((fields[0], int(fields[1])))
    return contigs


def _windows_from_lengths(contigs, block_size, window_size, bin_species):
    """
    Window index for a VCF-derived alignment, computed from contig lengths
    alone: vcf_to_maf writes blocks of `block_size` tiling every contig, so
    the windows the split will produce are known before the MAF exists.
    Mirrors scripts/split_maf.py binning (blocks assigned by start).
    """
    label = _split_maf_module().window_label
    windows = []
    for contig, length in contigs:
        chrom = f'{bin_species}.{contig}'
        bins = {}
        for start in range(0, length, block_size):
            bins[start // window_size] = bins.get(start // window_size, 0) + 1
        for bin_idx in sorted(bins):
            windows.append({'label': label(len(windows)), 'chrom': chrom,
                            'start': bin_idx * window_size,
                            'end': (bin_idx + 1) * window_size,
                            'n_blocks': bins[bin_idx]})
    return windows


def _window_index(maf_path, window_size, reference, species, cache_dir):
    """
    Windows of `maf_path` as [{'label', 'chrom', 'start', 'end', 'n_blocks'}],
    computed with the binning code in scripts/split_maf.py. The scan reads the
    whole alignment, so the result is cached in cache_dir keyed on the file's
    size and mtime — repeated `gwf status` calls don't rescan a large MAF.
    """
    stat = os.stat(maf_path)
    key = {'maf': os.path.abspath(maf_path), 'size': stat.st_size,
           'mtime_ns': stat.st_mtime_ns, 'window_size': window_size,
           'reference': reference, 'species': list(species) if species else None}
    cache_file = os.path.join(cache_dir, '.window_index.json')
    try:
        with open(cache_file) as f:
            cached = json.load(f)
        if cached.get('key') == key:
            return cached['windows']
    except (OSError, ValueError):
        pass

    module = _split_maf_module()
    blocks = module.scan_blocks(maf_path, reference=reference, species=species)
    windows = [{'label': w['label'], 'chrom': w['chrom'], 'start': w['start'],
                'end': w['end'], 'n_blocks': len(w['blocks'])}
               for w in module.assign_windows(blocks, window_size)]
    try:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache_file, 'w') as f:
            json.dump({'key': key, 'windows': windows}, f)
    except OSError:
        pass
    return windows


# task template function
def split_maf(maf_path, split_dir, window_size, window_mafs,
              reference=None, species=None, run_prefix='',
              cores=1, memory='8g', walltime='4:00:00', **extra_options):
    """
    Split a MAF alignment into fixed-size windows at block boundaries with
    scripts/split_maf.py. `window_mafs` is the list of files the split will
    produce, precomputed at workflow-definition time by _window_index with
    the same binning code the script uses.
    """
    inputs = [maf_path]
    outputs = {
        'windows': list(window_mafs),
        'manifest': os.path.join(split_dir, 'windows.tsv'),
        'blocks': os.path.join(split_dir, 'blocks.tsv'),
    }
    options = {'cores': cores, 'memory': memory, 'walltime': walltime,
               **extra_options}
    reference_flag = f'--reference {reference}' if reference else ''
    species_flag = f"--species {' '.join(species)}" if species else ''
    script = COMPONENT_DIR / 'scripts' / 'split_maf.py'
    spec = f"""
    mkdir -p {split_dir}
    {run_prefix} python {script} {maf_path} --window-size {window_size} \\
        --output-dir {split_dir} {reference_flag} {species_flag}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# task template function
def concat_windows(csv_paths, split_dir, decode_dir, suffix, output_path,
                   run_prefix='', cores=1, memory='8g', walltime='4:00:00',
                   **extra_options):
    """
    Stack per-window decode CSVs into one genome-wide CSV with
    scripts/concat_windows.py, prefixing each row with window label,
    chromosome, block start coordinate and source-block index.
    """
    inputs = [os.path.join(split_dir, 'windows.tsv'),
              os.path.join(split_dir, 'blocks.tsv')] + list(csv_paths)
    outputs = {'csv': output_path}
    options = {'cores': cores, 'memory': memory, 'walltime': walltime,
               **extra_options}
    script = COMPONENT_DIR / 'scripts' / 'concat_windows.py'
    spec = f"""
    mkdir -p {os.path.dirname(output_path)}
    {run_prefix} python {script} --split-dir {split_dir} --decode-dir {decode_dir} \\
        --suffix {suffix} --output {output_path}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# task template function
def vcf_to_maf(vcf_paths, fasta_path, samples, output_path,
               block_size=1_000_000, regions=None, haplotype=1,
               mask_uncovered=True, run_prefix='',
               cores=1, memory='16g', walltime='24:00:00', **extra_options):
    """
    Build the four-species MAF from phased (g)VCFs with
    scripts/vcf_to_maf.py: one haplotype per sample, variant alleles written
    onto the reference sequence, blocks of `block_size` reference bases.
    With `fasta_path=None`, sequences are reconstructed from the VCF records
    alone (all-sites genome VCFs; uncalled positions become N). The reserved
    sample name REF stands for the reference sequence (e.g. the outgroup
    genome the samples were mapped to).
    """
    inputs = list(vcf_paths) + ([fasta_path] if fasta_path else [])
    outputs = {'maf': output_path}
    options = {'cores': cores, 'memory': memory, 'walltime': walltime,
               **extra_options}
    fasta_flag = f'--fasta {fasta_path}' if fasta_path else ''
    regions_flag = f"--regions {' '.join(regions)}" if regions else ''
    mask_flag = '' if mask_uncovered else '--no-mask-uncovered'
    script = COMPONENT_DIR / 'scripts' / 'vcf_to_maf.py'
    spec = f"""
    mkdir -p {os.path.dirname(output_path)}
    {run_prefix} python {script} {' '.join(vcf_paths)} \\
        {fasta_flag} --samples {' '.join(samples)} \\
        --output {output_path} --block-size {block_size} \\
        --haplotype {haplotype} {regions_flag} {mask_flag}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# task template function
def itrails_optimize(maf_path, config_path, output_prefix, run_prefix='',
                     cores=32, memory='64g', walltime='24:00:00',
                     **extra_options):
    """
    Fit the TRAILS coalescent HMM to a MAF alignment with itrails-optimize.

    itrails caps its internal `n_cpu` setting at the Slurm allocation, so
    `cores` is the effective parallelism regardless of the value in the
    config file.
    """
    output_dir = os.path.dirname(output_prefix)

    inputs = [maf_path, config_path]
    # exact file names produced by itrails-optimize for a given output prefix
    outputs = {
        'best_model': f'{output_prefix}.best_model.yaml',
        'history': f'{output_prefix}.optimization_history.csv',
        'starting_params': f'{output_prefix}.starting_params.yaml',
    }
    options = {'cores': cores, 'memory': memory, 'walltime': walltime,
               **extra_options}

    # no tmp-file-then-move here: itrails writes best_model.yaml at iteration
    # zero and updates it as the optimization improves, so partial files are
    # expected and failure detection is left to the backend's job status
    spec = f"""
    mkdir -p {output_dir}
    {run_prefix} itrails-optimize {config_path} --input {maf_path} --output {output_prefix}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def _itrails_decode(subcommand, best_model_path, maf_path, output_prefix,
                    run_prefix='', reference=None,
                    cores=8, memory='32g', walltime='12:00:00',
                    **extra_options):
    kind = subcommand.split('-')[-1]  # 'viterbi' or 'posterior'

    output_dir = os.path.dirname(output_prefix)

    inputs = [best_model_path, maf_path]
    outputs = {
        kind: f'{output_prefix}.{kind}.csv',
        'hidden_states': f'{output_prefix}.hidden_states.csv',
    }
    options = {'cores': cores, 'memory': memory, 'walltime': walltime,
               **extra_options}

    reference_flag = f'--reference {reference}' if reference else ''

    # if {prefix}.hidden_states.csv already exists (a previous run), itrails
    # silently writes hidden_states_2.csv instead, leaving the declared output
    # stale and the target permanently out of date -- so remove both first
    spec = f"""
    mkdir -p {output_dir}
    rm -f {outputs['hidden_states']} {output_prefix}.hidden_states_2.csv
    {run_prefix} {subcommand} --config-file {best_model_path} --input {maf_path} \\
        --output {output_prefix} --n_cpu {cores} {reference_flag}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# task template function
def itrails_viterbi(best_model_path, maf_path, output_prefix, **kwargs):
    """
    Decode the most likely hidden-state (gene-tree) sequence along the
    alignment with itrails-viterbi, using the fitted model.
    """
    return _itrails_decode('itrails-viterbi', best_model_path, maf_path,
                           output_prefix, **kwargs)


# task template function
def itrails_posterior(best_model_path, maf_path, output_prefix, **kwargs):
    """
    Compute per-position posterior probabilities of every hidden state with
    itrails-posterior, using the fitted model.
    """
    return _itrails_decode('itrails-posterior', best_model_path, maf_path,
                           output_prefix, **kwargs)


def itrails_workflow(gwf=None, analyses=None, output_dir='steps/itrails',
                     account=None, run_prefix=None,
                     window_size=None, fit='genome',
                     optimize_options=None, viterbi_options=None,
                     posterior_options=None, split_options=None,
                     concat_options=None, vcf2maf_options=None):
    """
    Add iTRAILS targets to a gwf workflow and return (gwf, targets).

    Parameters
    ----------
    gwf :
        Existing Workflow to add targets to. If None, a new one is created
        with the current working directory as its working dir, so the module
        can also back a standalone `gwf run` in this repository.
    analyses :
        List of dicts, one per analysis. Input is either an alignment or
        phased variant calls (exactly one of `maf` / `vcf`):
          name        -- label used in target names and as output subdirectory
          maf         -- path to the MAF alignment (A, B, C, outgroup)
          vcf         -- path (or list of paths) to phased (g)VCF file(s),
                         all called against the same reference genome, so
                         that positions share that coordinate system (e.g.
                         baboon samples mapped to the macaque genome). An
                         alignment step reconstructs one haplotype per
                         sample and writes alignment/{name}.maf
          fasta       -- with `vcf`: the reference FASTA the VCFs were
                         called against. Optional: without it the sequences
                         are reconstructed from the VCF records alone, which
                         suits all-sites genome VCFs (a record with its REF
                         base at every callable position); positions without
                         a record become N, as do the interiors of banded
                         gVCF non-variant blocks (whose bases are not in the
                         file — give `fasta` for such files)
          samples     -- with `vcf`: exactly 4 sample IDs in species order
                         A, B, C, outgroup; the first haplotype of each is
                         used (override with `haplotype`). The reserved name
                         REF means the reference sequence itself (e.g. the
                         macaque genome as outgroup), taken from `fasta` or,
                         without it, assembled from the records' REF fields.
                         Sample IDs double as species names: use them in the
                         config's species_list, and avoid '.' in IDs
          block_size  -- with `vcf`: MAF block length in reference bases
                         (default 1_000_000); blocks delimit independent
                         HMM sequences and the units of windowing
          regions     -- with `vcf`: optional list of contigs to analyze
                         (default: every contig in the FASTA)
          haplotype   -- with `vcf`: 1-based haplotype to extract (default 1)
          mask_uncovered -- with `vcf`: positions without any VCF record
                         become N, the gVCF reading (default True); set
                         False for plain VCFs where no record means
                         reference-equal (needs `fasta`)
          config      -- path to the itrails-optimize YAML config
          reference   -- optional species name; polarizes decode coordinates
                         and, in windowed mode, defines the window bins
                         (default: the first s-line of each MAF block; with
                         `vcf`, must be one of `samples`)
          window_size -- optional; see below. Overrides the factory argument
          fit         -- 'genome' or 'window'; overrides the factory argument
          species     -- optional list; in windowed mode, drop MAF blocks
                         missing any of these species (set it to the iTRAILS
                         species_list so the split and the iTRAILS parser
                         agree on block indices; defaults to `samples` for
                         VCF input)
          optimize_options / viterbi_options / posterior_options /
          split_options / concat_options / vcf2maf_options
                      -- optional per-analysis resource overrides
        Paths are resolved relative to the directory `gwf run` is invoked
        from (the parent project root when used as a submodule).
    output_dir :
        Each analysis writes its files to {output_dir}/{name}/, with one
        subfolder per step: split/ (window MAFs and manifests), optimize/
        (model fits), viterbi/ and posterior/ (decode CSVs and state
        legends), concat/ (genome-wide concatenated CSVs).
    account :
        Slurm account added to every target's options.
    run_prefix :
        Command prefix putting the itrails executables on PATH. Defaults to
        running inside this repository's pixi `itrails` environment; a parent
        using its own environment can pass e.g.
        'eval "$(conda shell.bash activate itrails)" &&'-style prefix or ''.
    window_size :
        If set (e.g. 2_000_000), split each alignment into windows of that
        many reference bases at MAF-block boundaries and decode per window;
        window MAFs go to split/, per-window decode files to viterbi/ and
        posterior/, and the decode CSVs are concatenated into
        concat/{name}.viterbi.csv and concat/{name}.posterior.csv
        (annotated with window/chrom/block coordinates).
        Window boundaries are derived from the alignment at workflow-
        definition time (cached by file size/mtime), so **the MAF must exist
        when the workflow is defined** — if an upstream target produces it,
        run that part of the workflow first. VCF analyses are exempt: their
        block structure is deterministic, so windows are derived from contig
        lengths instead — the FASTA index ({fasta}.fai, created with
        samtools faidx), or without a fasta the ##contig header lines of
        the VCFs (which must then exist and carry length=).
    fit :
        'genome' (default): one itrails-optimize on the whole alignment, all
        windows decoded with that model. 'window': a separate optimize per
        window (each window decoded with its own model — for studying
        parameter variation along the genome; expect noisy estimates from
        small windows and one heavy job per window).
    optimize_options, viterbi_options, posterior_options, split_options,
    concat_options, vcf2maf_options :
        Dicts of gwf options (cores, memory, walltime, ...) applied to the
        corresponding step of every analysis. In windowed mode the decode
        (and, with fit='window', optimize) options apply to every per-window
        target, so consider sizing them for a single window.

    Returns
    -------
    (gwf, targets) where targets maps step names to lists of created targets:
    'optimize', 'viterbi', 'posterior', for VCF input also 'vcf2maf', and in
    windowed mode also 'split', 'concat_viterbi', 'concat_posterior'.
    Downstream targets in a parent workflow should depend on e.g.
    targets['viterbi'][0].outputs['viterbi'], or in windowed mode
    targets['concat_viterbi'][0].outputs['csv'].
    """
    if gwf is None:
        gwf = Workflow(working_dir=os.getcwd())
    if analyses is None:
        analyses = []
    if run_prefix is None:
        run_prefix = default_run_prefix()

    # definition-time file access (the window scan) must resolve relative
    # paths the same way gwf resolves target files: against the working dir
    working_dir = getattr(gwf, 'working_dir', None) or os.getcwd()

    common_options = {}
    if account is not None:
        common_options['account'] = account

    def merged(step_options, analysis, key):
        return {**common_options, **(step_options or {}),
                **(analysis.get(key) or {})}

    targets = defaultdict(list)
    for analysis in analyses:
        name = analysis['name']
        # target names must be valid identifiers
        safe_name = re.sub(r'\W+', '_', name)
        # one subfolder per step under {output_dir}/{name}/
        analysis_dir = os.path.join(output_dir, name)
        step_dir = {step: os.path.join(analysis_dir, step)
                    for step in ('alignment', 'split', 'optimize', 'viterbi',
                                 'posterior', 'concat')}
        reference = analysis.get('reference')

        if ('maf' in analysis) == ('vcf' in analysis):
            raise ValueError(f"analysis '{name}': give exactly one of 'maf' "
                             f"(alignment) or 'vcf' (phased (g)VCF input)")

        vcf2maf_target = None
        samples = None
        analysis_block_size = analysis.get('block_size', 1_000_000)
        if 'vcf' in analysis:
            vcfs = analysis['vcf']
            if isinstance(vcfs, str):
                vcfs = [vcfs]
            if not analysis.get('samples'):
                raise ValueError(f"analysis '{name}': VCF input also needs "
                                 f"'samples' (4 sample IDs). 'fasta' is "
                                 f"optional for all-sites genome VCFs")
            if not analysis.get('fasta') and analysis.get('mask_uncovered') is False:
                raise ValueError(f"analysis '{name}': mask_uncovered: false "
                                 f"needs 'fasta' (without the reference there "
                                 f"is nothing to fill uncovered positions "
                                 f"with)")
            samples = list(analysis['samples'])
            if len(samples) != 4:
                raise ValueError(f"analysis '{name}': 'samples' must list "
                                 f"exactly 4 sample IDs in species order "
                                 f"A, B, C, outgroup; got {len(samples)}")
            if reference is not None and reference not in samples:
                raise ValueError(f"analysis '{name}': 'reference' must be one "
                                 f"of the 4 samples with VCF input, got "
                                 f"{reference!r}")
            vcf2maf_target = gwf.target_from_template(
                f'itrails_vcf2maf_{safe_name}',
                vcf_to_maf(vcfs, analysis.get('fasta'), samples,
                           os.path.join(step_dir['alignment'], f'{name}.maf'),
                           block_size=analysis_block_size,
                           regions=analysis.get('regions'),
                           haplotype=analysis.get('haplotype', 1),
                           mask_uncovered=analysis.get('mask_uncovered', True),
                           run_prefix=run_prefix,
                           **merged(vcf2maf_options, analysis,
                                    'vcf2maf_options')))
            targets['vcf2maf'].append(vcf2maf_target)
            maf = vcf2maf_target.outputs['maf']
        else:
            maf = analysis['maf']

        opt_options = merged(optimize_options, analysis, 'optimize_options')
        vit_options = merged(viterbi_options, analysis, 'viterbi_options')
        post_options = merged(posterior_options, analysis, 'posterior_options')

        analysis_window_size = analysis.get('window_size', window_size)

        if not analysis_window_size:
            # unwindowed: one optimize + one decode pair on the whole
            # alignment; the per-step folders keep the two decode steps'
            # hidden_states.csv files apart
            optimize_target = gwf.target_from_template(
                f'itrails_optimize_{safe_name}',
                itrails_optimize(maf, analysis['config'],
                                 os.path.join(step_dir['optimize'], name),
                                 run_prefix=run_prefix, **opt_options))
            targets['optimize'].append(optimize_target)

            best_model = optimize_target.outputs['best_model']

            targets['viterbi'].append(gwf.target_from_template(
                f'itrails_viterbi_{safe_name}',
                itrails_viterbi(best_model, maf,
                                os.path.join(step_dir['viterbi'], name),
                                run_prefix=run_prefix, reference=reference,
                                **vit_options)))
            targets['posterior'].append(gwf.target_from_template(
                f'itrails_posterior_{safe_name}',
                itrails_posterior(best_model, maf,
                                  os.path.join(step_dir['posterior'], name),
                                  run_prefix=run_prefix, reference=reference,
                                  **post_options)))
            continue

        # windowed mode
        analysis_fit = analysis.get('fit', fit)
        if analysis_fit not in ('genome', 'window'):
            raise ValueError(f"analysis '{name}': fit must be 'genome' or "
                             f"'window', got {analysis_fit!r}")
        species = analysis.get('species') or samples

        if vcf2maf_target is not None:
            # the MAF does not exist yet, but vcf_to_maf writes blocks of
            # block_size tiling every contig, so the windows follow from the
            # contig lengths alone: from the FASTA index when a fasta is
            # given, else from the ##contig lines of the VCF headers
            if analysis.get('fasta'):
                fai = str(analysis['fasta']) + '.fai'
                fai_on_disk = (fai if os.path.isabs(fai)
                               else os.path.join(working_dir, fai))
                if not os.path.exists(fai_on_disk):
                    raise FileNotFoundError(
                        f"analysis '{name}': windowed VCF analyses derive "
                        f"their window boundaries from the FASTA index; "
                        f"create {fai} with 'samtools faidx "
                        f"{analysis['fasta']}'")
                contigs = _read_fai(fai_on_disk)
                length_source = fai
            else:
                vcfs_on_disk = [v if os.path.isabs(v)
                                else os.path.join(working_dir, v)
                                for v in vcfs]
                absent = [v for v in vcfs_on_disk if not os.path.exists(v)]
                if absent:
                    raise FileNotFoundError(
                        f"analysis '{name}': windowed VCF analyses without a "
                        f"fasta derive their window boundaries from the VCF "
                        f"headers, so the files must exist when the workflow "
                        f"is defined; missing: {absent}")
                try:
                    contigs = _script_module('vcf_to_maf').vcf_contigs(
                        vcfs_on_disk)
                except ValueError as error:
                    raise ValueError(f"analysis '{name}': {error}") from None
                length_source = 'the VCF headers'
            regions = analysis.get('regions')
            if regions:
                known = {contig for contig, _ in contigs}
                unknown = [region for region in regions if region not in known]
                if unknown:
                    raise ValueError(f"analysis '{name}': regions {unknown} "
                                     f"not in {length_source}")
                keep = set(regions)
                contigs = [c for c in contigs if c[0] in keep]
            no_length = [contig for contig, length in contigs
                         if length is None]
            if no_length:
                raise ValueError(
                    f"analysis '{name}': contig(s) {no_length} have no "
                    f"length in the VCF headers; windowed VCF analyses "
                    f"without a fasta need length= on every ##contig line")
            windows = _windows_from_lengths(contigs, analysis_block_size,
                                            analysis_window_size,
                                            reference or samples[0])
        else:
            maf_on_disk = (maf if os.path.isabs(maf)
                           else os.path.join(working_dir, maf))
            if not os.path.exists(maf_on_disk):
                raise FileNotFoundError(
                    f"analysis '{name}': windowed analyses derive their window "
                    f"boundaries from the alignment at workflow-definition "
                    f"time, so the MAF must already exist; {maf} not found. "
                    f"If an upstream target produces it, run that first.")
            cache_dir = (step_dir['split'] if os.path.isabs(step_dir['split'])
                         else os.path.join(working_dir, step_dir['split']))
            windows = _window_index(maf_on_disk, analysis_window_size,
                                    reference, species, cache_dir)
        if not windows:
            raise ValueError(f"analysis '{name}': no alignment blocks with "
                             f"the binning species found in {maf}")

        split_target = gwf.target_from_template(
            f'itrails_split_{safe_name}',
            split_maf(maf, step_dir['split'], analysis_window_size,
                      [os.path.join(step_dir['split'], f"{w['label']}.maf")
                       for w in windows],
                      reference=reference, species=species,
                      run_prefix=run_prefix,
                      **merged(split_options, analysis, 'split_options')))
        targets['split'].append(split_target)
        window_maf = dict(zip((w['label'] for w in windows),
                              split_target.outputs['windows']))

        if analysis_fit == 'genome':
            optimize_target = gwf.target_from_template(
                f'itrails_optimize_{safe_name}',
                itrails_optimize(maf, analysis['config'],
                                 os.path.join(step_dir['optimize'], name),
                                 run_prefix=run_prefix, **opt_options))
            targets['optimize'].append(optimize_target)
            best_model = {label: optimize_target.outputs['best_model']
                          for label in window_maf}
        else:
            best_model = {}
            for label in window_maf:
                optimize_target = gwf.target_from_template(
                    f'itrails_optimize_{safe_name}_{label}',
                    itrails_optimize(window_maf[label], analysis['config'],
                                     os.path.join(step_dir['optimize'], label),
                                     run_prefix=run_prefix, **opt_options))
                targets['optimize'].append(optimize_target)
                best_model[label] = optimize_target.outputs['best_model']

        decode_csvs = {'viterbi': [], 'posterior': []}
        for label in window_maf:
            for kind, template, opts in (
                    ('viterbi', itrails_viterbi, vit_options),
                    ('posterior', itrails_posterior, post_options)):
                decode_target = gwf.target_from_template(
                    f'itrails_{kind}_{safe_name}_{label}',
                    template(best_model[label], window_maf[label],
                             os.path.join(step_dir[kind], label),
                             run_prefix=run_prefix, reference=reference,
                             **opts))
                targets[kind].append(decode_target)
                decode_csvs[kind].append(decode_target.outputs[kind])

        cat_options = merged(concat_options, analysis, 'concat_options')
        for kind in ('viterbi', 'posterior'):
            targets[f'concat_{kind}'].append(gwf.target_from_template(
                f'itrails_concat_{kind}_{safe_name}',
                concat_windows(decode_csvs[kind], step_dir['split'],
                               step_dir[kind], f'.{kind}.csv',
                               os.path.join(step_dir['concat'],
                                            f'{name}.{kind}.csv'),
                               run_prefix=run_prefix, **cat_options)))

    return gwf, targets
