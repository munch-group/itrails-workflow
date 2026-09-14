"""
Split a MAF alignment into fixed-size genomic windows at block boundaries.

Alignment blocks are binned by the start coordinate of one species' s-line
(the *binning reference*: --reference, or the first s-line of each block —
the MAF convention for the reference genome). Blocks are copied verbatim and
never split, so windows are approximate at block granularity: a block that
straddles a boundary belongs to the window its start falls in. iTRAILS
treats each MAF block as an independent sequence in its HMM, so splitting at
block boundaries loses no information.

Outputs, in --output-dir:

    w0000.maf, w0001.maf, ...  one MAF per non-empty window
    windows.tsv                window, chrom, start, end, n_blocks
    blocks.tsv                 window, local_idx, src_block, chrom, start, length

`src_block` is the 0-based index of the block in the input file, `local_idx`
its 0-based position within the window's MAF — the Block_idx that iTRAILS
reports when decoding that window (provided iTRAILS keeps every block; pass
--species to drop blocks missing a species, mirroring iTRAILS' own parser,
so the two stay in step).

The functions `scan_blocks` and `assign_windows` are also imported by
itrails_workflow.py at workflow-definition time, so the gwf target list and
the files this script writes are derived from the same binning code.
"""

import argparse
import sys
from pathlib import Path


def window_label(idx):
    return f'w{idx:04d}'


def _sline_fields(line):
    # s <src> <start> <size> <strand> <srcSize> <text>
    fields = line.split()
    return {'src': fields[1], 'species': fields[1].split('.')[0],
            'start': int(fields[2]), 'size': int(fields[3])}


def scan_blocks(maf_path, reference=None, species=None):
    """
    Return one dict per alignment block in the file:

        {'src_block': int, 'chrom': str, 'start': int, 'length': int,
         'keep': bool}

    chrom/start/length come from the s-line of `reference` (or the first
    s-line if None). keep is False for blocks lacking the binning species
    or, if `species` is given, any of those species.
    """
    required = set(species) if species else None
    blocks = []
    slines = None

    def close_block():
        if slines is None:
            return
        by_species = {s['species']: s for s in slines}
        ref = by_species.get(reference) if reference else (slines[0] if slines else None)
        keep = ref is not None and (required is None or required <= by_species.keys())
        blocks.append({
            'src_block': len(blocks),
            'chrom': ref['src'] if ref else None,
            'start': ref['start'] if ref else None,
            'length': ref['size'] if ref else None,
            'keep': keep,
        })

    with open(maf_path) as f:
        for line in f:
            if line.startswith('a'):
                close_block()
                slines = []
            elif line.startswith('s') and slines is not None:
                slines.append(_sline_fields(line))
        close_block()
    return blocks


def assign_windows(blocks, window_size):
    """
    Bin kept blocks into windows of `window_size` reference bases and return
    an ordered list of dicts:

        {'label': 'w0000', 'chrom': str, 'start': int, 'end': int,
         'blocks': [block dicts in file order]}

    start/end are the nominal bin bounds (a block starting inside the bin may
    extend past `end`). Only non-empty windows are returned, ordered by
    chromosome first appearance in the file, then bin start.
    """
    bins = {}
    chrom_order = []
    for block in blocks:
        if not block['keep']:
            continue
        if block['chrom'] not in chrom_order:
            chrom_order.append(block['chrom'])
        bins.setdefault((block['chrom'], block['start'] // window_size), []).append(block)
    windows = []
    for chrom in chrom_order:
        for chrom_, bin_idx in sorted(k for k in bins if k[0] == chrom):
            windows.append({
                'label': window_label(len(windows)),
                'chrom': chrom,
                'start': bin_idx * window_size,
                'end': (bin_idx + 1) * window_size,
                'blocks': bins[(chrom, bin_idx)],
            })
    return windows


def split(maf_path, output_dir, window_size, reference=None, species=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    blocks = scan_blocks(maf_path, reference=reference, species=species)
    windows = assign_windows(blocks, window_size)
    if not windows:
        sys.exit(f'error: no alignment blocks with the binning species found in {maf_path}')
    label_of = {b['src_block']: w['label'] for w in windows for b in w['blocks']}

    # remove window files from a previous split so a re-split with fewer
    # windows leaves no stale w*.maf behind
    for old in output_dir.glob('w[0-9]*.maf'):
        old.unlink()

    header = '##maf version=1 scoring=none\n'
    for window in windows:
        (output_dir / f"{window['label']}.maf").write_text(header)

    # second pass: append each kept block verbatim to its window's file
    with open(maf_path) as f:
        src_block, out, buffer = -1, None, []

        def flush():
            if out is not None and buffer:
                with open(output_dir / f'{out}.maf', 'a') as w:
                    w.write('\n')
                    w.writelines(buffer)

        for line in f:
            if line.startswith('a'):
                flush()
                src_block += 1
                out, buffer = label_of.get(src_block), []
            if out is not None and line.strip():
                buffer.append(line)
        flush()

    with open(output_dir / 'windows.tsv', 'w') as f:
        f.write('window\tchrom\tstart\tend\tn_blocks\n')
        for w in windows:
            f.write(f"{w['label']}\t{w['chrom']}\t{w['start']}\t{w['end']}\t{len(w['blocks'])}\n")

    with open(output_dir / 'blocks.tsv', 'w') as f:
        f.write('window\tlocal_idx\tsrc_block\tchrom\tstart\tlength\n')
        for w in windows:
            for local_idx, b in enumerate(w['blocks']):
                f.write(f"{w['label']}\t{local_idx}\t{b['src_block']}\t"
                        f"{b['chrom']}\t{b['start']}\t{b['length']}\n")

    dropped = sum(not b['keep'] for b in blocks)
    print(f"{len(blocks)} blocks ({dropped} dropped) -> "
          f"{len(windows)} windows in {output_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('maf', help='input MAF alignment')
    parser.add_argument('--window-size', type=int, required=True,
                        help='window size in reference bases')
    parser.add_argument('--output-dir', required=True,
                        help='directory for window MAFs and manifests')
    parser.add_argument('--reference', default=None,
                        help='species whose coordinates define the windows '
                             '(default: first s-line of each block)')
    parser.add_argument('--species', nargs='+', default=None,
                        help='drop blocks missing any of these species '
                             '(mirror of the iTRAILS species_list)')
    args = parser.parse_args()
    split(args.maf, args.output_dir, args.window_size,
          reference=args.reference, species=args.species)


if __name__ == '__main__':
    main()
