"""
Concatenate per-window iTRAILS decode CSVs into one genome-wide CSV.

Reads windows.tsv and blocks.tsv written by split_maf.py from --split-dir,
then stacks {decode-dir}/{window}{suffix} for every window in manifest
order. iTRAILS reports a window-local block index (column 0 of both the
viterbi and posterior CSVs) and block-local positions, so each row is
prefixed with columns that restore genome context:

    window       window label (w0000, ...)
    chrom        binning-reference chromosome of the block (e.g. hg38.chr1)
    block_start  reference start coordinate of the block
    src_block    0-based index of the block in the original, unsplit MAF

followed by the original columns unchanged. The mapping assumes iTRAILS
kept every block in the window MAF (guaranteed when split_maf.py was run
with --species matching the iTRAILS species_list); a block index outside
the manifest is an error.
"""

import argparse
import csv
import sys
from pathlib import Path


def read_manifests(split_dir):
    split_dir = Path(split_dir)
    with open(split_dir / 'windows.tsv') as f:
        windows = [row['window'] for row in csv.DictReader(f, delimiter='\t')]
    blocks = {}
    with open(split_dir / 'blocks.tsv') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            blocks[(row['window'], int(row['local_idx']))] = row
    return windows, blocks


def concat(split_dir, decode_dir, suffix, output_path):
    windows, blocks = read_manifests(split_dir)
    header = None
    with open(output_path, 'w', newline='') as out:
        writer = csv.writer(out)
        for window in windows:
            path = Path(decode_dir) / f'{window}{suffix}'
            with open(path, newline='') as f:
                reader = csv.reader(f)
                file_header = next(reader)
                if header is None:
                    header = file_header
                    writer.writerow(['window', 'chrom', 'block_start', 'src_block']
                                    + header)
                elif file_header != header:
                    sys.exit(f'error: column mismatch in {path}: '
                             f'{file_header} vs {header}')
                for row in reader:
                    local_idx = int(float(row[0]))
                    block = blocks.get((window, local_idx))
                    if block is None:
                        sys.exit(
                            f'error: {path} refers to block {local_idx}, not in '
                            f'blocks.tsv for {window}. iTRAILS probably skipped '
                            f'blocks missing a species; re-split with --species '
                            f'(the `species` analysis key) so both drop the '
                            f'same blocks.')
                    writer.writerow([window, block['chrom'], block['start'],
                                     block['src_block']] + row)
    print(f'wrote {output_path} ({len(windows)} windows)')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--split-dir', required=True,
                        help='directory with the windows.tsv and blocks.tsv '
                             'manifests written by split_maf.py')
    parser.add_argument('--decode-dir', required=True,
                        help='directory with the per-window decode CSVs')
    parser.add_argument('--suffix', required=True,
                        help="per-window file suffix after the window label, "
                             "e.g. '.viterbi.csv'")
    parser.add_argument('--output', required=True, help='concatenated CSV path')
    args = parser.parse_args()
    concat(args.split_dir, args.decode_dir, args.suffix, args.output)


if __name__ == '__main__':
    main()
