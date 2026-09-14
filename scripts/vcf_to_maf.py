"""
Build a MAF alignment for iTRAILS from phased (g)VCFs, with or without the
reference FASTA.

Use case: closely related species that were all mapped and called against the
same reference genome (e.g. baboon species against macaque), so that variant
coordinates place every sample in the reference coordinate system and the
genomes are implicitly aligned. For each of the requested samples, the first
haplotype (or --haplotype H) is reconstructed; the reserved sample name REF
denotes the reference sequence itself (e.g. the macaque genome as the
outgroup).

Two modes:

  With --fasta, the reference sequence provides the bases everywhere and
  variant alleles are written on top of it.

  Without --fasta, sequences are reconstructed from the VCF records alone.
  This is meant for all-sites "genome" VCFs, where every callable position
  has a record carrying its REF base: any position without a record for a
  sample is written as N, and the REF pseudo-sample is assembled from the
  REF fields of the records. Contig names, lengths and order come from the
  ##contig header lines (lengths are required). Note that banded gVCF
  non-variant blocks (INFO END=...) only carry the base at their anchor
  position, so in this mode the rest of each block is N — supply --fasta
  for such files.

Reconstruction rules per sample and position:

  - SNP / same-length substitution: the chosen allele replaces the REF bases.
  - length-changing or symbolic allele (<NON_REF>, <*>, *): the alignment
    cannot represent it in reference coordinates, so the REF span becomes N.
  - missing allele (.): the REF span becomes N.
  - gVCF non-variant blocks (INFO END=...) with allele 0 mark the span as
    reference-covered (with --fasta) or contribute their anchor base only
    (without --fasta).
  - positions not covered by any record are N. With --fasta, pass
    --no-mask-uncovered for plain VCFs where absence of a record means
    reference-equal (impossible to honor without the reference).
  - FILTER is ignored; pre-filter the VCF if needed.

Samples may be spread over several files (per-sample or per-region VCFs);
each file contributes genotypes for whichever requested samples its header
contains. Files must be coordinate-sorted with contigs in the same order as
the FASTA / the first file's header (the standard reference-dict order). If
two files carry the same sample, later files win where they overlap.

The alignment is written in blocks of --block-size reference bases:

    s {sample}.{contig} {start} {size} + {contig_length} {sequence}

so downstream windowing splits at these block boundaries. Sample names are
used as MAF species names and must not contain '.' (iTRAILS takes the
species as the part before the first dot).
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np

GT_SPLIT = re.compile(r'[|/]')
CONTIG_RE = re.compile(r'##contig=<(.*)>')


def open_text(path):
    path = str(path)
    if path.endswith('.gz') or path.endswith('.bgz'):
        return gzip.open(path, 'rt')
    return open(path)


def read_fasta(path):
    """Yield (name, sequence) one contig at a time (memory-light)."""
    name, chunks = None, []
    with open_text(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, ''.join(chunks)


def fasta_contig_names(path):
    """Contig names in file order, from {path}.fai if present, else a scan."""
    fai = Path(str(path) + '.fai')
    if fai.exists():
        with open(fai) as f:
            return [line.split('\t')[0] for line in f if line.strip()]
    return [name for name, _ in read_fasta(path)]


def vcf_contigs(vcf_paths):
    """
    [(contig, length_or_None)] from the ##contig header lines, in processing
    order: the first file's header order, then contigs later files add.
    Raises ValueError on inconsistent lengths.
    """
    order, length_of = [], {}
    for path in vcf_paths:
        with open_text(path) as f:
            for line in f:
                if not line.startswith('##'):
                    break
                match = CONTIG_RE.match(line)
                if not match:
                    continue
                fields = dict(part.split('=', 1)
                              for part in match.group(1).split(',')
                              if '=' in part)
                name = fields.get('ID')
                if name is None:
                    continue
                length = int(fields['length']) if 'length' in fields else None
                if name not in length_of:
                    order.append(name)
                    length_of[name] = length
                elif length is not None:
                    if length_of[name] is None:
                        length_of[name] = length
                    elif length_of[name] != length:
                        raise ValueError(
                            f'contig {name} has inconsistent lengths across '
                            f'VCF headers ({length_of[name]} vs {length} '
                            f'in {path})')
    return [(name, length_of[name]) for name in order]


class VCFReader:
    """Streams records of one VCF, exposing the requested samples it holds."""

    def __init__(self, path, samples):
        self.path = path
        self.file = open_text(path)
        self._buffer = None
        self.done_contigs = set()
        for line in self.file:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                columns = line.rstrip('\n').split('\t')
                header_samples = columns[9:]
                self.sample_column = {s: 9 + header_samples.index(s)
                                      for s in samples if s in header_samples}
                return
            break
        sys.exit(f'error: {path} has no #CHROM header line')

    def _next(self):
        if self._buffer is not None:
            record, self._buffer = self._buffer, None
            return record
        for line in self.file:
            if line.strip():
                return line.rstrip('\n').split('\t')
        return None

    def records_for(self, contig, future_contigs):
        """Yield this file's records for `contig`, assuming coordinate-sorted
        input; records for non-target contigs are dropped, records for a
        later target contig end the iteration (kept for later)."""
        while True:
            record = self._next()
            if record is None:
                return
            chrom = record[0]
            if chrom == contig:
                yield record
            elif chrom in future_contigs:
                self._buffer = record
                return
            elif chrom in self.done_contigs:
                sys.exit(f'error: {self.path} is not sorted in the same '
                         f'contig order as the FASTA/first VCF header '
                         f'(saw {chrom} again)')
            # else: contig not requested -- drop the record


def apply_records(reader, contig, future_contigs, haplotypes, covered,
                  ref_seq, haplotype, stats):
    """
    Write one file's variants for `contig` into the haplotype arrays.
    `covered` is a dict of coverage masks (FASTA mode) or None (VCF-only
    mode, where the background is already N); `ref_seq` is an array to
    accumulate REF bases into (VCF-only mode) or None.
    """
    for record in reader.records_for(contig, future_contigs):
        pos0 = int(record[1]) - 1
        ref = record[3]
        alts = record[4].split(',') if record[4] != '.' else []
        info = record[7]
        fmt = record[8].split(':')
        try:
            gt_idx = fmt.index('GT')
        except ValueError:
            continue
        end_match = re.search(r'(?:^|;)END=(\d+)', info)
        end0 = int(end_match.group(1)) if end_match else None
        if ref_seq is not None:
            ref_seq[pos0:pos0 + len(ref)] = np.frombuffer(
                ref.upper().encode(), dtype='S1')

        for sample, column in reader.sample_column.items():
            gt = record[column].split(':')[gt_idx]
            alleles = GT_SPLIT.split(gt)
            if len(alleles) > 1 and '/' in gt and len(set(alleles)) > 1:
                stats['unphased_het'] += 1
            allele = alleles[haplotype - 1] if haplotype <= len(alleles) else '.'
            hap = haplotypes[sample]

            if end0 is not None and allele == '0':
                # gVCF non-variant block: sample matches the reference
                if covered is not None:
                    covered[sample][pos0:end0] = True
                else:
                    # only the anchor base is in the file
                    hap[pos0:pos0 + len(ref)] = np.frombuffer(
                        ref.upper().encode(), dtype='S1')
                    stats['no_base'] += max(0, end0 - pos0 - len(ref))
                continue
            span = end0 - pos0 if end0 is not None else len(ref)
            if allele == '.':
                hap[pos0:pos0 + span] = b'N'
                if covered is not None:
                    covered[sample][pos0:pos0 + span] = True
                stats['missing'] += 1
                continue
            chosen = ref if allele == '0' else alts[int(allele) - 1]
            if chosen.startswith('<') or chosen == '*' or len(chosen) != len(ref) \
                    or end0 is not None:
                # not representable in reference coordinates
                hap[pos0:pos0 + span] = b'N'
                stats['masked'] += 1
            else:
                hap[pos0:pos0 + span] = np.frombuffer(
                    chosen.upper().encode(), dtype='S1')
                if chosen != ref:
                    stats['applied'] += 1
            if covered is not None:
                covered[sample][pos0:pos0 + span] = True


def vcf_to_maf(vcf_paths, fasta_path, samples, output_path,
               block_size=1_000_000, regions=None, haplotype=1,
               mask_uncovered=True):
    for sample in samples:
        if '.' in sample:
            sys.exit(f"error: sample name {sample!r} contains '.', which "
                     f"iTRAILS uses to separate species from chromosome")
    if fasta_path is None and not mask_uncovered:
        sys.exit('error: --no-mask-uncovered needs --fasta (without the '
                 'reference there is nothing to fill uncovered positions '
                 'with)')
    vcf_samples = [s for s in samples if s != 'REF']
    readers = [VCFReader(path, vcf_samples) for path in vcf_paths]
    found = set().union(*(r.sample_column for r in readers)) if readers else set()
    missing = [s for s in vcf_samples if s not in found]
    if missing:
        sys.exit(f'error: sample(s) {missing} not found in any VCF')

    if fasta_path is not None:
        contigs = [(name, None) for name in fasta_contig_names(fasta_path)]
    else:
        try:
            contigs = vcf_contigs(vcf_paths)
        except ValueError as error:
            sys.exit(f'error: {error}')
    if regions:
        known = {name for name, _ in contigs}
        unknown = [region for region in regions if region not in known]
        if unknown:
            sys.exit(f'error: region(s) {unknown} not found in the '
                     f'{"FASTA" if fasta_path else "VCF headers"}')
        keep = set(regions)
        contigs = [c for c in contigs if c[0] in keep]
    if fasta_path is None:
        no_length = [name for name, length in contigs if length is None]
        if no_length:
            sys.exit(f'error: contig(s) {no_length} have no length in the '
                     f'VCF headers; without --fasta, ##contig lines must '
                     f'carry length=')
    contig_names = [name for name, _ in contigs]
    target_set = set(contig_names)

    if fasta_path is not None:
        sequence_source = ((name, seq) for name, seq in read_fasta(fasta_path)
                           if name in target_set)
    else:
        sequence_source = ((name, length) for name, length in contigs)

    stats = {'applied': 0, 'masked': 0, 'missing': 0, 'unphased_het': 0,
             'no_base': 0}
    n_blocks = 0
    with open(output_path, 'w') as out:
        out.write('##maf version=1 scoring=none\n')
        out.write(f'# built by scripts/vcf_to_maf.py from '
                  f'{len(vcf_paths)} VCF(s), haplotype {haplotype}'
                  + (', no reference FASTA' if fasta_path is None else '')
                  + '\n')
        for name, payload in sequence_source:
            future = set(contig_names[contig_names.index(name) + 1:])
            haplotypes, covered, ref_seq = {}, None, None
            if fasta_path is not None:
                reference = np.frombuffer(payload.upper().encode(), dtype='S1')
                length = len(reference)
                covered = {}
                for sample in vcf_samples:
                    haplotypes[sample] = reference.copy()
                    covered[sample] = np.zeros(length, dtype=bool)
            else:
                length = payload
                ref_seq = np.full(length, b'N', dtype='S1')
                for sample in vcf_samples:
                    haplotypes[sample] = np.full(length, b'N', dtype='S1')

            for reader in readers:
                apply_records(reader, name, future, haplotypes, covered,
                              ref_seq, haplotype, stats)
                reader.done_contigs.add(name)
            if covered is not None and mask_uncovered:
                for sample in vcf_samples:
                    haplotypes[sample][~covered[sample]] = b'N'
            if 'REF' in samples:
                haplotypes['REF'] = reference if fasta_path is not None else ref_seq

            for start in range(0, length, block_size):
                end = min(start + block_size, length)
                out.write('\na score=0\n')
                for sample in samples:
                    seq = haplotypes[sample][start:end].tobytes().decode()
                    out.write(f's {sample}.{name} {start} {end - start} + '
                              f'{length} {seq}\n')
                n_blocks += 1
        out.write('\n')

    print(f'wrote {output_path}: {len(contig_names)} contig(s), {n_blocks} '
          f"blocks; {stats['applied']} alleles applied, {stats['masked']} "
          f"spans N-masked (indel/symbolic), {stats['missing']} missing "
          f"genotypes"
          + (f", {stats['unphased_het']} unphased heterozygous genotypes "
             f"(first listed allele used)" if stats['unphased_het'] else ''))
    if stats['no_base']:
        print(f"warning: {stats['no_base']} positions sit inside gVCF "
              f"non-variant blocks whose bases are not in the file and were "
              f"written as N -- pass --fasta to fill them from the reference",
              file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('vcfs', nargs='+', help='phased (g)VCF file(s), '
                        'plain or gzipped, coordinate-sorted')
    parser.add_argument('--fasta', default=None,
                        help='reference FASTA the VCFs were called against; '
                             'omit for all-sites genome VCFs to reconstruct '
                             'purely from the records (uncalled positions '
                             'become N)')
    parser.add_argument('--samples', nargs='+', required=True,
                        help='sample IDs to extract, in species order A B C '
                             'outgroup; the reserved name REF means the '
                             'reference sequence (from --fasta, or assembled '
                             'from the records\' REF fields without it)')
    parser.add_argument('--output', required=True, help='output MAF path')
    parser.add_argument('--block-size', type=int, default=1_000_000,
                        help='MAF block length in reference bases '
                             '(default 1000000)')
    parser.add_argument('--regions', nargs='+', default=None,
                        help='contigs to convert (default: all in the FASTA '
                             'or VCF headers)')
    parser.add_argument('--haplotype', type=int, default=1,
                        help='1-based haplotype to extract from each '
                             'genotype (default 1 = first)')
    parser.add_argument('--no-mask-uncovered', action='store_true',
                        help='with --fasta only: treat positions without any '
                             'VCF record as reference-equal (plain VCF) '
                             'instead of N (gVCF)')
    args = parser.parse_args()
    vcf_to_maf(args.vcfs, args.fasta, args.samples, args.output,
               block_size=args.block_size, regions=args.regions,
               haplotype=args.haplotype,
               mask_uncovered=not args.no_mask_uncovered)


if __name__ == '__main__':
    main()
