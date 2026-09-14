"""
Simulate a dummy four-species MAF alignment for the iTRAILS workflow.

One haploid genome is sampled per species on the tree
(((hg38, panTro5), gorGor5), ponAbe2) with split times and ancestral
population sizes mirroring human/chimp/gorilla/orangutan (generation
time 25 years):

    human-chimp split      240,000 generations (~6 Mya),   N_AB   =  65,000
    +gorilla               350,000 generations (~8.75 Mya), N_ABC = 100,000
    +orangutan             700,000 generations (~17.5 Mya), N_root = 100,000

Extant population sizes do not affect the data (a single lineage per
species cannot coalesce before the splits) but are set to plausible
values anyway. The mutation rate matches the fixed `mu` in
config/example_config.yaml rather than the literature human rate, since
all iTRAILS estimates are scaled relative to that fixed value.

Regenerate with:  pixi run simulate-example
Output:           data/simulated_alignment.maf (1 Mb in 50 kb MAF blocks)
"""

from pathlib import Path

import msprime
import numpy as np

SEED = 42
SEQ_LEN = 1_000_000
BLOCK_LEN = 50_000
MU = 2e-8    # per site per generation; must match the config's fixed mu
RHO = 1e-8   # per site per generation

# species order: A, B, C, outgroup
SPECIES = ['hg38', 'panTro5', 'gorGor5', 'ponAbe2']

T_AB = 240_000     # human-chimp split (generations)
T_ABC = 350_000    # gorilla split
T_ROOT = 700_000   # orangutan split

N_EXTANT = {'hg38': 15_000, 'panTro5': 35_000, 'gorGor5': 25_000, 'ponAbe2': 25_000}
N_AB = 65_000
N_ABC = 100_000
N_ROOT = 100_000

OUT_PATH = Path(__file__).parent.parent / 'data' / 'simulated_alignment.maf'


def simulate():
    demography = msprime.Demography()
    for species in SPECIES:
        demography.add_population(name=species, initial_size=N_EXTANT[species])
    demography.add_population(name='AB', initial_size=N_AB)
    demography.add_population(name='ABC', initial_size=N_ABC)
    demography.add_population(name='ABCD', initial_size=N_ROOT)
    demography.add_population_split(time=T_AB, derived=['hg38', 'panTro5'], ancestral='AB')
    demography.add_population_split(time=T_ABC, derived=['AB', 'gorGor5'], ancestral='ABC')
    demography.add_population_split(time=T_ROOT, derived=['ABC', 'ponAbe2'], ancestral='ABCD')

    ts = msprime.sim_ancestry(
        samples={species: 1 for species in SPECIES},
        demography=demography,
        ploidy=1,
        sequence_length=SEQ_LEN,
        recombination_rate=RHO,
        random_seed=SEED,
    )
    return msprime.sim_mutations(ts, rate=MU, random_seed=SEED + 1)


def alignment(ts):
    """Return one sequence per species: a shared random ancestral background
    with the simulated alleles written in at variant sites."""
    rng = np.random.default_rng(SEED + 2)
    ancestral = rng.choice(np.array(list('ACGT')), size=SEQ_LEN)
    seqs = {species: ancestral.copy() for species in SPECIES}

    sample_order = [ts.population(ts.node(u).population).metadata['name']
                    for u in ts.samples()]
    n_variants = 0
    for variant in ts.variants():
        pos = int(variant.site.position)
        for species, genotype in zip(sample_order, variant.genotypes):
            seqs[species][pos] = variant.alleles[genotype]
        n_variants += 1
    print(f'{n_variants} variant sites in {SEQ_LEN} bp')
    return seqs


def write_maf(seqs, path):
    with open(path, 'w') as f:
        f.write('##maf version=1 scoring=none\n')
        f.write('# simulated with msprime by scripts/simulate_data.py\n')
        for start in range(0, SEQ_LEN, BLOCK_LEN):
            end = min(start + BLOCK_LEN, SEQ_LEN)
            f.write('\na score=0\n')
            for species in SPECIES:
                block = ''.join(seqs[species][start:end])
                f.write(f's {species}.chr1 {start} {end - start} + {SEQ_LEN} {block}\n')
        f.write('\n')
    print(f'wrote {path}')


if __name__ == '__main__':
    write_maf(alignment(simulate()), OUT_PATH)
