
import msprime

# Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.
ts = msprime.sim_ancestry(
    samples=10,
    recombination_rate=1e-8, # as in humans
    sequence_length=10_000,
    population_size=10_000,
    random_seed=1234)

# we can add mutations
mutated_ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=1234)
print(mutated_ts.tables.sites)

for variant in mutated_ts.variants():
    print(variant)

