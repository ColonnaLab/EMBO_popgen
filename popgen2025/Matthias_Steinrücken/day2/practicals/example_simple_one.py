
import msprime

# Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.

# 100_000 individuals at present
# 10_000 individuals before 1000 generations ago
demography = msprime.Demography()
demography.add_population(name="A", initial_size=100_000)
# instantaneous reduction of size to 10_000
demography.add_population_parameters_change(population="A", time=1000, initial_size=10_000)


ts = msprime.sim_ancestry(
    samples=10,
    recombination_rate=1e-8, # as in humans
    sequence_length=10_000,
    demography=demography,
    random_seed=1234)

# we can add mutations
mutated_ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=1234)
print(mutated_ts.tables.sites)

for variant in mutated_ts.variants():
    print(variant)

