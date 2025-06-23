
import os
import csv
import msprime
import tskit
import numpy as np
import tqdm

# Function for simulating data under an IM model with parameters:
# Nanc, T_split, N1, N2, mig

def im(params, sample_sizes, seed):
    """Simulate data for 2 populations."""
    assert len(sample_sizes) == 2

    # Extract parameters
    N1 = params.get("N1")
    N2 = params.get("N2")
    T_split = params.get("T_split")
    N_anc = params.get("N_anc")

    # Define population configurations
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=sample_sizes[0], initial_size=N1),
        msprime.PopulationConfiguration(sample_size=sample_sizes[1], initial_size=N2)
    ]

    # Define migration events
    mig = params.get("mig")
    mig_time = T_split / 2  # no migration initially
    if mig >= 0:            # directional (pulse)
        mig_event = msprime.MassMigration(time=mig_time, source=1, destination=0, proportion=abs(mig)) # migration from pop 1 into pop 0 (back in time)
    else:
        mig_event = msprime.MassMigration(time=mig_time, source=0, destination=1, proportion=abs(mig)) # migration from pop 0 into pop 1 (back in time)

    # Define demographic events
    demographic_events = [
        mig_event,
        msprime.MassMigration(time=T_split, source=1, destination=0, proportion=1.0), # move all in deme 1 to deme 0
        msprime.PopulationParametersChange(time=T_split, initial_size=N_anc, population_id=0) # change to ancestral size
    ]

    # Simulate tree sequence
    ts = msprime.simulate(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        mutation_rate=params.get("mut"),
        length=params.get("length"),
        recombination_rate=params["reco"],
        random_seed=seed
    )

    return ts

# Define some initial parameters
params = {
    "mut": 3.5e-9,      # Mutation rate, fixed
    "length": 1e4,      # Sequence length, fixed
    "reco": 8.4e-9,     # recombination rate, fixed
}

sample_sizes = [50, 50]  # Sample sizes for two populations
seed = None               # Random seed

# Output directory
output_directory = "."
# Output file name
output_file = os.path.join(output_directory, "mosquito-task2.csv")


# Perform simulations
numSegSites = []
for i in tqdm.tqdm(range(100)):

    # normally distributed around 8000 with std dev 2000, minimmum is 1000
    params["T_split"] = max (1000,int(np.random.normal(loc=8000, scale=2000, size=1)[0]))

    # normally distributed around 150,000 with std dev 10,000, minimum 5,000
    params["N1"] = max (5000,int(np.random.normal(loc=8000, scale=2000, size=1)[0]))

    # normally distributed around N1/30 with std dev 333, minimum 160
    params["N1"] = max (params["N1"]/30,int(np.random.normal(loc=8000, scale=2000, size=1)[0]))

    # coin flip, either 0 or 0.1, both with prob 0.5
    coin = np.random.randint(2)
    if (coin == 0):
        params["mig"] = 0
    else:
        params["mig"] = 0.1

    # normally distributed around 7,000,000 with std dev 100, minimum 1000
    params["N_anc"] = max (1000,int(np.random.normal(loc=7_000_000, scale=100, size=1)[0]))

    # do the simulations
    ts = im(params, sample_sizes, seed)

    # store the numSegSites from this run
    ssites = ts.segregating_sites (sample_sets=[ts.samples(population=0), ts.samples(population=1)])
    numSegSites.append (ssites[0])

print(np.mean(numSegSites))
