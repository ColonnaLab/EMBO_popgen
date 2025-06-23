
import os
import csv
import msprime
import tskit
import numpy as np
import tqdm

# Function for simulating data under an IM model with parameters:
# Nanc, T_split, N1, N2, mig

def im(params, sample_sizes, seed, reco):
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
        recombination_rate=reco,
        random_seed=seed
    )

    return ts

# Define some initial parameters
params = {
    "N1": 150000,    # Population 1 size 
    "N2": 5000,     # Population 2 size 
    "T_split": 5000,    # Time of population split
    "N_anc": 7148911,   # Ancestral population size (7,148,911)
    "mut": 3.5e-9,      # Mutation rate, fixed
    "length": 1e4,      # Sequence length, fixed
    "reco": 8.4e-9,     # recombination rate, fixed
    "mig": 0            # migration rate, fixed
}

sample_sizes = [50, 50]  # Sample sizes for two populations
seed = None               # Random seed

# Output directory
output_directory = "."
# Output file name
output_file = os.path.join(output_directory, "mosquito-task2.csv")

# Open the output file in write mode
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter=",")

    # Write header
    writer.writerow(["N1", "N2", "T_split", "MigRate", "Fst", "dxy", "segsites1", "segsites2", "pi1", "pi2", "tajima1", "tajima2"])

    # Perform simulations
    for i in tqdm.tqdm(range(1000)):

        params["N1"] = 150_000
        params["N2"] = 5_000
        params["mig"] = 0
        params["T_split"] = max(1000,int(np.random.normal(loc=8000, scale=2000, size=1)[0]))
                        
        ts = im(params, sample_sizes, seed, params["reco"])

        dxy = ts.divergence(sample_sets=[ts.samples(population=0), ts.samples(population=1)])

        Fst = ts.Fst(sample_sets=[ts.samples(population=0), ts.samples(population=1)])

        ssites = ts.segregating_sites(sample_sets=[ts.samples(population=0), ts.samples(population=1)])
        div = ts.diversity(sample_sets=[ts.samples(population=0), ts.samples(population=1)])
        tajima = ts.Tajimas_D(sample_sets=[ts.samples(population=0), ts.samples(population=1)])
    
        # Write data to file or print data
        writer.writerow([params["N1"], params["N2"], params["T_split"], params["mig"], Fst, dxy, ssites[0], ssites[1], div[0], div[1], tajima[0], tajima[1]])


