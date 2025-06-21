import msprime
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def repeat_simulations(mut, sample_sizes, length, reco, pop_size, num_simulations, seed=None):
    results = []
    for i in tqdm(range(num_simulations), desc="Running simulations"): 
        if seed is not None:
            np.random.seed(seed + i) 
        # Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.
        ts = msprime.sim_ancestry(
            samples=sum(sample_sizes),
            recombination_rate=reco,
            sequence_length=length,
            population_size=pop_size,
            random_seed=np.random.randint(99999999))
        
        # we can add mutations
        mutated_ts = msprime.sim_mutations(ts, rate=mut, random_seed=np.random.randint(99999999))

        diversity = mutated_ts.diversity()
        tajimas_d = mutated_ts.Tajimas_D()
        allele_frequency_spectrum = mutated_ts.allele_frequency_spectrum(polarised=True)
        results.append((mutated_ts, None, diversity, tajimas_d, allele_frequency_spectrum))
    return results

mut= 1e-8
sample_sizes = [20]
length = 10_000
seed = 4711
reco = 1e-8
pop_size = 10_000
num_simulations = 100

results = repeat_simulations(mut, sample_sizes, length, reco, pop_size, num_simulations, seed=seed)


# plot 
diversities = [result[2] for result in results]
tajimas_ds = [result[3] for result in results]
allele_frequency_spectra = [result[4] for result in results]


plt.clf()
plt.figure(figsize=(10, 5))
plt.hist(diversities, bins=10, color='skyblue', edgecolor='black', alpha=0.7)
plt.xlabel("Nucleotide Diversity (π)")
plt.ylabel("Frequency")
plt.title("Histogram of Nucleotide Diversity Across Simulations")
plt.show()


plt.clf()
plt.figure(figsize=(10, 5))
plt.hist(tajimas_ds, bins=10, color='pink', edgecolor='black', alpha=0.7)
plt.xlabel("Tajima's D")
plt.ylabel("Frequency")
plt.title("Distribution of Tajima's D Across Simulations")
plt.show()


plt.clf()
plt.figure(figsize=(10, 5))
bar_width = 0.8 / num_simulations 
colors = plt.cm.tab20(np.linspace(0, 1, num_simulations)) 
for i, afs in enumerate(allele_frequency_spectra):
    x_positions = np.arange(len(afs)) + i * bar_width  
    plt.bar(x_positions, afs, width=bar_width, color=colors[i], label=f'Simulation {i+1}')

plt.xlabel("Frequency")
plt.ylabel("Number of Sites")
plt.title("Allele Frequency Spectrum Across Simulations")
plt.show()

combined_afs = np.sum(allele_frequency_spectra, axis=0)
normalized_afs = combined_afs / np.sum(combined_afs)

plt.clf()
plt.figure(figsize=(10, 5))
plt.bar(range(len(normalized_afs)), normalized_afs, color='green', edgecolor='black', alpha=0.7)
plt.xlabel("Derived Allele Frequency")
plt.ylabel("Proportion of Sites")
plt.title("Allele Frequency Spectrum Across Simulations")
plt.show()

# Summary stats
print(f"Nucleotide Diversity: mean = {np.mean(diversities)}, std = {np.std(diversities)}")
print(f"Tajima's D: mean = {np.mean(tajimas_ds)}, std = {np.std(tajimas_ds)}")
