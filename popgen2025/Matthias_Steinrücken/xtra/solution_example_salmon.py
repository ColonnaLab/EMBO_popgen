import msprime
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def fish(params, sample_sizes, length, seed, reco, generation_time=6):
    assert len(sample_sizes) == 1

    demographic_events = []
    mutation_rate_per_year = params.get("mut") / generation_time
    recombination_rate_per_year = reco / generation_time

    #initial population size
    demographic_events.append(msprime.PopulationParametersChange(time=0, initial_size=100_000, growth_rate=0))

    def add_demographic_event(start_size, end_size, start_time, end_time):
        rate = np.log(end_size / start_size) / (end_time - start_time)
        for t in range(start_time, end_time + 1):
            population_size = int(start_size * np.exp(rate * (t - start_time)))
            demographic_events.append(msprime.PopulationParametersChange(time=t / generation_time, initial_size=population_size))

    add_demographic_event(100_000, 800_000, 25, 45)

    for t in range(46, 50):
        demographic_events.append(msprime.PopulationParametersChange(time=t / generation_time, initial_size=800_000))

    ts = msprime.simulate(
        sample_size=sum(sample_sizes),
        demographic_events=demographic_events,
        mutation_rate=params.get("mut"),
        length=length,
        recombination_rate=reco,
        random_seed=seed
    )

    return ts, demographic_events

def save_results_to_csv(results, filename):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Simulation', 'Diversity', 'Tajimas_D', 'Allele_Frequency_Spectrum'])
        for i, result in enumerate(results):
            writer.writerow([i+1, result[2], result[3], result[4]])

def repeat_simulations(params, sample_sizes, length, reco, generation_time, num_simulations, seed=None):
    results = []
    for i in tqdm(range(num_simulations), desc="Running simulations"): 
        if seed is not None:
            np.random.seed(seed + i) 
        ts, demographic_events = fish(params, sample_sizes, length, seed, reco, generation_time)
        diversity = ts.diversity()
        tajimas_d = ts.Tajimas_D()
        allele_frequency_spectrum = ts.allele_frequency_spectrum(polarised=True)
        results.append((ts, demographic_events, diversity, tajimas_d, allele_frequency_spectrum))
    return results

params = {"mut": 1e-8}
sample_sizes = [20]
length = 100_000
seed = None
reco = 1e-8
generation_time = 6
num_simulations = 100

results = repeat_simulations(params, sample_sizes, length, reco, generation_time, num_simulations, seed=seed)

times = []
sizes = []
for event in results[0][1]:
    times.append(event.time * generation_time)
    sizes.append(event.initial_size)

times = np.array(times)
sizes = np.array(sizes)


plt.clf()
plt.plot(times, sizes)
plt.xlabel("Time (years)")
plt.ylabel("Population Size")
plt.title("Population Size Changes Over Time")
plt.gca().invert_xaxis()
plt.show()

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
