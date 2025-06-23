library(ggplot2);

rm(list=ls());
# First Session 14:00-15:30 Practical session: Maths & Stats practical

# Exercise 1. 
# a. Given an effective population size Ne, and a given frequency of the derived allele for a SNP,
# create a forward simulator that works under neutrality and under selection

next.gen <- function(Ne, f)
{
  coin <- rep(0,Ne);
  coin[sample(1:Ne,rbinom(1,Ne,f))] <- 1;
  return(coin);
}

next.gen.sel <- function(Ne, f, s)
{
  coin <- rep(0,Ne);
  coin[sample(1:Ne,rbinom(1,Ne,f))] <- 1;
  id.0 <- which(coin==0);
  id.0.to.change <- sample(id.0,sum(rbinom(length(id.0),1,s)));
  coin[id.0.to.change] <- 1;
  return(coin);
}

Ne = 1000;
pop <- next.gen(Ne, 1/Ne);
f <- mean(pop);
gen <- 1;
record.f <- c(1/Ne);
s <- 0.01

while(gen < 500)
{
  # fixation. Re-start
  if(f==0)
  {
    gen = 1;
    pop <- first.gen(Ne);
    record.f <- c(1/Ne);
  }
  else
  {
    pop <- next.gen.sel(Ne,f,s);
    gen = gen + 1;
  }
  f = mean(pop);
  record.f <- c(record.f,f);
}

plot(record.f,type="l",xlab = "number of generations since mutation starts", ylab = "frequency in pop");

# b. Simulate the TMRCA from two sampled chromosomes in a panmictic population with a given Ne. Estimate it for K independent loci

# A very silly way of doing it.
tmrca <- function(Ne, l)
{
  tmrca <- rep(NA,l);
  for(i in 1:l)
  {
    g = 1;
    while(sample(Ne,1)!=sample(Ne,1))
    {
      g = g+1;
    }
    tmrca[i] <- g;
  }
  return(tmrca);
}

# run it with Ne = 100 and Ne = 1000 on 5000 loci. Plot the log(TMRCA) of both. What is happening?

ne100 <- log(tmrca(100,5000));
ne1000 <- log(tmrca(1000,5000));
df <- data.frame(logTMRCA = c(ne100,ne1000), model = c(rep("ne100",length(ne100)),rep("ne1000",length(ne1000))))

ggplot(df, aes(x=logTMRCA, color=model)) +
  geom_density()

rm(list=ls());

folder.fastSimcoal2 <- "/home/lao/";
exe <- "fsc28"
args <- c("-i", paste(folder.fastSimcoal2,"DemographicModelSplitR.par", sep=""), "-x", "-s0", "-d", "-n", "1", "-q", "-G")

setwd(folder.fastSimcoal2);

# Exercise 2)
# Produce the par file for fastSimcoal2. Generate number_of_blocks 
# independent genomic regions, each of 1 Mb. Assume a mutation rate of 1.6*10^-7 and recombination rate of 10^-8.
# a) assuming a constant population size A. Sample 100 chromosomes

model.single.pop <- function(effective_population_size_1, number_of_blocks)
{
  lines <- c(
    "//Number of population samples (demes)",
    "1",
    "//Population effective sizes (number of genes)",
    effective_population_size_1,
    "//Sample sizes",
    "100",
    "//Growth rates\t: negative growth implies population expansion",
    "0",
    "//Number of migration matrices : 0 implies no migration between demes",
    "0",
    "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
    "0 historical event",
    "//Number of independent loci [chromosome]",
    paste(number_of_blocks,"1",sep=" ")
  )
  
  # Repeating the block 22 times
  block <- c(
    "//Per chromosome: Number of linkage blocks",
    "1",
    "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
    "DNA 1000000 1.0E-8 1.855284327902964E-7 0.0"
  )
  
  # Append number_of_blocks blocks
  for (i in 1:number_of_blocks) {
    lines <- c(lines, block)
  }
  
  # Write to file
  writeLines(lines, "DemographicModelSplitR.par")
  # Execute
  system2(exe, args = args);
  
  data.t <- read.table(file=paste(folder.fastSimcoal2,"DemographicModelSplitR/DemographicModelSplitR_1_1.gen", sep=""), header = T);
  
  # First four columns are snp info
  # haplotype matrix. Rows are haplotypes, columns are positions
  H <- t(as.matrix(data.t[,5:ncol(data.t)]));
  rownames(H) <- rep("A",100);
  # return a list with the chromosomal positions and the haplotype matrix
  return(list(position = data.t[,1:2],haplotype_matrix = H));
}

ne <-1000
d <- model.single.pop(ne, 5)$haplotype_matrix;


# b) assuming a constant population size A with a recent population change
model.single.pop.with.recent.population.change <- function(population_size_in_present, effective_population_size_1, time_change_population_size, number_of_blocks)
{
  r <- (effective_population_size_1/population_size_in_present);
  lines <- c(
    "//Number of population samples (demes)",
    "1",
    "//Population effective sizes (number of genes)",
    effective_population_size_1,
    "//Sample sizes",
    "100",
    "//Growth rates\t: negative growth implies population expansion",
    "0",
    "//Number of migration matrices : 0 implies no migration between demes",
    "0",
    "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
    "1 historical event",
    paste(time_change_population_size, "0 0 1", r,"0 0",sep = " "),
    "//Number of independent loci [chromosome]",
    paste(number_of_blocks,"1",sep=" ")
  )
  
  # Repeating the block 22 times
  block <- c(
    "//Per chromosome: Number of linkage blocks",
    "1",
    "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
    "DNA 1000000 1.0E-8 1.855284327902964E-8 0.0"
  )
  
  # Append number_of_blocks blocks
  for (i in 1:number_of_blocks) {
    lines <- c(lines, block)
  }
  
  # Write to file
  writeLines(lines, "DemographicModelSplitR.par")
  # Execute
  system2(exe, args = args);
  
  data.t <- read.table(file=paste(folder.fastSimcoal2,"DemographicModelSplitR/DemographicModelSplitR_1_1.gen", sep=""), header = T);
  
  # First four columns are snp info
  # haplotype matrix. Rows are haplotypes, columns are positions
  H <- t(as.matrix(data.t[,5:ncol(data.t)]));
  rownames(H) <- rep("A",100);
  # return a list with the chromosomal positions and the haplotype matrix
  return(list(position = data.t[,1:2],haplotype_matrix = H));
}

ne <-1000
ne.past <- 10000;
time.change <- 10;
d <- model.single.pop.with.recent.population.change(ne, ne.past, time.change, 5)$haplotype_matrix;

# c) The demographic model has two populations (A and B).
# Population A has 20000 chromosomes. Population B has 1000 chromosomes. Both populations split
# 1000 generations ago. Sample 100 chromosomes from A and 100 chromosomes from B.

model.two.pops <- function(effective_population_size_1, effective_population_size_2, time_split, number_of_blocks)
{
  lines <- c(
    "//Number of population samples (demes)",
    "2",
    "//Population effective sizes (number of genes)",
    effective_population_size_1,
    effective_population_size_2,
    "//Sample sizes",
    "100",
    "100",
    "//Growth rates\t: negative growth implies population expansion",
    "0",
    "0",
    "//Number of migration matrices : 0 implies no migration between demes",
    "0",
    "//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix",
    "1 historical event",
    paste(time_split,"0 1 1 1 0 0",sep = " "),
    "//Number of independent loci [chromosome]",
    paste(number_of_blocks,"1",sep=" ")
  )
  
  # Repeating the block 22 times
  block <- c(
    "//Per chromosome: Number of linkage blocks",
    "1",
    "//per Block: data type, num loci, rec. rate and mut rate + optional parameters",
    "DNA 1000000 1.0E-8 1.855284327902964E-8 0.0"
  )
  
  # Append number_of_blocks blocks
  for (i in 1:number_of_blocks) {
    lines <- c(lines, block)
  }
  
  # Write to file
  writeLines(lines, "DemographicModelSplitR.par")
  
  # Execute
  system2(exe, args = args)
  
  data.t <- read.table(file=paste(folder.fastSimcoal2,"DemographicModelSplitR/DemographicModelSplitR_1_1.gen", sep=""), header = T);
  
  # First four columns are snp info
  # haplotype matrix. Rows are haplotypes, columns are positions
  
  H <- t(as.matrix(data.t[,5:ncol(data.t)]));
  # First 100 haplotypes is pop 1, next pop 2
  rownames(H) <- c(rep("A",100),rep("B",100));
  return(list(position = data.t[,1:2],haplotype_matrix = H));
}

# Exercise 3. Very basic single biallelic marker summary statistics
# 1) Implement for a given marker a classical Fst implementation.

classical_fst <- function(s)
{
  a <- s[1:length(s)/2];
  b <- s[(1+length(s)/2):length(s)];
  pa <- mean(a);
  pb <- mean(b);
  p <- (pa+pb)/2;
  fst <- ((pa-p)^2 + (pb-p)^2)/(4*p*(1-p));
  return(fst);
}

# 2) Implement for a given marker the informativeness of ancestry.
rosenberg_informativeness_of_ancestry <- function(s)
{
  info <- 0;
  a <- s[1:length(s)/2];
  b <- s[(1+length(s)/2):length(s)];
  for(al in 1:2)
  {
    if(al==2)
    {
      a <- 1-a;
      b <- 1-b;
    }
    pa <- mean(a);
    pb <- mean(b);
    p <- (pa+pb)/2;
    if(pa==0)
    {
      pa = 1;
    }
    if(pb==0)
    {
      pb = 1;
    }
    info <- info + -p*log(p) + 0.5*(pa*log(pa) + pb*log(pb));    
  }
  return(info);
}

# 3) Estimate the SFS of a single population

sfs <- function(m)
{
  alleles.by.snp <- colSums(m);
  min <- 1;
  max <- nrow(m) - 1;
  return(tabulate(factor(alleles.by.snp, levels = min:max))/ncol(m));
}

# 5) Estimate the HWE for each marker. Make a simulation with the single population.
# compute the genotypes and estimate for each genotype the pvalue

HWE <- function(genotype)
{
  n <- length(genotype);
  p <- mean(genotype)/2;
  G <- matrix(nrow=3,ncol=2,0);
  G[1,1] <- round(n*p^2);
  G[1,2] <- sum(genotype==2);
  G[2,1] <- round(n*2*p*(1-p));
  G[2,2] <- sum(genotype==1);
  G[3,1] <- round(n*(1-p)^2);
  G[3,2] <- sum(genotype==0);  
  return(log(fisher.test(G)$p.value));
}

#4) Implement a pipeline that
# a) considers the model of a single population, with ne changing from 100 to 200000 chromosomes at random.
# estimate the mean distance between all the haplotypes. Use five blocks of 1Mb. Repeat 100 times. Plot the Ne vs the mean heterozygosity.

ne.parameter <- rep(0,100);
heterozygosity.stat <- rep(0,100);
for(rep in 1:100)
{
  ne <- round(runif(1,100,20000));
  d <- model.single.pop(ne, 5)$haplotype_matrix;
  n <- nrow(d);
  nsnps <- ncol(d);
  distance <- sum(as.matrix(dist(d)))/(n*(n-1));
  heterozygosity.stat[rep] <- distance;
  ne.parameter[rep] <- ne;
}

plot(ne.parameter,heterozygosity.stat, xlab = "ne", ylab = "Mean Heterozygosity");

# b) considers that popA and popB has Ne = 1000
# with a time of split between the two populations that range between 10 to 20000 generations.
# For each simulated dataset, compute the mean Fst among all the SNPs.
# Repeat 1000 simulations. Plot the time of split vs the Fst.

effective_population_size_1 <- 1000;
effective_population_size_2 <- 1000;
time.parameter <- rep(0,100);
fst.stat <- rep(0,100);
for(rep in 1:100)
{
  time_split <- round(runif(1,10,20000));
  d <- model.two.pops(effective_population_size_1, effective_population_size_2, time_split,2)$haplotype_matrix;
  fs <- apply(d,2,rosenberg_informativeness_of_ancestry);
  time.parameter[rep] <- time_split;
  fst.stat[rep] <- mean(fs);
}

plot(time.parameter,fst.stat, xlab = "time of split", ylab = "Mean Genetic distance")


# run with one population. Make a histogram. What is happening?

ne <- round(runif(1,100,20000));
d <- model.single.pop(ne, 5)$haplotype_matrix;

genotypes <- d[seq(from=1,to=nrow(d),by=2),] + d[seq(from=2,to=nrow(d),by=2),];
hwe <- apply(genotypes,2,HWE);
length(hwe)
hist(hwe)

# repeat but with two populations. Make a histogram. What is happening?

effective_population_size_1 <- 1000;
effective_population_size_2 <- 1000;
time_split <- round(runif(1,1000,20000));
d <- model.two.pops(effective_population_size_1, effective_population_size_2, time_split,2)$haplotype_matrix;
genotypes <- d[seq(from=1,to=nrow(d),by=2),] + d[seq(from=2,to=nrow(d),by=2),];
hwe.two.pop <- apply(genotypes,2,HWE);
hist(hwe.two.pop)

#################################################################
#
# HWE in two pops changing the time of split
#
#
#################################################################

t.split <- rep(0,100);
proportion.out.HWE <- rep(0,100);

for(rep in 1:100)
{
  effective_population_size_1 <- 1000;
  effective_population_size_2 <- 1000;
  time_split <- round(runif(1,1,20000));
  d <- model.two.pops(effective_population_size_1, effective_population_size_2, time_split,2)$haplotype_matrix;
  genotypes <- d[seq(from=1,to=nrow(d),by=2),] + d[seq(from=2,to=nrow(d),by=2),];
  hwe.two.pop <- apply(genotypes,2,HWE);
  t.split[rep] <- time_split; 
  proportion.out.HWE[rep] <- mean(hwe.two.pop < log(0.05));
}

plot(t.split, proportion.out.HWE);


# Compute the SFS of a single population with Ne = 10000

ne <- 10000;
d <- model.single.pop(ne, 5)$haplotype_matrix;
sfs.single.pop <- sfs(d);

plot(sfs.single.pop);

# Compute the SFS of a single population with Ne = 10000, and a recent bottleneck with Ne = 1000, 100 generations ago

population_size_in_present <- 1000;
effective_population_size_1 <- 10000;
time_change_population_size <- 100;
number_of_blocks <- 5;

d <- model.single.pop.with.recent.population.change(population_size_in_present, effective_population_size_1, time_change_population_size, number_of_blocks)$haplotype_matrix;
sfs.single.pop.recent.change <- sfs(d);

# Compare both. What is happening?

ma <- data.frame(counts = c(seq(from=1,to = (nrow(d)-1),by=1),seq(from=1,to = (nrow(d)-1),by=1)), sfs = c(sfs.single.pop, sfs.single.pop.recent.change), model = c(rep("constant"),rep("recent")))

ggplot(ma, aes(x=counts, y=sfs, color=model)) +
  geom_point()

# Repeat using a time of 10 generations. Compare both

population_size_in_present <- 1000;
effective_population_size_1 <- 10000;
time_change_population_size <- 10;
number_of_blocks <- 5;

d <- model.single.pop.with.recent.population.change(population_size_in_present, effective_population_size_1, time_change_population_size, number_of_blocks)$haplotype_matrix;
sfs.single.pop.recent.change <- sfs(d);
# Compare both. What is happening?

ma <- data.frame(counts = c(seq(from=1,to = (nrow(d)-1),by=1),seq(from=1,to = (nrow(d)-1),by=1)), sfs = c(sfs.single.pop, sfs.single.pop.recent.change), model = c(rep("constant"),rep("recent")))

ggplot(ma, aes(x=counts, y=sfs, color=model)) +
  geom_point()

# LD

maf <- function(h)
{
  f <- mean(h);
  return(min(f,1-f));
}

ld <- function(a,b)
{
  return(cor(a,b)^2);
}

ne <- 1000;
sim <- model.single.pop(ne, 1);
d <- sim$haplotype_matrix;
select.maf <- apply(d,2,maf) > 0.1;
d <- d[,select.maf];
po <- sim$position[select.maf];
ld.m <- c();
di.m <- c();
ma.m <- c();
for(s1 in 1:(nrow(d)-1))
{
  m1 <- maf(d[,s1]);
  for(s2 in (s1+1):nrow(d))
  {
    di <- abs(po[s1,2] - po[s2,2]);
    m2 <- maf(d[,s2]);
    di.m <- c(di.m,di);
    ld.m <- c(ld.m, ld(d[,s1],d[,s2]));
    ma.m <- c(ma.m, abs(m1-m2));
  }
}

plot(log(di.m), ld.m);



#######################################################
#
# Working with real files.
# Plink files. BED, BIM, FAM
#
#######################################################

library(BEDMatrix)

sim <- model.single.pop(ne, 5);
# first 5 columns are information of the SNPs
haplotypes <- sim$haplotype_matrix;
# Genotypes are computed by collapsing two contiguous haplotypes
genotypes <- haplotypes[seq(from=1,to=nrow(haplotypes),by=2),] + haplotypes[seq(from=2,to=nrow(haplotypes),by=2),];
# Provide labels
pop <- rownames(haplotypes)[seq(from=1,to=nrow(haplotypes),by=2)];
rownames(genotypes) <- pop;
# Compute the standard deviation of each marker to decide if it is worth being included
sd.snps <- apply(genotypes,2,sd);
consider.snps <- which(sd.snps > 0.1);
genotypes.clean <- genotypes[,consider.snps];

# store the genotype matrix in Plink format
bim.data <- data.frame(sim$position[consider.snps,1], paste("rs",sim$position[consider.snps,2],sep=""), rep(0,length(consider.snps)), sim$position[consider.snps,2], sim$position[consider.snps,3],sim$position[consider.snps,4]);
colnames(bim.data) <- c("chr", "id", "posg", "pos", "alt", "ref");
fam.data <- data.frame(pop,paste(pop,seq(from=1,to=length(pop),by=1),sep="_"),rep(0,length(pop)),rep(0,length(pop)),rep(0,length(pop)),rep(-9,length(pop)));
colnames(fam.data) <- c("fam", "id", "pat", "mat", "sex", "pheno");
bed.matrix <- t(genotypes.clean);
colnames(bed.matrix) <- fam.data$id;
write_plink("model_single_pop", bed.matrix, bim = bim.data, fam = fam.data);
# load BED
# Lets check that the dataset was properly written
x <- BEDMatrix("model_single_pop");
# Mean minor allele frequency (MAF) by SNP
maf <- colMeans(X, na.rm = TRUE) / 2
hist(maf, main = "Minor Allele Frequency Distribution", xlab = "MAF")

