# Inference of admixture and population structure



## PCAngsd and selection

For very recent selection we can look within closely related individuals for example with in Europeans

**data:**

 - Genotype likelihoods in Beagle format
 - ~150k random SNPs with maf > 5%
 - Four EU populations with ~100 individuals in each
 - whole genome sequencing
 - depth 2-9X (1000 genome project)

 ```
CEU | Europeans in Utah (British)
GBR | Great Britain
IBS | Iberian/Spain
TSI | Italien
```

First lets set the paths

```
 # NB this must be done every time you open a new terminal
ThePath=/home/albrechtsen/embo2022/

## copy positions and sample information
cp $ThePath/PCangsd/data/eu1000g.sample.Info .

## load the python module
## PCAngsd
PCANGSD='/home/albrechtsen/pcangsd_vir/bin/pcangsd'

#set path to data
EU1000=$ThePath/PCangsd/data/eu1000g.small.beagle.gz
wc eu1000g.sample.Info
N=424 #one line for header
```


## Explore the data


Take a quick look at the samples data

First try to get an overview of the dataset by copying the information file and making a summary using the following:



```
 # view first lines of sample file
head eu1000g.sample.Info
## cut first column | sort | count
cut -f 2 -d " " eu1000g.sample.Info | sed 1d| sort | uniq -c
```

 - How many samples from each country?

Now let's have a look at the GL file that you have created with ANGSD. It is a "beagle format" file called all.beagle.gz - and will be the input file to PCAangsd.
The first line in this file is a header line and after that it contains a line for each locus with GLs. By using the unix command wc we can count the number of lines in the file:

```
gunzip -c $EU1000 | wc -l
```

 - Use this to find out how many loci there are GLs for in the data set?

Next, to get an idea of what the GL file contains try from the command line to print the first 9 columns of the first 7 lines of the file:


```
gunzip -c $EU1000 | head -n 7 | cut -f1-9 | column -t
```

In general, the first three columns of a beagle file contain marker name and the two alleles, allele1 and allele2, present in the locus (in beagle A=0, C=1, G=2, T=3).

All following columns contain genotype likelihoods (three columns for each individual: first GL for homozygote for allele1,
then GL for heterozygote and then GL for homozygote for allele2). Note that the GL values sum to one per site for each individuals. This is just a normalization of the genotype likelihoods in order to avoid underflow problems in the beagle software it does not mean that they are genotype probabilities.

 - Based on this, what is the most likely genotype of Ind0 in the first locus and the locus six?


## PCAngsd

Run PCangsd with to estimate the covariance matrix while jointly estimating the individuals allele frequencies

```
$PCANGSD -b $EU1000 -o EUsmall -t 10
```

The program estimates the covariance matrix that can then be used for PCA. look at the output from the program

 - The algorithm might only need to low number of PCs to estimate the allele freuqencies. How many significant PCs (see MAP test in output)?

Plot the results in R

```
## R
 cov <- as.matrix(read.table("EUsmall.cov"))

 e<-eigen(cov)
 ID<-read.table("eu1000g.sample.Info",head=T,as.is=F)
 plot(e$vectors[,1:2],col=ID$POP)

 legend("topleft",fill=1:4,levels(ID$POP))
## close R after view plot
```

 - Does the plot look like you expected? Which populations are close and distant to each other?


Since the European individuals in 1000G are not simple homogeneous disjoint populations it is hard to use PBS/FST or similar statistics to infer selection based on populating differences. However, PCA offers a good description of the differences between individuals which out having the define disjoint groups.

Now let try to use the PC to infer selection along the genome based on the PCA

```
$PCANGSD -b $EU1000 -o EUsmall --selection --sites_save --minMaf 0 -t 10
# crate file with position and chromosome
 paste <(zcat $ThePath/PCangsd/data/eu1000g.small.beagle.gz| cut -f 1 | sed 's/\_/\t/g' | sed 1d ) EUsmall.sites  > EUsmall.sites.info
```

View the SNP location info that you will need to plot the results (the third column indicate if the site is used=1 or not =0)

```
head EUsmall.sites.info 
```


plot the results of the selection scan 

```
library(RcppCNPy,lib="/home/albrechtsen/R/x86_64-pc-linux-gnu-library/4.1/") # Numpy library for R

## function for QQplot
qqchi<-function(x,...){
lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
legend("topleft",paste("lambda=",lambda))
}

### read in seleciton statistics (chi2 distributed)
s<-npyLoad("EUsmall.selection.npy")
## make QQ plot to QC the test statistics
qqchi(s)

# convert test statistic to p-value
pval<-1-pchisq(s,1)

## read positions (hg38)
p<-read.delim("EUsmall.sites.info",colC=c("factor","integer","integer"),head=F,as.is=F)

names(p)<-c("chr","pos","keep")

## make manhatten plot
plot(-log10(pval),col=p$chr[p$keep==1],xlab="Chromosomes",main="Manhattan plot")

## zoom into region
 w<-range(which(pval<1e-7)) + c(-100,100)
 keep<-w[1]:w[2]
 plot(p$pos[keep],-log10(pval[keep]),col=p$chr[keep],xlab="HG38 Position chr2")

## see the position of the most significant SNP
 p$pos[which.max(s)]
```

see if you can make sense of the top hit based on the genome.
 - Look in [UCSC browser](http://genome.ucsc.edu/cgi-bin/hgGateway)
 - Choose human GRCh38/hg38
 - search for the position of the top hit and identify the genes at that loci


# Bonus exercises
If you want to try more PCA and admixture proportions from low depth sequencing then try these [Bonus exercises](scanPCA_bonus.md)

