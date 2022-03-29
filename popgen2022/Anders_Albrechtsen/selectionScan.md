# Scan statistics

  # <contents>

## Variability statistics

 




We will use individuals sequencing in the 1000 genomes project. Even thought the individuals we sequenced at low/medium depth they were able to obtain good genotype calls for most of the SNPs in the genome (using imputation + external information). Because of computational demands we will only be looking at a 20Mb region for a small subset of the individuals. We will use 56 African individuals (YRI) and 60 European individuals (CEU). These are the individuals what overlap with the HapMap project where many selection scans in Humans have been performed. We will try to explore the LCT gene located at position 136.6 Mb.  

**Aim**: Locate the LCT selection signal using phased genotype data

We will use the program selscan with contains many of the haplotype based methods used for scan statistics. These methods are based on phased genotpe data.

In the terminal set the paths to the data and program
```console
## path for this exercise
ThePath=/home/albrechtsen/embo2022/

#VCF files
ceuVCF=$ThePath/EHH/data/ceuLCT.recode.vcf
yriVCF=$ThePath/EHH/data/yriLCT.recode.vcf

#genetic map (positions in centimorgans)
MAP=$ThePath/EHH/data/geneticV2.map
```


Look inside one of the VCF files e.g.

```console
less -S $yriVCF
```

 - How can you tell that it is phased data?
 - The VCF file has been filters for certain types of sites. Which ones?



### variability/Tajimas pi

First lets look at the variability of the data in the CEU.

 - If there has been positive selection at the LCT loci what do you expect?

Estimate tajimas theta (pi) in 10k windows

```console
selscan --pi --vcf $ceuVCF --map $MAP --out ceuLCT --threads 8 --pi-win 10000
```

You can plot the results in R e.g.

```console
#in R (open R by typing R in the terminal)
r<-read.table("ceuLCT.pi.10000bp.out",head=F)
r<-subset(r,V3!=0)
plot(r[,1],r[,3],ylab="pi",xlab="position")
causalSNP <- 136608646
abline(v=136.5e6,col="red")
legend("topleft",fill="red","LCT")
## don't close R
```

<!--- If you are having problems with graphic then you can find the plots [[plots][here]] --->

 - Can you see the reduced variability?

To determine whether it is extreme we can compare with the rest of the region
```console
#continue in R
hist(r[,3],br=100)
(causalWin<-subset(r,V1<causalSNP & V2>causalSNP))
abline(v=causalWin[,3],col="red")
```
 

 - is the variability in the LCT region extreme?

 - Try different windows sizes to see if it will change your results (NB! you have to change the name when reading into R)


NB!. This is based on genotypes. I will get back to how to do it for low depth data using genotype likelihoood


[BONUS: If you want to see tajimas D instead]( https://github.com/populationgenetics/exercises/blob/master/NaturalSelection/selectionAA.md )

Use the browser to see how Tajima's D perform on the same data

## Haplotype based ( EHH ) scan statistics

### IHS

Lets see if the haplotype homozygosity does a better job. Run IHS using the command

```console
selscan --ihs --vcf $ceuVCF --map $MAP --out ceuLCT --threads 8
```

The analysis will take a couple of minutes. If you cannot wait then you can copy the results with the command

```console
cp $ThePath/run/ceuLCT.ihs* .
```

The output colums are: 

```
<locusID> <physicalPos> <'1' freq> <ihh1> <ihh0> <unstandardized iHS>
 ```


This statistics will be affected by the frequency of the SNPs therefore we have to normalize in frequency bins. The default in 100 bins 

```console
norm --ihs --files ceuLCT.ihs.out 
```

The number of bins in too high for this data set since we do not have enough SNPs for each bin of allele frequencies. Therefore, redo the analysis where  you  reduce the number of bins to 20 with the - -bins 20 option.

Lets plot the results in R
```console
r <- read.table("ceuLCT.ihs.out.20bins.norm",as.is=T,head=F)
names(r) <- c("locusID", "physicalPos","freq","ihh1","ihh0","unstandardizediHS")
r[which.max(r$ihh1/r$ihh0),]
causalSNP <- 136608646
#plot without frequency standardization
plot(r$physicalPos,r$unstandardizediHS);

## with standardiztion IHS=ihh0/ihh1
r$iHS<-log(r$ihh1/r$ihh0)
plot(r$physicalPos,r$iHS,ylab="iHs");
abline(v=causalSNP,col="red")

## causal SNP test statistics vs. rest of region
(causalSite<-r[which(r$physicalPos==causalSNP),])
hist(r$iHS)
abline(v=causalSite$iHS,col="red")
# close R when you are done
```

Lets try to use the West Africans (YRI) to normalize the IHS using XpEHH

```console
selscan --xpehh --vcf $ceuVCF --map $MAP --vcf-ref $yriVCF   --out ceuLCT --threads 8
```
It will take about 5 mins and you can choose to copy the results instead

```console
#copy results instaed of waiting                                  
cp $ThePath/run/ceuLCT.xpehh* .
```


We also have to normalize this

```console
norm --xpehh --files ceuLCT.xpehh.out
```

again you can plot the results in R.

```console
r<-read.table("ceuLCT.xpehh.out.norm",head=T,as.is=T,row.names=NULL)
causalSNP <- 136608646
(causalSite<-r[which(r$pos==causalSNP),])                 
plot(r$pos,r$normxpehh)
abline(v=causalSNP,col="red")

#print the site with maximum statistic 
r[which.max(r$normxpehh),]

#plot the distribution
hist(r$normxpehh)
abline(v=causalSite$normxpehh,col="red")

#get the quantile
mean(causalSite$normxpehh>r$normxpehh)

```

 - Are you more convinced that the site is under selection?

