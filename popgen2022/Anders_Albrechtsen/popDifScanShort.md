





# SFS, Fst and PBS

The data is from the 1000 genomes project which included the populations:

```
CEU     | Europeans (mostly of British ancestry)
JPT     | East Asian - Japanese individuals
YRI     | West African - Nigerian Yoruba individuals
```



Due to computation we will use a very reduced data set:

 - Input data: bam files
 - 10 individuals from each population
 - A very reduced genome 30 x 100k random regions across the autosomes + a non-random region
 - Each individual is sequenced at 2-6X



## Aims:

 - To reconstruct the SFS (1D and 2D)
 - To estimate Fst between pairs of popualtions
 - To perform a scan statistics using PBS to detect signs of positive selection


## Set up

First set some paths

```
# NB this must be done every time you open a new terminal
ThePath=/home/albrechtsen/embo2022/

# Set path to ANGSD program
ANGSD=$ThePath/prog/angsd/angsd

#realSFS
REAL=$ThePath/prog/angsd/misc/realSFS

#ancestral fasta file (chimp)
ANC=$ThePath/sfs/data/hg19ancNoChr.fa.gz

#reference genome for human 
REF=$ThePath/sfs/data/hg19.fa.gz

# a bam filelist for a several bam files
BAMFOLDER=$ThePath/sfs/data/smallerbams
BAMFOLDERchr5=$ThePath/sfs/data/chr5_33M_v2

#copy R plot function to folder
cp $ThePath/sfs/plot2dSFS.R .
```

Make some file lists of bam files

```
#a African population
find $BAMFOLDER | grep bam$ | grep YRI > YRI.filelist
#a Asian population
find $BAMFOLDER | grep bam$ | grep JPT > JPT.filelist
#a European population
find $BAMFOLDER | grep bam$ | grep CEU > CEU.filelist
```

                                                                                                 
     

## Reconstructing the site frequency spectrum



First lets set some filter to remove the worst reads (minMapQ), remove the worst of the bases (minQ). 

```
FILTERS="-minMapQ 30 -minQ 20"
```

Lets set some options that means we will calculate genotype likelihoods using the GATK  model (gl) and calculate the site allele frequency likelihoods (saf)
```
OPT=" -dosaf 1 -gl 2"
```

Generate site frequency likelihoods using  ANGSD  
```
$ANGSD -b  YRI.filelist  -anc $ANC -out yri $FILTERS $OPT -ref $REF &
$ANGSD -b  JPT.filelist  -anc $ANC -out jpt $FILTERS $OPT -ref $REF &
$ANGSD -b  CEU.filelist  -anc $ANC -out ceu $FILTERS $OPT -ref $REF
```
The run time is a couple of minutes

If it takes too long then you can copy the results using this command:

```
cp $ThePath/run/yri.saf* .
cp $ThePath/run/ceu.saf* .
cp $ThePath/run/jpt.saf* .
```

Estimate the site frequency spectrum for each of the 3 populations without having to call genotypes or variable sites directly from the site frequency likelihoods

```
#calculate the 1 dimensional SFS
$REAL yri.saf.idx > yri.sfs
$REAL jpt.saf.idx > jpt.sfs
$REAL ceu.saf.idx > ceu.sfs
```


In order to plot the results open R and make a barplot
```
 ##run in R                      
#plot the results
nnorm <- function(x) x/sum(x)
#expected number of sites with 1:20 derived alleles
res <- rbind(
  YRI=scan("yri.sfs")[-1],
  JPI=scan("jpt.sfs")[-1],
  CEU=scan("ceu.sfs")[-1]
)
colnames(res) <- 1:20

# density instead of expected counts
res <- t(apply(res,1,nnorm))


#plot the polymorphic sites. 
resPoly <- t(apply(res[,-20],1,nnorm))
barplot(resPoly,beside=T,legend=c("YRI","JPT","CEU"),names=1:19,main="realSFS polymorphic sites")

#due the very limited amount of sites
#downsample to 5 individuals (10 chromosome) and exclude fixed derived
downsampleSFS <- function(x,chr){ #x 1:2n , chr < 2n
    n<-length(x)
    mat <- sapply(1:chr,function(i) choose(1:n,i)*choose(n- (1:n),chr-i)/choose(n,chr))
    nnorm( as.vector(t(mat) %*% x)[-chr] )
}
resDown <- t(apply(res,1,downsampleSFS,chr=10))
barplot(resDown,beside=T,legend=c("YRI","JPT","CEU"),names=1:9,main="realSFS downsampled polymorphic sites")

```

 - Which population has the largest population size?
 - The analysed whole chromosome for the 1000G individual look [like this](http://popgen.dk/albrecht/phdcourse/sfs/Moltke5V2.pdf)





## Fst and PBS

In order to estimate Fst between two population we will need to estimate the 2-dimensional frequency spectrum from the site allele frequency likelihoods 

```
#calculate the 2D SFS 
$REAL yri.saf.idx ceu.saf.idx >yri.ceu.ml &
$REAL yri.saf.idx jpt.saf.idx >yri.jpt.ml &
$REAL jpt.saf.idx ceu.saf.idx >jpt.ceu.ml
```


In order to get a measure of this populations are most closely related we willl estimate the pairwise Fst from the 2D SFS

```
#first will will index the sample so the same sites are analysed for each population
$REAL fst index jpt.saf.idx ceu.saf.idx -sfs jpt.ceu.ml -fstout jpt.ceu
$REAL fst index yri.saf.idx ceu.saf.idx -sfs yri.ceu.ml -fstout yri.ceu
$REAL fst index yri.saf.idx jpt.saf.idx -sfs yri.jpt.ml -fstout yri.jpt

#get the global estimate
$REAL fst stats jpt.ceu.fst.idx
$REAL fst stats yri.jpt.fst.idx
$REAL fst stats yri.ceu.fst.idx 
```

look at the weigthed Fst (Fst.Weight).
 - which two populations are most closely related?
 - which two populations are most distantly related?	
 


Lets see how the Fst and PBS varies between different regions of the genome my using a sliding windows approach (windows site of 50kb)

```
$REAL fst index yri.saf.idx jpt.saf.idx ceu.saf.idx -fstout yri.jpt.ceu -sfs yri.jpt.ml -sfs yri.ceu.ml -sfs jpt.ceu.ml
$REAL fst stats2 yri.jpt.ceu.fst.idx -win 50000 -step 10000 >slidingwindowBackground
```


read the data into R

```
 ##run in R                      
r<-read.delim("slidingwindowBackground",as.is=T,head=T)
names(r)[-c(1:4)] <- c("wFst_YRI_JPT","wFst_YRI_CEU","wFst_JPT_CEU","PBS_YRI","PBS_JPT","PBS_CEU")


head(r) #print the results to the screen

#plot the distribution of Fst
mmax<-max(c(r$wFst_YRI_JPT,r$wFst_YRI_CEU,r$wFst_JPT_CEU),na.rm=T)
par(mfcol=c(3,2))
hist(r$wFst_YRI_JPT,col="lavender",xlim=c(0,mmax),br=20)
hist(r$wFst_YRI_CEU,col="mistyrose",xlim=c(0,mmax),br=20)
hist(r$wFst_JPT_CEU,col="hotpink",xlim=c(0,mmax),br=20)

mmax<-max(c(r$PBS_CEU,r$PBS_YRI,r$PBS_JPT),na.rm=T)

#plot the distribution of PBS
mmax<-max(c(r$PBS_CEU,r$PBS_YRI,r$PBS_JPT),na.rm=T)
hist(r$PBS_YRI,col="lavender",xlim=c(0,mmax),br=20)
hist(r$PBS_CEU,col="mistyrose",xlim=c(0,mmax),br=20)
hist(r$PBS_JPT,col="hotpink",xlim=c(0,mmax),br=20)


```

note the maximum observed values for both the pairwise fst and the PBS




Lets do the same for not so randomly selection 1Mb region of on chr 5. 
Remember to close R

```                                                                                                                       
#a African population for a region on chr 5                                                 
find $BAMFOLDERchr5 | grep bam$ | grep YRI > YRIchr5.filelist
#a Asian population for a region on chr 5                                                                                       
find $BAMFOLDERchr5 | grep bam$ | grep JPT > JPTchr5.filelist
#a European population for a region on chr 5                                                                                        
find $BAMFOLDERchr5 |  grep bam$ | grep CEU > CEUchr5.filelist

#use the same filters and options as before
FILTERS="-minMapQ 30 -minQ 20 -baq 1 -C 50 -minInd 8"
OPT=" -dosaf 1 -gl 2"

#get site frequency likelihoods
$ANGSD -b  YRIchr5.filelist  -anc $ANC -out yriChr5 $FILTERS $OPT -ref $REF
$ANGSD -b  JPTchr5.filelist  -anc $ANC -out jptChr5 $FILTERS $OPT -ref $REF
$ANGSD -b  CEUchr5.filelist  -anc $ANC -out ceuChr5 $FILTERS $OPT -ref $REF

#estimate the 1D SFS
$REAL yriChr5.saf.idx ceuChr5.saf.idx >yri.ceuChr5.ml
$REAL yriChr5.saf.idx jptChr5.saf.idx >yri.jptChr5.ml
$REAL jptChr5.saf.idx ceuChr5.saf.idx >jpt.ceuChr5.ml

#get FST and PBS in sliding window
$REAL fst index yriChr5.saf.idx jptChr5.saf.idx ceuChr5.saf.idx -fstout yri.jpt.ceuChr5 -sfs yri.jptChr5.ml -sfs yri.ceuChr5.ml -sfs jpt.ceuChr5.ml
$REAL fst stats2 yri.jpt.ceuChr5.fst.idx -win 50000 -step 10000 >slidingwindowChr5

```



Lets view how it looks in this region

```
#run in R
r<-read.delim("slidingwindowChr5",as.is=T,head=T)
names(r)[-c(1:4)] <- c("wFst_YRI_JPT","wFst_YRI_CEU","wFst_JPT_CEU","PBS_YRI","PBS_JPT","PBS_CEU")


par(mfrow=1:2)
plot(r$midPos,r$wFst_YRI_CEU,ylim=c(0,max(r$wFst_YRI_CEU)),type="b",pch=18,ylab="Fst",xlab="position on Chr 5")
points(r$midPos,r$wFst_YRI_JPT,col=2,type="b",pch=18)
points(r$midPos,r$wFst_JPT_CEU,col=3,type="b",pch=18)
legend("topleft",fill=1:3,c("YRI vs. CEU","YRI vs. JPT","JPT vs CEU"))

plot(r$midPos,r$PBS_YRI,ylim=c(0,max(r$PBS_CEU)),type="b",pch=18,ylab="PBS",xlab="position on Chr 5")
points(r$midPos,r$PBS_JPT,col=2,type="b",pch=18)
points(r$midPos,r$PBS_CEU,col=3,type="b",pch=18)
legend("topleft",fill=1:3,c("YRI","JPT","CEU"))
```

 - Compare the values you observed on this part of the genome with the random pars of the genome you looked at earlier [PDF](http://popgen.dk/albrecht/oulu2016/web/PBS.pdf). Is this region extreme?
 - Why is there two peak for the Fst and only one for the PBS?
 - In which of the populations are this loci under selection?


Find out what genes is in this region by going to the [UCSC browser](https://genome.ucsc.edu/index.html). Choose Genome browser. Choose human GRCh37/hg19 and find the region. Read about this gene on wikipedia and see if this fits PBS results. 

