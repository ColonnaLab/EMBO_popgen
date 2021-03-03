
##########################################
# practical exercises f3 / f4 / f4 ratio #
##########################################

## In this practical we will explore some basic f3 and f4/D tests on an examples dataset of modern and ancient humans from the literature (Human origins array and selected ancient individuals)

## --------------------------------------------------------------------------------
## examine the dataset

## number of individuals and SNPs
wc -l ho_anc.indiv
wc -l ho_anc.snp

## populations and sample sizes
awk '{print $3}' ho_anc.indiv | sort | uniq -c
    

## --------------------------------------------------------------------------------
## admixture f3 statistics - modern source populations

## examine par file
cat master.f3.par


## African - American admixture
qp3Pop -p aa.f3.par | tee aa.f3.log

## Q: Do we have evidence for admixture?
## A: 


## --------------------------------------------------------------------------------
## admixture f3 statistics - ancient source populations

## Admixture in two East African target populations (Somali from Somalia, Dinka from Sudan); Source populations are Mota (~4,000 year old individual from Ethiopia) and different modern and ancient (LBK) Eurasian populations

qp3Pop -p eafr.f3.par | tee eafr.f3.log

## Q: which combinations of source / target populations show evidence for admixture?
## A:

## Q: based on the values of the f3 statistics, which source population is the best proxy for the true admixing population?
## A:

## Admixture in modern Europeans from three possible sources: Loschbour (Mesolithic hunter-gatherer); LBK (Early Neolithic farmers); Yamnaya (Steppe pastoralists)
    
qp3Pop -p eur.f3.par | tee eur.f3.log

## Q: which populations show neolithic farmer / hunter-gatherer admixture?
## A:

## Q: which populations have evidence for Yamnaya admixture
## A:

## outgroup f3
qp3Pop -p outgroup.f3.par | tee outgroup.f3.log 

## Q: which population shares the most drift with LBK?
grep result: outgroup.f3.log | sort -rnk5,5
## A:


## --------------------------------------------------------------------------------
## f4/D statistics

## examine par file
cat master.f4.par

## European ancient admixture, here we are testing whether a modern European test population (p3) forms a clade with LBK to the exclusion of Yamnaya 
qpDstat -p eur.f4.par | tee eur.f4.log

## Q: which configuration is consistent with a simple tree?
## A:

## Q: how can we interpret the configurations that fail the test?
## A:

## AA admixture
qpDstat -p aa.f4.par | tee aa.f4.log

## Q: Do we have evidence for admixture?
## A:

## Q: Which of the two source population (French, Yoruba) shares more drift with African Americans?
## A:

## Q: Can we conclude that African Americans have higher proportions of European admixture
## A:


## --------------------------------------------------------------------------------
## f4 ratio estimation

## examine par file
cat master.f4ratio.par

## AA admixture
qpF4ratio -p aa.f4ratio.par | tee aa.f4ratio.log

## Q: what are the admixture proportions for the African Americans in the dataset?
## A:

## Q: compare the results to the f4 results. which source population contributed higher ancestry fraction? which one shares more drift?
## A:
