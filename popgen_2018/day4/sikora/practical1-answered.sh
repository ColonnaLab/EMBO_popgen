
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
## A: Yes, we have a significantly negative f3 value (f3 =-0.017457, Z = -43)


## --------------------------------------------------------------------------------
## admixture f3 statistics - ancient source populations

## Admixture in two East African target populations (Somali from Somalia, Dinka from Sudan); Source populations are Mota (~4,000 year old individual from Ethiopia) and different modern and ancient (LBK) Eurasian populations

qp3Pop -p eafr.f3.par | tee eafr.f3.log

## Q: which combinations of source / target populations show evidence for admixture?
## A: Only the Somali show evidence of admixture; for all Source population combinations

## Q: based on the values of the f3 statistics, which source population is the best proxy for the true admixing population?
## A: The combination with the lowest f3 value is Mota/LBK, incolving early neolithic European farmers (LBK), suggesting that admixture involved groups associated with farming-related migrations during the neolithic 

## Admixture in modern Europeans from three possible sources: Loschbour (Mesolithic hunter-gatherer); LBK (Early Neolithic farmers); Yamnaya (Steppe pastoralists)
    
qp3Pop -p eur.f3.par | tee eur.f3.log

## Q: which populations show neolithic farmer / hunter-gatherer admixture?
## A: All four target populations show evidence for farmer/hunter-gatherer admixture, weakest in Sardinians (Z = -2.8)

## Q: which populations have evidence for Yamnaya admixture
## A: Only French and Spanish show evidence for admixture with Yamnaya (using LBK as the other source population)

## outgroup f3
qp3Pop -p outgroup.f3.par | tee outgroup.f3.log 

## Q: which population shares the most drift with LBK?
grep result: outgroup.f3.log | sort -rnk5,5
## A: the Iceman, an ancient chalcolithic individual from Northern Italy


## --------------------------------------------------------------------------------
## f4/D statistics

## examine par file
cat master.f4.par

## European ancient admixture, here we are testing whether a modern European test population (p3) forms a clade with LBK to the exclusion of Yamnaya 
qpDstat -p eur.f4.par | tee eur.f4.log

## Q: which configuration is consistent with a simple tree?
## A: Only Sardinians show D=0 with LBK, consistent with the two forming a clade in a tree

## Q: how can we interpret the configurations that fail the test?
## A: All other modern test populations share more alleles with Yamnaya than LBK does, indicating admixture from the Steppe

## AA admixture
qpDstat -p aa.f4.par | tee aa.f4.log

## Q: Do we have evidence for admixture?
## A: Yes, each of the three configurations is rejected, and therefore not consistent with a simple tree-like relationship

## Q: Which of the two source population (French, Yoruba) shares more drift with African Americans?
## A: The French, D = -0.006 / Z = 45.49 in the configuration D(Mbuti,AA;French,Yoruba)

## Q: Can we conclude that African Americans have higher proportions of European admixture
## A: No, we can't, as the expected value of f4/D depends both on admixture proportions and the amount of genetic drift shared with each of the source populations


## --------------------------------------------------------------------------------
## f4 ratio estimation

## examine par file
cat master.f4ratio.par

## AA admixture
qpF4ratio -p aa.f4ratio.par | tee aa.f4ratio.log

## Q: what are the admixture proportions for the African Americans in the dataset?
## A: ~16% for all European sources tested

## Q: compare the results to the f4 results. which source population contributed higher ancestry fraction? which one shares more drift?
## A: Yoruba contributed ~84% of ancestry, yet French share more drift with AA
