
##########################################
# practical exercises f3 / f4 / f4 ratio #
##########################################

## --------------------------------------------------------------------------------
## examine the dataset

## number of individuals and SNPs
wc -l ho_anc.indiv
wc -l ho_anc.snp

## populations and sample sizes
awk '{print $3}' ho_anc.indiv | sort | uniq -c
    

## --------------------------------------------------------------------------------
## admixture f3 statistics

## examine par file
cat master.f3.par


## modern source populations

## African - American admixture
cat master.f3.par | perl -p -e 's/<POPFILE>/aa.f3.popfile/' > aa.f3.par
qp3Pop -p aa.f3.par | tee aa.f3.log

## Q: Do we have evidence for admixture?


## Native American admixture
cat master.f3.par | perl -p -e 's/<POPFILE>/natam.f3.popfile/' > natam.f3.par
qp3Pop -p natam.f3.par | tee natam.f3.log

## Q: which combinations of source / target populations show evidence for admixture?
## Q: are results qualitatively different between the Mayan and Surui target populations? if so why?


## ancient source populations

## East African admixture
cat master.f3.par | perl -p -e 's/<POPFILE>/eafr.f3.popfile/' > eafr.f3.par
qp3Pop -p eafr.f3.par | tee eafr.f3.log

## Q: which combinations of source / target populations show evidence for admixture?
## Q: can we say anything about the likely source population for admixture?


## European ancient admixture
cat master.f3.par | perl -p -e 's/<POPFILE>/eur.f3.popfile/' > eur.f3.par
qp3Pop -p eur.f3.par | tee eur.f3.log

## Q: which populations show neolithic farmer / hunter-gatherer admixture?
## Q: which populations have Yamnaya ancestry?


## outgroup f3
cat master.f3.par | perl -p -e 's/<POPFILE>/outgroup.f3.popfile/' > outgroup.f3.par
qp3Pop -p outgroup.f3.par | tee outgroup.f3.log 

## Q: which population shares the most drift with LBK?
grep result: outgroup.f3.log | sort -rnk5,5


## --------------------------------------------------------------------------------
## f4/D statistics

## examine par file
cat master.f4.par

## AA admixture
cat master.f4.par | perl -p -e 's/<POPFILE>/aa.f4.popfile/' > aa.f4.par
qpDstat -p aa.f4.par | tee aa.f4.log

## Q: do we have evidence for admixture?
## Q: which source population shares more drift with AA?
## Q: how can we interpret this result?


## NatAm admixture
cat master.f4.par | perl -p -e 's/<POPFILE>/natam.f4.popfile/' > natam.f4.par
qpDstat -p natam.f4.par | tee natam.f4.log

## Q: do we have evidence for admixture?
## Q: which population p2 shares more drift with native americans?
## Q: how can we interpret these results?


## European ancient admixture
cat master.f4.par | perl -p -e 's/<POPFILE>/eur.f4.popfile/' > eur.f4.par
qpDstat -p eur.f4.par | tee eur.f4.log

## Q: which configuration is consistent with an unrooted tree?
## Q: how can we interpret the configurations that fail the test?


## --------------------------------------------------------------------------------
## f4 ratio statistics

## examine par file
cat master.f4ratio.par

## AA admixture
cat master.f4ratio.par | perl -p -e 's/<POPFILE>/aa.f4ratio.popfile/' > aa.f4ratio.par
qpF4ratio -p aa.f4ratio.par | tee aa.f4ratio.log

## Q: what are the admixture proportions for the African Americans in the dataset?
## Q: compare the results to the f4 results. which source population contributed higher ancestry fraction? which one shares more drift?


## NatAm admixture 
cat master.f4ratio.par | perl -p -e 's/<POPFILE>/natam.f4ratio.popfile/' > natam.f4ratio.par
qpF4ratio -p natam.f4ratio.par | tee natam.f4ratio.log

## Q: what are the admixture proportions for the Native American test populations?


## --------------------------------------------------------------------------------
## to further explore if time permits

## Neandertal admixture into Eurasians
## African admixture into Southern Europe
