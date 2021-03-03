
#######################################
# practical exercises qpGraph / alder #
#######################################


## In this practical we will explore a specific topic of human population history (ancestral populations of modern Europeans) in greater detail using the full suite of more advanced f-statistic methods (qpAdm, qpGraph)

## --------------------------------------------------------------------------------
## Recapitulate f3 and f4 results for modern Europeans using ancient sources: Loschbour (Mesolithic hunter-gatherer); LBK (Early Neolithic farmers); Yamnaya (Steppe pastoralists)

cat eur.f3.log
cat eur.f4.log

## Q: Do any combinations of source / target populations show evidence for admixture?
## A: None of the tested configurations show evidence for admixture 

## Q: Can we conclude that our modern Native Americans are not admixed between the source populations?
## A: No we cannot. A positive f3 can occur even for an admixed population if genetic drift in the target population post-admixture was sufficinetly high 


## --------------------------------------------------------------------------------
## qpWave

## here we first use qpwave to confirm that our three source populations can be distinguished with the outgroup set (i.e. are related through at least 3 stream of ancestry)

## examine par file, outgroups and output 
cat ref.qpadm.par
cat right.txt
qpAdm -p ref.qpwave.par | tee ref.qpwave.log


## --------------------------------------------------------------------------------
## qpAdm

## examine par file and output on spanish as example
cat spanish.qpadm.par
qpAdm -p spanish.qpadm.par | tee spanish.qpadm.log


## basque
qpAdm -p basque.qpadm.par | tee basque.qpadm.log

## French
qpAdm -p french.qpadm.par | tee french.qpadm.log

## Sardinian
qpAdm -p sardinian.qpadm.par | tee sardinian.qpadm.log


## Q: Which of the modern Europeans shows Yamnaya admixture? What fraction is estimated to originate from Yamnaya?
## A: Spanish (30%, pattern 100); Basque (39%, pattern 100); French (45%, pattern 100)

## Q: Which population requires ancestry from hunter-gatherers (Loschbour) to fit?
## A: None of the test populations requires additional HG ancestry (all p > 0.05 for pattern 100); However for Basque the fit improves marginally when adding Loschbour (11.% ancestry, nested p-value 0.06)

## Q: Do Sardinians show evidence for admixture?
## A: No, they are consistent with being a clade with LBK (early neolithic farmers), pattern 101 not rejected (p = 0.13)



## --------------------------------------------------------------------------------
## qpGraph

## examine par files
cat qpgraph.par

## base graph 
cat base.qpgraph.graph
qpGraph -p qpgraph.par -g base.qpgraph.graph -d base.qpgraph.dot -o base.qpgraph.graphout | tee base.qpgraph.log
dot -Tpdf base.qpgraph.dot > base.qpgraph.pdf

## Q: does the graph fit the observed f-statistics?
## A: Yes it does

## add Sardinian unadmixed to graph
cat n01.sar.qpgraph.graph
qpGraph -p qpgraph.par -g n01.sar.qpgraph.graph -d n01.sar.qpgraph.dot -o n01.sar.qpgraph.graphout | tee n01.sar.qpgraph.log
dot -Tpdf n01.sar.qpgraph.dot > n01.sar.qpgraph.pdf

## Q: does the graph fit?
## A: No it does not, worst Z = -5.2

## add admixture from Yamnaya 
qpGraph -p qpgraph.par -g x01.sar.qpgraph.graph -d x01.sar.qpgraph.dot -o x01.sar.qpgraph.graphout | tee x01.sar.qpgraph.log
dot -Tpdf x01.sar.qpgraph.dot > x01.sar.qpgraph.pdf

## Q: has the fit improved?
## A: Only marginally, worst Z-score is now -4.2

## Q: how much Yamnaya admixture does the fit show
## A: only 1%, indicating that this model is not promising for further exploration

