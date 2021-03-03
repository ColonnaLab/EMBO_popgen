
#######################################
# practical exercises qpGraph / alder #
#######################################


## In this practical we will explore a specific topic of human population history (ancestral populations of modern Europeans) in greater detail using the full suite of more advanced f-statistic methods (qpAdm, qpGraph)

## --------------------------------------------------------------------------------
## Recapitulate f3 and f4 results for modern Europeans using ancient sources: Loschbour (Mesolithic hunter-gatherer); LBK (Early Neolithic farmers); Yamnaya (Steppe pastoralists)

cat eur.f3.log
cat eur.f4.log

## Q: Do any combinations of source / target populations show evidence for admixture?
## A:

## Q: Can we conclude that our modern Native Americans are not admixed between the source populations?
## A:

## --------------------------------------------------------------------------------
## qpWave

## here we first use qpwave to confirm that our three source populations can be distinguished with the outgroup set (i.e. are related through at least 3 stream of ancestry)

## examine par file, outgroups and output 
cat ref.qpadm.par
cat right.txt
qpAdm -p ref.qpwave.par | tee ref.qpwave.log


## --------------------------------------------------------------------------------
## qpAdm

## examine par file, outgroups and output on spanish as example
cat spanish.qpadm.par
cat right.txt
qpAdm -p spanish.qpadm.par | tee spanish.qpadm.log


## basque
qpAdm -p basque.qpadm.par | tee basque.qpadm.log

## French
qpAdm -p french.qpadm.par | tee french.qpadm.log

## Sardinian
qpAdm -p sardinian.qpadm.par | tee sardinian.qpadm.log


## Q: Which of the modern Europeans shows Yamnaya admixture? What fraction is estimated to originate from Yamnaya?
## A:

## Q: Which population requires ancestry from hunter-gatherers (Loschbour) to fit?
## A:

## Q: Do Sardinians show evidence for admixture?
## A:



## --------------------------------------------------------------------------------
## qpGraph

## examine par files
cat qpgraph.par

## base graph 
cat base.qpgraph.graph
qpGraph -p qpgraph.par -g base.qpgraph.graph -d base.qpgraph.dot -o base.qpgraph.graphout | tee base.qpgraph.log
dot -Tpdf base.qpgraph.dot > base.qpgraph.pdf

## Q: does the graph fit the observed f-statistics?
## A:

## add Sardinian unadmixed to graph
cat n01.sar.qpgraph.graph
qpGraph -p qpgraph.par -g n01.sar.qpgraph.graph -d n01.sar.qpgraph.dot -o n01.sar.qpgraph.graphout | tee n01.sar.qpgraph.log
dot -Tpdf n01.sar.qpgraph.dot > n01.sar.qpgraph.pdf

## Q: does the graph fit?
## A:

## add admixture from Yamnaya 
qpGraph -p qpgraph.par -g x01.sar.qpgraph.graph -d x01.sar.qpgraph.dot -o x01.sar.qpgraph.graphout | tee x01.sar.qpgraph.log
dot -Tpdf x01.sar.qpgraph.dot > x01.sar.qpgraph.pdf

## Q: has the fit improved?
## A:

## Q: how much Yamnaya admixture does the fit show
## A:

