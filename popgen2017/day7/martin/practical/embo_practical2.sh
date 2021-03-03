
#######################################
# practical exercises qpGraph / alder #
#######################################


## --------------------------------------------------------------------------------
## qpGraph

## examine par files
cat qpgraph.par

## base graph (MA1 / Han)
cat base.qpgraph.graph
qpGraph -p qpgraph.par -g base.qpgraph.graph -d base.qpgraph.dot -o base.qpgraph.graphout | tee base.qpgraph.log
dot -Tpdf base.qpgraph.dot > base.qpgraph.pdf

## add Surui unadmixed to graph (MA1 / Han)
cat n01_surui.qpgraph.graph
cat n02_surui.qpgraph.graph

qpGraph -p qpgraph.par -g n01_surui.qpgraph.graph -d n01_surui.qpgraph.dot -o n01_surui.qpgraph.graphout | tee n01_surui.qpgraph.log
dot -Tpdf n01_surui.qpgraph.dot > n01_surui.qpgraph.pdf
qpGraph -p qpgraph.par -g n02_surui.qpgraph.graph -d n02_surui.qpgraph.dot -o n02_surui.qpgraph.graphout | tee n02_surui.qpgraph.log
dot -Tpdf n02_surui.qpgraph.dot > n02_surui.qpgraph.pdf

## Q: does the graph fit the observed f-statistics?
## Q: which configurations show the worst fit?

## add Surui admixed to graph (MA1 / Han)
cat x01_surui.qpgraph.graph
qpGraph -p qpgraph.par -g x01_surui.qpgraph.graph -d x01_surui.qpgraph.dot -o x01_surui.qpgraph.graphout | tee x01_surui.qpgraph.log
dot -Tpdf x01_surui.qpgraph.dot > x01_surui.qpgraph.pdf

## Q: does the graph fit?
## Q: what are the admixture proportions?

## add Mayan unadmixed to graph
cat n01_mayan.x01_surui.qpgraph.graph
qpGraph -p qpgraph.par -g n01_mayan.x01_surui.qpgraph.graph -d n01_mayan.x01_surui.qpgraph.dot -o n01_mayan.x01_surui.qpgraph.graphout | tee n01_mayan.x01_surui.qpgraph.log
dot -Tpdf n01_mayan.x01_surui.qpgraph.dot > n01_mayan.x01_surui.qpgraph.pdf

## Q: does the graph fit?
## Q: how can we intrepret the outlier f-statistics?


## --------------------------------------------------------------------------------
## to further explore if time permits

## qpGraph for modern Europeans as mixtures of ancient sources (Loschbour / LBK / Yamnaya)
cat sard.qpgraph.graph
qpGraph -p qpgraph.par -g sard.qpgraph.graph -d sard.qpgraph.dot -o sard.qpgraph.graphout | tee sard.qpgraph.log


## --------------------------------------------------------------------------------
## ALDER

## African - American admixture
cat master.alder.par | perl -p -e 's/<ADMIXPOP>/AA/' | perl -p -e 's/<REFPOPS>/Yoruba;French/' | perl -p -e 's/<RAW_OUT>/aa.alder.ld/' > aa.alder.par
alder -p aa.alder.par | tee aa.alder.log

## Q: did alder detect admixture?
## Q: what admxixture time (in generations) did we get?
## Q: do we get consistent decay rates in all curves? 

cat master.alder.par | perl -p -e 's/<ADMIXPOP>/AA/' | perl -p -e 's/<REFPOPS>/Yoruba;Mbuti;Dinka;French;Spanish;Basque/' | perl -p -e 's/<RAW_OUT>/aa.alder.ld/' > aa1.alder.par
alder -p aa1.alder.par | tee aa1.alder.log

## Q: which combination of source populations shwos admixture?
## Q: do we find consistent admixture time?


## mayan / surui admixture
cat master.alder.par | perl -p -e 's/<ADMIXPOP>/Mayan/' | perl -p -e 's/<REFPOPS>/Han;Spanish/' | perl -p -e 's/<RAW_OUT>/mayan.alder.ld/' > mayan.alder.par
alder -p mayan.alder.par | tee mayan.alder.log
cat master.alder.par | perl -p -e 's/<ADMIXPOP>/Surui/' | perl -p -e 's/<REFPOPS>/Han;Spanish/' | perl -p -e 's/<RAW_OUT>/surui.alder.ld/' > surui.alder.par
alder -p surui.alder.par | tee surui.alder.log

## Q: which native americans show admixture?
## Q: what is different between the results for the two target populations?

## eafr admixture
cat master.alder.par | perl -p -e 's/<ADMIXPOP>/Somali/' | perl -p -e 's/<REFPOPS>/Mota;LBK;Sardinian;French;Han/' | perl -p -e 's/<RAW_OUT>/eafr.alder.ld/' > eafr.alder.par
alder -p eafr.alder.par | tee eafr.alder.log

## Q: what admixture time to we get?
## Q: which source population is the best proxy for the actual source?

