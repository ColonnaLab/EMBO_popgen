setwd('~/Documents/workshop/01vcftools/EMBOtest_run01')

# read vcftools files in
# FST AFR/EUR
FST_AFR_EUR_in = read.csv("EMBOrun01_FST_AfrEur.weir.fst", sep='\t')

POS_LCT = 136608646
POS_EDAR = 109513601

FST_AfrEur_data = FST_AFR_EUR_in[-which(is.na(FST_AFR_EUR_in[,3])),]
FST_AfrEur_data[which(FST_AfrEur_data[,3]<0),3] = 0

FST_AfrEur_distr = sort(FST_AfrEur_data[,3])
FST_AfrEur_distrQT = quantile(FST_AfrEur_distr, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
#        1%         5%        10%        25%        50%        75%        90%        95%        99% 
#0.00000000 0.00000000 0.00000000 0.00000000 0.00112564 0.01383270 0.06640710 0.13275700 0.32299400 

FSTdata = FST_AfrEur_data
SNP_POS = POS_LCT
nb_SNP_randomDISTR = 28440
SNPregion = 20000

# IDENTIFY THE REGION OF INTEREST ±SNPregion
# i. the start and end position in BP around the
# SNP of interest
SNPfrom_BP = SNP_POS - SNPregion
SNPto_BP = SNP_POS + SNPregion
# and the corresponding index (row number) in the
# data
SNPfrom_id = max(which(FSTdata[,2]<=SNPfrom_BP))
SNPto_id = min(which(FSTdata[,2]>=SNPto_BP))
length(FSTdata[SNPfrom_id:SNPto_id, 2])
# ii. select these rows from the data table
FSTdata_SNP = FSTdata[SNPfrom_id:SNPto_id, ]

plot(ylim=c(0,1), x=FSTdata_SNP[,2], y=FSTdata_SNP[,3], xlab='pos', ylab='FST', pch=20, cex=0.2)
points(x=FSTdata_SNP[which(FSTdata_SNP[,2]==SNP_POS),2],  y=FSTdata_SNP[which(FSTdata_SNP[,2]==SNP_POS),3], col='blue')
abline(h=FST_AfrEur_distrQT[[9]], lty=2)

#

# Randomly sample FSToutside the SNP region
fst_random = sample(FSTdata[-(SNPfrom_id:SNPto_id), 3], nb_SNP_randomDISTR, replace=FALSE)
wilcox.test(fst_random, FSTdata_SNP[, 3], paired = FALSE, alternative = "greater")
#	Wilcoxon rank sum test with continuity correction
#
#data:  fst_random and FSTdata_SNP[, 3]
#W = 10618000, p-value = 0.003298
#alternative hypothesis: true location shift is greater than 0
for (s in 1:100){
	fst_random = sample(FSTdata[-(SNPfrom_id:SNPto_id), 3], nb_SNP_randomDISTR, replace=FALSE)
	Pvalue = wilcox.test(fst_random, FSTdata_SNP[, 3], paired = FALSE, alternative = "greater")$p.value
	if (Pvalue<=0.05){
		cat("pvalue=", wilcox.test(fst_random, FSTdata_SNP[, 3], paired = FALSE, alternative = "greater")$p.value, sep='', '\n')
	}	
}

# Try with EDAR:
SNP_POS = POS_EDAR

# all pvalues > 0.05

# Try with EDAR and AFR_EAS
# FST AFR/EAS
FST_AFR_EAS_in = read.csv("EMBOrun01_FST_AfrEas.weir.fst", sep='\t')

FST_AfrEas_data = FST_AFR_EAS_in[-which(is.na(FST_AFR_EAS_in[,3])),]
FST_AfrEas_data[which(FST_AfrEas_data[,3]<0),3] = 0

FST_AfrEas_distr = sort(FST_AfrEas_data[,3])
FST_AfrEas_distrQT = quantile(FST_AfrEas_distr, c(0.01, 0.05, 0.1, .25, .50,  .75, .90, 0.95, .99))
#         1%         5%        10%        25%        50%        75%        90%        95%        99% 
#0.00000000 0.00000000 0.00000000 0.00000000 0.00112834 0.01248580 0.06984910 0.14816420 0.38567400 

FSTdata = FST_AfrEas_data
SNP_POS = POS_EDAR
nb_SNP_randomDISTR = 28440
SNPregion = 20000

# IDENTIFY THE REGION OF INTEREST ±SNPregion
# i. the start and end position in BP around the
# SNP of interest
SNPfrom_BP = SNP_POS - SNPregion
SNPto_BP = SNP_POS + SNPregion
# and the corresponding index (row number) in the
# data
SNPfrom_id = max(which(FSTdata[,2]<=SNPfrom_BP))
SNPto_id = min(which(FSTdata[,2]>=SNPto_BP))
length(FSTdata[SNPfrom_id:SNPto_id, 2])
# ii. select these rows from the data table
FSTdata_SNP = FSTdata[SNPfrom_id:SNPto_id, ]

plot(ylim=c(0,1), x=FSTdata_SNP[,2], y=FSTdata_SNP[,3], xlab='pos', ylab='FST', pch=20, cex=0.2)
points(x=FSTdata_SNP[which(FSTdata_SNP[,2]==SNP_POS),2],  y=FSTdata_SNP[which(FSTdata_SNP[,2]==SNP_POS),3], col='blue')
abline(h=FST_AfrEas_distrQT[[9]], lty=2)

for (s in 1:100){
	fst_random = sample(FSTdata[-(SNPfrom_id:SNPto_id), 3], nb_SNP_randomDISTR, replace=FALSE)
	Pvalue = wilcox.test(fst_random, FSTdata_SNP[, 3], paired = FALSE, alternative = "greater")$p.value
	if (Pvalue<=0.05){
		cat("pvalue=", wilcox.test(fst_random, FSTdata_SNP[, 3], paired = FALSE, alternative = "greater")$p.value, sep='', '\n')
	}
}
# not significant for EDAR
