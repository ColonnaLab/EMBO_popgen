setwd('~/Documents/workshop/01vcftools/EMBOtest_run01')

# read TajD files
TajD_AFR = read.csv("EMBOrun01_TajD_AFR.Tajima.D", sep='\t')
TajD_EAS = read.csv("EMBOrun01_TajD_EAS.Tajima.D", sep='\t')
TajD_EUR = read.csv("EMBOrun01_TajD_EUR.Tajima.D", sep='\t')

nEUR = 503
nAFR = 661
nEAS = 504

POS_LCT = 136608646
POS_EDAR = 109513601

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# 		LCT						#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
LCTindEUR = max(which(TajD_EUR[,2]<=POS_LCT))
LCTindEAS = max(which(TajD_EAS[,2]<=POS_LCT))
LCTindAFR = max(which(TajD_AFR[,2]<=POS_LCT))
LCTind = LCTindEUR

SNP_FROM = LCTind - 100
SNP_TO = LCTind + 100

plot(x=TajD_EUR[SNP_FROM:SNP_TO,2], y=TajD_EUR[SNP_FROM:SNP_TO,4], type='l', lty=2, xlab="pos", ylab="TajD")
lines(x=TajD_EAS[SNP_FROM:SNP_TO,2], y=TajD_EAS[SNP_FROM:SNP_TO,4], type='l', col='grey', lty=2)
lines(x=TajD_AFR[SNP_FROM:SNP_TO,2], y=TajD_AFR[SNP_FROM:SNP_TO,4], type='l', col='yellow', lty=2)
points(x=TajD_EUR[LCTind,2], y=TajD_EUR[LCTind,4], pch=20, cex=1)
points(x=TajD_EAS[LCTind,2], y=TajD_EAS[LCTind,4], pch=20, cex=1, col='grey')
points(x=TajD_AFR[LCTind,2], y=TajD_AFR[LCTind,4], pch=20, cex=1, col='orange')

TajD_EUR[LCTind,4]
TajD_EAS[LCTind,4]
TajD_AFR[LCTind,4]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# 		EDAR					#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
EDARindEUR = max(which(TajD_EUR[,2]<=POS_EDAR))
EDARindEAS = max(which(TajD_EAS[,2]<=POS_EDAR))
EDARindAFR = max(which(TajD_AFR[,2]<=POS_EDAR))
EDARind = EDARindEAS

SNP_FROM = EDARind - 100
SNP_TO = EDARind + 100

plot(x=TajD_EUR[SNP_FROM:SNP_TO,2], y=TajD_EUR[SNP_FROM:SNP_TO,4], type='l', lty=2, xlab="pos", ylab="TajD")
lines(x=TajD_EAS[SNP_FROM:SNP_TO,2], y=TajD_EAS[SNP_FROM:SNP_TO,4], type='l', col='grey', lty=2)
lines(x=TajD_AFR[SNP_FROM:SNP_TO,2], y=TajD_AFR[SNP_FROM:SNP_TO,4], type='l', col='yellow', lty=2)
points(x=TajD_EUR[EDARind,2], y=TajD_EUR[EDARind,4], pch=20, cex=1)
points(x=TajD_EAS[EDARind,2], y=TajD_EAS[EDARind,4], pch=20, cex=1, col='grey')
points(x=TajD_AFR[EDARind,2], y=TajD_AFR[EDARind,4], pch=20, cex=1, col='orange')

TajD_EUR[EDARind,4]
TajD_EAS[EDARind,4]
TajD_AFR[EDARind,4]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# 		DISTRIBUTION			#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
TajD_DISTR = sort(c(TajD_AFR[,4], TajD_EAS[,4], TajD_EUR[,4]), decreasing = FALSE)
plot(density(TajD_DISTR))
points(x=TajD_EUR[LCTind,4] , y=0, pch=20)
points(x=TajD_EAS[LCTind,4] , y=0)
points(x=TajD_EUR[EDARind,4] , y=0, col='blue')
points(x=TajD_EAS[EDARind,4] , y=0, pch=20, col='blue')

