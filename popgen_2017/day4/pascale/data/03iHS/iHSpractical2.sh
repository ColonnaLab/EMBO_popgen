setwd("~/Documents/workshop/04vcflib/EMBOrun04")

# EDAR: rs3827760
# location: 109,513,601
# 109513601-1000000 = 108513601
# 109513601+1000000 = 110513601
POS_EDAR = 109513601
SNP_EDAR = 'rs3827760'

# LCT (−13910*T): rs4988235
# location: 136,608,646
#136608646 - 1000000 = 135608646
#136608646 + 1000000 = 137608646
POS_LCT = 136608646
SNP_LCT = 'rs4988235'

EDAR_CHB_ihs_in = read.csv(file=paste("CHB.EDAR.ihs.ihs.out"), sep="\t", header=FALSE)
#For iHS the file will be named <outfile>.ihs[.alt].out and formatted as
#<locusID > <physicalPos > <’1’ freq > <ihh1 > <ihh0 > < unstandardized iHS >
colnames(EDAR_CHB_ihs_in) = c('locusID', 'physicalPos', '1freq', 'ihh1', 'ihh0', 'unstandardized_iHS')
EDAR_GBR_ihs_in = read.csv(file=paste("GBR.EDAR.ihs.ihs.out"), sep="\t", header=FALSE)
colnames(EDAR_GBR_ihs_in) = c('locusID', 'physicalPos', '1freq', 'ihh1', 'ihh0', 'unstandardized_iHS')

plot(x=EDAR_CHB_ihs_in[,2], y=EDAR_CHB_ihs_in[,5], type='l') 	# pch=20, cex=0.5, 
lines(x=EDAR_GBR_ihs_in[,2], y=EDAR_GBR_ihs_in[,5], type='l', col='grey') 	# pch=20, cex=0.5, 
points(x=EDAR_CHB_ihs_in[which(EDAR_CHB_ihs_in[,2]==POS_EDAR),2], y=EDAR_CHB_ihs_in[which(EDAR_CHB_ihs_in[,2]==POS_EDAR),5], cex=0.5, col='blue')

FST_EDAR_CHB_GBR=read.csv("EDAR.FST.GBR.CHB.weir.fst", sep='\t', header=TRUE)
plot(x=FST_EDAR_CHB_GBR[,2], y=FST_EDAR_CHB_GBR[,3], pch=20, cex=0.2)
points(x=FST_EDAR_CHB_GBR[which(FST_EDAR_CHB_GBR[,2]==POS_EDAR),2], y=FST_EDAR_CHB_GBR[which(FST_EDAR_CHB_GBR[,2]==POS_EDAR),3], col='blue')


EDAR_CHB_ihs_norm = read.csv("~/Documents/workshop/PascalePractical/03iHS/CHB.EDAR.ihs.norm.out", header=TRUE, sep='\t')
normEDAR_CHB_ihs = matrix(0, nrow=length(EDAR_CHB_ihs_in[,1]), ncol=1)
head(EDAR_CHB_ihs_norm)

f=1
while (f<=length(EDAR_CHB_ihs_norm[,1])){
	if (f==1){
		normEDAR_CHB_ihs[which((EDAR_CHB_ihs_in[,3]>0) & (EDAR_CHB_ihs_in[,3]<=EDAR_CHB_ihs_norm$bin[f])),1] = (EDAR_CHB_ihs_in[which((EDAR_CHB_ihs_in[,3]>0) & (EDAR_CHB_ihs_in[,3]<=EDAR_CHB_ihs_norm$bin[f])),6]*EDAR_CHB_ihs_norm$mean[which(EDAR_CHB_ihs_norm$bin==EDAR_CHB_ihs_norm$bin[f])])/sqrt(EDAR_CHB_ihs_norm$variance[which(EDAR_CHB_ihs_norm$bin==EDAR_CHB_ihs_norm$bin[f])])
	} else {
		normEDAR_CHB_ihs[which((EDAR_CHB_ihs_in[,3]>EDAR_CHB_ihs_norm$bin[f-1]) & (EDAR_CHB_ihs_in[,3]<=EDAR_CHB_ihs_norm$bin[f])),1] = (EDAR_CHB_ihs_in[which((EDAR_CHB_ihs_in[,3]>EDAR_CHB_ihs_norm$bin[f-1]) & (EDAR_CHB_ihs_in[,3]<=EDAR_CHB_ihs_norm$bin[f])),6]*EDAR_CHB_ihs_norm$mean[which(EDAR_CHB_ihs_norm$bin==EDAR_CHB_ihs_norm$bin[f])])/sqrt(EDAR_CHB_ihs_norm$variance[which(EDAR_CHB_ihs_norm$bin==EDAR_CHB_ihs_norm$bin[f])])
	}
	f=f+1
}

EDAR_GBR_ihs_norm = read.csv("~/Documents/workshop/PascalePractical/03iHS/GBR.EDAR.ihs.norm.out", header=TRUE, sep='\t')
normEDAR_GBR_ihs = matrix(0, nrow=length(EDAR_GBR_ihs_in[,1]), ncol=1)
head(EDAR_GBR_ihs_norm)

f=1
while (f<=length(EDAR_GBR_ihs_norm[,1])){
	if (f==1){
		normEDAR_GBR_ihs[which((EDAR_GBR_ihs_in[,3]>0) & (EDAR_GBR_ihs_in[,3]<=EDAR_GBR_ihs_norm$bin[f])),1] = (EDAR_GBR_ihs_in[which((EDAR_GBR_ihs_in[,3]>0) & (EDAR_GBR_ihs_in[,3]<=EDAR_GBR_ihs_norm$bin[f])),6]*EDAR_GBR_ihs_norm$mean[which(EDAR_GBR_ihs_norm$bin==EDAR_GBR_ihs_norm$bin[f])])/sqrt(EDAR_GBR_ihs_norm$variance[which(EDAR_GBR_ihs_norm$bin==EDAR_GBR_ihs_norm$bin[f])])
	} else {
		normEDAR_GBR_ihs[which((EDAR_GBR_ihs_in[,3]>EDAR_GBR_ihs_norm$bin[f-1]) & (EDAR_GBR_ihs_in[,3]<=EDAR_GBR_ihs_norm$bin[f])),1] = (EDAR_GBR_ihs_in[which((EDAR_GBR_ihs_in[,3]>EDAR_GBR_ihs_norm$bin[f-1]) & (EDAR_GBR_ihs_in[,3]<=EDAR_GBR_ihs_norm$bin[f])),6]*EDAR_GBR_ihs_norm$mean[which(EDAR_GBR_ihs_norm$bin==EDAR_GBR_ihs_norm$bin[f])])/sqrt(EDAR_GBR_ihs_norm$variance[which(EDAR_GBR_ihs_norm$bin==EDAR_GBR_ihs_norm$bin[f])])
	}
	f=f+1
}

plot(x=EDAR_CHB_ihs_in[,2], y=normEDAR_CHB_ihs[,1], pch=20, cex=0.5) 
points(x=EDAR_GBR_ihs_in[,2], y=normEDAR_GBR_ihs[,1], pch=20, cex=0.5, col='grey') 	#  
points(x=EDAR_CHB_ihs_in[which(EDAR_CHB_ihs_in[,2]==POS_EDAR),2], y=normEDAR_CHB_ihs[which(EDAR_CHB_ihs_in[,2]==POS_EDAR),1], cex=0.5, col='blue')





LCT_GBR_ihs_in = read.csv(file=paste("GBR.LCT.ihs.ihs.out"), sep="\t", header=FALSE)
colnames(LCT_GBR_ihs_in) = c('locusID', 'physicalPos', '1freq', 'ihh1', 'ihh0', 'unstandardized_iHS')
LCT_CHB_ihs_in = read.csv(file=paste("CHB.LCT.ihs.ihs.out"), sep="\t", header=FALSE)
colnames(LCT_GBR_ihs_in) = c('locusID', 'physicalPos', '1freq', 'ihh1', 'ihh0', 'unstandardized_iHS')

plot(x=LCT_GBR_ihs_in[,2], y=LCT_GBR_ihs_in[,4], type='l') 	# pch=20, cex=0.5, 
lines(x=LCT_CHB_ihs_in[,2], y=LCT_CHB_ihs_in[,4], type='l', col='grey') 	# pch=20, cex=0.5, 
points(x=LCT_GBR_ihs_in[which(LCT_GBR_ihs_in[,2]==POS_LCT),2], y=LCT_GBR_ihs_in[which(LCT_GBR_ihs_in[,2]==POS_LCT),4], cex=0.5, col='blue')

FST_LCT_CHB_GBR=read.csv("LCT.FST.GBR.CHB.weir.fst", sep='\t', header=TRUE)
plot(x=FST_LCT_CHB_GBR[,2], y=FST_LCT_CHB_GBR[,3], pch=20, cex=0.2)
points(x=FST_LCT_CHB_GBR[which(FST_LCT_CHB_GBR[,2]==POS_LCT),2], y=FST_LCT_CHB_GBR[which(FST_LCT_CHB_GBR[,2]==POS_LCT),3], col='blue')


LCT_GBR_ihs_norm = read.csv("~/Documents/workshop/PascalePractical/03iHS/GBR.LCT.ihs.norm.out", header=TRUE, sep='\t')
normLCT_GBR_ihs = matrix(0, nrow=length(LCT_GBR_ihs_in[,1]), ncol=1)
head(LCT_GBR_ihs_norm)

f=1
while (f<=length(LCT_GBR_ihs_norm[,1])){
	if (f==1){
		normLCT_GBR_ihs[which((LCT_GBR_ihs_in[,3]>0) & (LCT_GBR_ihs_in[,3]<=LCT_GBR_ihs_norm$bin[f])),1] = (LCT_GBR_ihs_in[which((LCT_GBR_ihs_in[,3]>0) & (LCT_GBR_ihs_in[,3]<=LCT_GBR_ihs_norm$bin[f])),6]*LCT_GBR_ihs_norm$mean[which(LCT_GBR_ihs_norm$bin==LCT_GBR_ihs_norm$bin[f])])/sqrt(LCT_GBR_ihs_norm$variance[which(LCT_GBR_ihs_norm$bin==LCT_GBR_ihs_norm$bin[f])])
	} else {
		normLCT_GBR_ihs[which((LCT_GBR_ihs_in[,3]>LCT_GBR_ihs_norm$bin[f-1]) & (LCT_GBR_ihs_in[,3]<=LCT_GBR_ihs_norm$bin[f])),1] = (LCT_GBR_ihs_in[which((LCT_GBR_ihs_in[,3]>LCT_GBR_ihs_norm$bin[f-1]) & (LCT_GBR_ihs_in[,3]<=LCT_GBR_ihs_norm$bin[f])),6]*LCT_GBR_ihs_norm$mean[which(LCT_GBR_ihs_norm$bin==LCT_GBR_ihs_norm$bin[f])])/sqrt(LCT_GBR_ihs_norm$variance[which(LCT_GBR_ihs_norm$bin==LCT_GBR_ihs_norm$bin[f])])
	}
	f=f+1
}

normLCT_CHB_ihs = matrix(0, nrow=length(LCT_CHB_ihs_in[,1]), ncol=1)
LCT_CHB_ihs_norm = read.csv("~/Documents/workshop/PascalePractical/03iHS/CHB.LCT.ihs.norm.out", header=TRUE, sep='\t')
head(LCT_CHB_ihs_norm)

f=1
while (f<=length(LCT_CHB_ihs_norm[,1])){
	if (f==1){
		normLCT_CHB_ihs[which((LCT_CHB_ihs_in[,3]>0) & (LCT_CHB_ihs_in[,3]<=LCT_CHB_ihs_norm$bin[f])),1] = (LCT_CHB_ihs_in[which((LCT_CHB_ihs_in[,3]>0) & (LCT_CHB_ihs_in[,3]<=LCT_CHB_ihs_norm$bin[f])),6]*LCT_CHB_ihs_norm$mean[which(LCT_CHB_ihs_norm$bin==LCT_CHB_ihs_norm$bin[f])])/sqrt(LCT_CHB_ihs_norm$variance[which(LCT_CHB_ihs_norm$bin==LCT_CHB_ihs_norm$bin[f])])
	} else {
		normLCT_CHB_ihs[which((LCT_CHB_ihs_in[,3]>LCT_CHB_ihs_norm$bin[f-1]) & (LCT_CHB_ihs_in[,3]<=LCT_CHB_ihs_norm$bin[f])),1] = (LCT_CHB_ihs_in[which((LCT_CHB_ihs_in[,3]>LCT_CHB_ihs_norm$bin[f-1]) & (LCT_CHB_ihs_in[,3]<=LCT_CHB_ihs_norm$bin[f])),6]*LCT_CHB_ihs_norm$mean[which(LCT_CHB_ihs_norm$bin==LCT_CHB_ihs_norm$bin[f])])/sqrt(LCT_CHB_ihs_norm$variance[which(LCT_CHB_ihs_norm$bin==LCT_CHB_ihs_norm$bin[f])])
	}
	f=f+1
}

plot(x=LCT_GBR_ihs_in[,2], y=normLCT_GBR_ihs[,1], type='l') 	# pch=20, cex=0.5, 
lines(x=LCT_CHB_ihs_in[,2], y=normLCT_CHB_ihs[,1], type='l', col='grey') 	# pch=20, cex=0.5, 
points(x=LCT_GBR_ihs_in[which(LCT_GBR_ihs_in[,2]==POS_LCT),2], y=normLCT_GBR_ihs[which(LCT_GBR_ihs_in[,2]==POS_LCT),1], cex=0.5, col='blue')

