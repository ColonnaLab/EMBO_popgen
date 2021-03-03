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

