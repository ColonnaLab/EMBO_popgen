#Number of features
zcat /home/delaneau/chr20.RAW.vcf.gz | grep -v "#" | wc -l
zcat /home/delaneau/chr20.RAW.vcf.gz | grep "CHROM" | cut -f10- | wc -w

#Compute missing data rates
vcftools --gzvcf /home/delaneau/chr20.RAW.vcf.gz --missing-site --stdout > chr20.RAW.missing.txt

#Plot missing data rates
R
DATA = read.table("chr20.RAW.missing.txt", head=TRUE)
pdf("chr20.RAW.missing.pdf", 12, 4)
plot(DATA$F_MISS, xlab="Variant index", ylab="Missing data rate", col=ifelse(DATA$F_MISS >= 0.1, "red", "grey"))
abline(h=0.1, col="red")
dev.off()
quit(save="no")

#Remove poorly called variants
vcftools --gzvcf /home/delaneau/chr20.RAW.vcf.gz --max-missing 0.9 --recode --stdout | bgzip -c > chr20.STEP1.vcf.gz
tabix -p vcf chr20.STEP1.vcf.gz

#Compute variant frequencies for European samples
vcftools --gzvcf /home/delaneau/reference/chr20.EUR.vcf.gz --freq2 --stdout | sed '1d' | awk '{print $2, $4, $5, $6}' > chr20.EUR.freq
vcftools --gzvcf chr20.STEP1.vcf.gz --freq2 --stdout | sed '1d' | awk '{print $2, $4, $5, $6}' > chr20.STEP1.freq

#Plot the two sets of frequencies
R
EXP = read.table("chr20.EUR.freq", head=FALSE)
colnames(EXP) = c("pos", "tot_exp", "ref_exp", "alt_exp")
OBS = read.table("chr20.STEP1.freq", head=FALSE)
colnames(OBS) = c("pos", "tot_obs", "ref_obs", "alt_obs")
M = merge(EXP, OBS, by="pos")
M$pvalue = apply(M, 1, FUN=function(x) fisher.test(matrix(round(c(x[2]*x[3], x[2]*x[4], x[5]*x[6], x[5]*x[7])), ncol=2))$p.value)
pdf("chr20.STEP1.frequencies.pdf")
plot(M$alt_exp, M$alt_obs, xlab="ALT frequency in Reference", ylab="ALT frequency in Observed", col=ifelse(M$pvalue < 1e-10, "red", "black"))
legend("bottomright", legend=c("pvalue > 1e-10", "pvalue < 1e-10"), fill=c("black","red"), bg="white")
dev.off()
write.table(cbind(rep(20, sum(M$pvalue < 1e-10)), M$pos[M$pvalue < 1e-10]), "chr20.STEP1.filtered.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
quit(save="no")

#Remove frequency discordant variants
vcftools --gzvcf chr20.STEP1.vcf.gz --exclude-positions chr20.STEP1.filtered.txt --recode --stdout | bgzip -c > chr20.STEP2.vcf.gz

#Compute HWE
vcftools --gzvcf chr20.STEP2.vcf.gz --hardy --stdout | sed '1d' | awk '{ print $2, $3, $4, $6 }' > chr20.STEP2.hwe.txt

#Plot per variant HWE test
R
DATA = read.table("chr20.STEP2.hwe.txt", head=FALSE)
pdf("variant_non_hwe.pdf", 12, 4)
plot(-log10(DATA$V4), xlab="Variant index", ylab="-log10(HWE test)", col=ifelse(-log10(DATA$V4) > 5, "red", "grey"))
abline(h=5, col="red")
dev.off()
quit(save="no")

#Remove variant violating HWE
cat chr20.STEP2.hwe.txt | awk '{ if ($4 < 1e-5) print "20", $1 }' > chr20.STEP2.filtered.txt
vcftools --gzvcf chr20.STEP2.vcf.gz --exclude-positions chr20.STEP2.filtered.txt --recode --stdout | bgzip -c > chr20.STEP3.vcf.gz

#Compute the per sample missing rates
vcftools --gzvcf chr20.STEP3.vcf.gz --missing-indv --stdout > chr20.STEP3.missing.txt

#Plot the per sample missing rates
R
DATA=read.table("chr20.STEP3.missing.txt", head=TRUE)
pdf("chr20.STEP3.missing.pdf")
plot(DATA$F_MISS, xlab="Sample Index", ylab="Missing data rate", col=ifelse(DATA$F_MISS > 0.05, "red", "black"), main="Missing data report per individual")
abline(h=0.05, col="red")
dev.off()
write.table( DATA$INDV[DATA$F_MISS > 0.05], "chr20.STEP3.filtered.txt", quote=FALSE, col.names=FALSE, row.names=FALSE)
quit(save="no")

#Remove poorly called samples
vcftools --gzvcf chr20.STEP3.vcf.gz --remove chr20.STEP3.filtered.txt --recode --stdout | bgzip -c > chr20.STEP4.vcf.gz
$TABIX -p vcf chr20.STEP4.vcf.gz

#Find related samples
plink --vcf chr20.STEP4.vcf.gz --genome --ppc-gap 100 --out chr20.STEP4

#Plot relatedness between samples
R
DATA = read.table("chr20.STEP4.genome", head=TRUE)
pdf("chr20.STEP4.genome.pdf", 12, 4)
par(mfrow=c(1,3))
plot(DATA$Z0, DATA$Z1, xlab="P(IBD=0)", ylab="P(IBD=1)", main="IBD0 versus IBD1")
plot(DATA$Z0, DATA$Z2, xlab="P(IBD=0)", ylab="P(IBD=2)", main="IBD0 versus IBD2")
plot(DATA$Z1, DATA$Z2, xlab="P(IBD=1)", ylab="P(IBD=2)", main="IBD1 versus IBD2")
dev.off()
quit(save="no")

#Remove duplicate and related samples
vcftools --gzvcf chr20.STEP4.vcf.gz --remove-indv HG01606dup0 --remove-indv HG01770HG02238 --recode --stdout | bgzip -c > chr20.STEP5.vcf.gz
tabix -p vcf chr20.STEP5.vcf.gz

#Merge our data with 1000 Genomes for PCA
bcftools merge -m id -Oz -o chr20.MERGED.vcf.gz /home/delaneau/reference/chr20.ALL.vcf.gz chr20.STEP5.vcf.gz
tabix -p vcf chr20.MERGED.vcf.gz

#Perform PCA
QTLtools_1.1_Ubuntu16.04_x86_64 pca --vcf chr20.MERGED.vcf.gz --scale --center --distance 50000 --maf 0.05 --out chr20.MERGED

#Plot PCA and extract outliers
R
PCA = read.table("chr20.MERGED.pca", head = TRUE)
PCA = data.frame(V1=colnames(PCA)[2:ncol(PCA)], t(PCA[1,2:ncol(PCA)]), t(PCA[2,2:ncol(PCA)]))
POP = read.table("/home/delaneau/reference/populations.txt", head=FALSE)
DATA = merge(POP, PCA, by="V1")
head(DATA)
pdf("chr20.MERGED.pca.pdf", 10, 5)
par(mfrow=c(1,2))
plot(DATA$X1[DATA$V3 != "SPA"], DATA$X2[DATA$V3 != "SPA"], xlab="PC1", ylab="PC2", pch=20, col="grey")
points(DATA$X1[DATA$V3 == "SPA"], DATA$X2[DATA$V3 == "SPA"], pch=20, col="red")
legend("topleft", legend=c("Our samples", "1000 genomes samples"), fill=c("red","grey"))
plot(DATA$X1[DATA$V3 != "SPA"], DATA$X2[DATA$V3 != "SPA"], xlab="PC1", ylab="PC2", pch=20, col=DATA$V3[DATA$V3 != "SPA"])
legend("topleft", legend=unique(DATA$V3[DATA$V3 != "SPA"]), fill=unique(DATA$V3[DATA$V3 != "SPA"]))
dev.off()
write.table(DATA[DATA$X2 > -5 & DATA$V3 == "SPA", 1], "chr20.STEP5.filtered.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
quit(save="no")

#Remove outliers comming from PCA to get the final data set
vcftools --gzvcf chr20.STEP5.vcf.gz --remove chr20.STEP5.filtered.txt --recode --stdout | bgzip -c > chr20.FINAL.vcf.gz
tabix -p vcf chr20.FINAL.vcf.gz

