#Get the input data
cp /home/delaneau/chr20.FINAL.vcf.gz .
tabix –fp vcf chr20.FINAL.vcf.gz

#Extract data for the example chunk
bcftools view -r 20:800000-3200000 -Oz -o chr20.chunk1.vcf.gz chr20.FINAL.vcf.gz
tabix -p vcf chr20.chunk1.vcf.gz

#Convert VCF file in GEN/SAMPLE format
bcftools convert -g chr20.chunk1 chr20.chunk1.vcf.gz

#Phase the genotype data
shapeit --input-gen chr20.chunk1.gen.gz chr20.chunk1.samples --input-map /home/delaneau/reference/chr20.genetic_map.txt.gz --output-max chr20.chunk1.phased

#Compute haplotype frequencies
cat chr20.chunk1.phased.haps | grep '876578_T_C\|877205_G_T\|880964_A_C\|882429_A_G' | cut -d' ' -f6- | /home/delaneau/transpose | sort | uniq -c

#APPROACH1: Pre-phasing
impute2 -use_prephased_g -phase -known_haps_g chr20.chunk1.phased.haps -m /home/delaneau/reference/chr20.genetic_map.txt.gz -h /home/delaneau/reference/chr20.EUR.hap.gz -l /home/delaneau/reference/chr20.EUR.legend.gz -int 1e6 3e6 -buffer 200 -o chr20.chunk1.imputed.approach1

#APPROACH1: plot performance
R
DATA1 = read.table("chr20.chunk1.imputed.approach1_info", head=TRUE)
DATA1 = DATA1[DATA1$snp_id == "---", ]
DATA1$maf = ifelse(DATA1$exp_freq_a1 < 0.5, DATA1$exp_freq_a1, 1- DATA1$exp_freq_a1)
DATA1$bin = cut(DATA1$maf, breaks = seq(0,0.5,0.05), labels = 1:10, include.lowest = TRUE)
pdf("chr20.chunk1.imputed.approach1.pdf")
plot(by(DATA1$maf, DATA1$bin, mean), by(DATA1$info, DATA1$bin, mean), type="l", xlab="MAF", ylab="Mean info score", xlim=c(0,0.5), ylim=c(0.6, 1.0))
dev.off()
quit(save="no")

#APPROACH2: Full stuff
impute2 -g chr20.chunk1.gen.gz -m /home/delaneau/reference/chr20.genetic_map.txt.gz -h /home/delaneau/reference/chr20.EUR.hap.gz -l /home/delaneau/reference/chr20.EUR.legend.gz -int 1e6 3e6 -buffer 200 -o chr20.chunk1.imputed.approach2

#APPROACH2: plot performances
R
DATA1 = read.table("chr20.chunk1.imputed.approach1_info", head=TRUE)
DATA2 = read.table("chr20.chunk1.imputed.approach2_info", head=TRUE)
DATA1 = DATA1[DATA1$snp_id == "---", ]
DATA2 = DATA2[DATA2$snp_id == "---", ]
DATA1$maf = ifelse(DATA1$exp_freq_a1 < 0.5, DATA1$exp_freq_a1, 1- DATA1$exp_freq_a1)
DATA2$maf = ifelse(DATA2$exp_freq_a1 < 0.5, DATA2$exp_freq_a1, 1- DATA2$exp_freq_a1)
DATA1$bin = cut(DATA1$maf, breaks = seq(0,0.5,0.05), labels = 1:10, include.lowest = TRUE)
DATA2$bin = cut(DATA2$maf, breaks = seq(0,0.5,0.05), labels = 1:10, include.lowest = TRUE)
pdf("chr20.chunk1.imputed.approach2.pdf")
plot(by(DATA1$maf, DATA1$bin, mean), by(DATA1$info, DATA1$bin, mean), type="l", xlab="MAF", ylab="Mean info score", xlim=c(0,0.5), ylim=c(0.6, 1.0), col="grey")
points(by(DATA2$maf, DATA2$bin, mean), by(DATA2$info, DATA2$bin, mean), type="l", col="red")
legend("bottomright", legend=c("Approach 1","Approach 2"), fill=c("grey","red"))
dev.off()
quit(save="no")

#APPROACH3: Pre-phasing
impute2 -use_prephased_g -phase -known_haps_g chr20.chunk1.phased.haps -m /home/delaneau/reference/chr20.genetic_map.txt.gz -h /home/delaneau/reference/chr20.ALL.hap.gz -l /home/delaneau/reference/chr20.ALL.legend.gz -int 1e6 3e6 -buffer 200 -o chr20.chunk1.imputed.approach3

#APPROACH3: plot performances
R
DATA1 = read.table("chr20.chunk1.imputed.approach1_info", head=TRUE)
DATA2 = read.table("chr20.chunk1.imputed.approach2_info", head=TRUE)
DATA3 = read.table("chr20.chunk1.imputed.approach3_info", head=TRUE)
DATA1 = DATA1[DATA1$snp_id == "---", ]
DATA2 = DATA2[DATA2$snp_id == "---", ]
DATA3 = DATA3[DATA3$snp_id == "---", ]
DATA1$maf = ifelse(DATA1$exp_freq_a1 < 0.5, DATA1$exp_freq_a1, 1- DATA1$exp_freq_a1)
DATA2$maf = ifelse(DATA2$exp_freq_a1 < 0.5, DATA2$exp_freq_a1, 1- DATA2$exp_freq_a1)
DATA3$maf = ifelse(DATA3$exp_freq_a1 < 0.5, DATA3$exp_freq_a1, 1- DATA3$exp_freq_a1)
DATA1$bin = cut(DATA1$maf, breaks = seq(0,0.5,0.05), labels = 1:10, include.lowest = TRUE)
DATA2$bin = cut(DATA2$maf, breaks = seq(0,0.5,0.05), labels = 1:10, include.lowest = TRUE)
DATA3$bin = cut(DATA3$maf, breaks = seq(0,0.5,0.05), labels = 1:10, include.lowest = TRUE)
pdf("chr20.chunk1.imputed.approach3.pdf")
plot(by(DATA1$maf, DATA1$bin, mean), by(DATA1$info, DATA1$bin, mean), type="l", xlab="MAF", ylab="Mean info score", xlim=c(0,0.5), ylim=c(0.6, 1.0), col="black")
points(by(DATA2$maf, DATA2$bin, mean), by(DATA2$info, DATA2$bin, mean), type="l", col="grey")
points(by(DATA3$maf, DATA3$bin, mean), by(DATA3$info, DATA3$bin, mean), type="l", col="red")
legend("bottomright", legend=c("Approach 1","Approach 2", "Approach 3"), fill=c("black", "grey","red"))
dev.off()
quit(save="no")

#Proceed with the next chunk of data
bcftools view -r 20:2800000-5200000 -Oz -o chr20.chunk2.vcf.gz chr20.FINAL.vcf.gz
bcftools convert -g chr20.chunk2 chr20.chunk2.vcf.gz
shapeit --input-gen chr20.chunk2.gen.gz chr20.chunk2.samples --input-map /home/delaneau/reference/chr20.genetic_map.txt.gz --output-max chr20.chunk2.phased
impute2 -use_prephased_g -phase -known_haps_g chr20.chunk2.phased.haps -m /home/delaneau/reference/chr20.genetic_map.txt.gz -h /home/delaneau/reference/chr20.ALL.hap.gz -l /home/delaneau/reference/chr20.ALL.legend.gz -int 3e6 5e6 -buffer 200 -o chr20.chunk2.imputed.approach3
cat chr20.chunk1.imputed.approach3 chr20.chunk2.imputed.approach3 > chr20.chunkALL.imputed.approach3
cat chr20.chunk1.imputed.approach3_info chr20.chunk2.imputed.approach3_info | grep -v "position" > chr20.chunkALL.imputed.approach3_info

#Convert to VCF and filtering
bcftools convert -G chr20.chunkALL.imputed.approach3,chr20.chunk1.samples | bgzip -c > chr20.chunkALL.imputed.approach3.vcf.gz
tabix -p vcf chr20.chunkALL.imputed.approach3.vcf.gz
cat chr20.chunkALL.imputed.approach3_info | awk '{ if ($7 >= 0.4) print "20", $3; }' > chr20.chunkALL.imputed.approach3.filtering.txt
vcftools --gzvcf chr20.chunkALL.imputed.approach3.vcf.gz --positions chr20.chunkALL.imputed.approach3.filtering.txt --recode --stdout | bgzip -c > chr20.chunkALL.imputed.filtered.approach3.vcf.gz






