#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Where does the data come from?
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Retrieve the data:
# 1000 Genomes website: http://www.internationalgenome.org/
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# download chr2
ALL.chr2.*
# file > README_vcf_info_annotation.20141104
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/README_vcf_info_annotation.20141104
# We added AF for each super-population to the INFO field.
#East Asian		EAS
#South Asian	SAS
#African		AFR
#European		EUR
#American		AMR
# says "The super population assignment for each sample can be found in integrated_call_samples_v3.20130502.ALL.panel"
# file > integrated_call_samples_v3.20130502.ALL.panel

# We are interested in distinct population sets:
# EUR
# GBR
# EAS
# CHB
# AFR
# YRI

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Calculate Tajima's D and FST 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# NOTE: I have used vcftools, but other software can compute these and other statistics (e.g. bcftools)

# Using a set of unix commands: get sample IDs for each population set of interest

# Calculate FST and Tajima's D on Chr 2 for EUR, AFR and EAS:
# using vcftools, what do the followings command mean?

vcftools --gzvcf ${INFILE} --remove-indels --recode --recode-INFO-all --out EMBOrun01_chr2SNPs_only

vcftools --vcf EMBOrun01_chr2SNPs_only.recode.vcf --out EMBOrun01_FST_EasEur --chr 2 --weir-fst-pop ${EAScodes} --weir-fst-pop ${EURcodes} 

vcftools --vcf EMBOrun01_chr2SNPs_only.recode.vcf --out EMBOrun01_TajD_EUR --chr 2 --keep ${EURcodes} --TajimaD 1000

# Produce the following files:
#EMBOrun01_FST_EasEur.weir.fst
#EMBOrun01_FST_AfrEas.weir.fst
#EMBOrun01_FST_AfrEur.weir.fst
#EMBOrun01_TajD_EUR.Tajima.D
#EMBOrun01_TajD_AFR.Tajima.D
#EMBOrun01_TajD_EAS.Tajima.D

# Read *Tajima.D files with TajimasDpractical.R.sh
# Read *.weir.fst files with Bersaglieri_etal_2004_FSTpractical.sh

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Calculate iHS 
# SELSCAN
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# selscan requires data in a certain format
# In selscan manual, find the input file format requirements

# What do the following commands do?
gunzip ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
bgzip ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf

# What are the loci of interest?
# LCT and EDAR, specific SNP?
# example:
# LCT (−13910*T): rs4988235
# location: 136,608,646
#136608646 - 1000000 = 135608646
#136608646 + 1000000 = 137608646
# How do you know the location?

vcffilter -f "VT = SNP" -f "AC > 1" -f "QUAL > 20" -f "DP > 40" -f "DP < 40000" -r 2:135545410-137545410 ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > ALL.chr2.LCT.tmp.vcf
bcftools view -h ALL.chr2.LCT.tmp.vcf
grep -v MULTI_ALLELIC ALL.chr2.LCT.tmp.vcf > ALL.chr2.LCT.vcf
bcftools view -h ALL.chr2.LCT.vcf

vcfkeepsamples ALL.chr2.LCT.vcf NA18525 NA18526 NA18528 NA18530 NA18531 NA18532 NA18533 NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 NA18557 NA18558 NA18559 NA18560 NA18561 NA18562 NA18563 NA18564 NA18565 NA18566 NA18567 NA18570 NA18571 NA18572 NA18573 NA18574 NA18577 NA18579 NA18582 NA18591 NA18592 NA18593 NA18595 NA18596 NA18597 NA18599 NA18602 NA18603 NA18605 NA18606 NA18608 NA18609 NA18610 NA18611 NA18612 NA18613 NA18614 NA18615 NA18616 NA18617 NA18618 NA18619 NA18620 NA18621 NA18622 NA18623 NA18624 NA18625 NA18626 NA18627 NA18628 NA18629 NA18630 NA18631 NA18632 NA18633 NA18634 NA18635 NA18636 NA18637 NA18638 NA18639 NA18640 NA18641 NA18642 NA18643 NA18644 NA18645 NA18646 NA18647 NA18648 NA18740 NA18745 NA18747 NA18748 NA18749 NA18757 > CHB.chr2.LCT.tmp.vcf
bcftools view -h CHB.chr2.LCT.tmp.vcf
sed -i -e 's/\t\./\t0|0/g' CHB.chr2.LCT.vcf
bcftools view -h CHB.chr2.LCT.vcf

mkdir CHBMap
tar -C CHBMap -xvf CHB_omni_recombination_20130507.tar
zcat CHBMap/CHB/CHB-2-final.txt.gz > CHBgenetic_map_chrom2.map
rm -rf CHBMap
rm CHB_omni_recombination_20130507.tar
#------------
# What follows is an R code to keep only positions of interest from a recombination map file 
#------------
vcf=read.table(fvcf, stringsAsFactors=F)
rs=vcf[,3]
pos=as.numeric(vcf[,2])
rm(vcf)
map=read.table(fmap, head=T)
ind=match(pos, map$Position)
maps=map$Map[ind]
app=approx(x=c(pos[1]-200,pos,pos[length(pos)]+200), y=c(min(maps,na.rm=T)-0.001, maps, max(maps,na.rm=T)+0.001), xout=c(pos[1]-200,pos,pos[length(pos)]+200))
newmap=cbind(rep(2, length(rs)), rs, app$y[2:(length(app$y)-1)], pos)
write.table(newmap, file='CHBgenetic_map_GRCh37_chr2LCT.map', quote=F, sep=" ", col.names=F, row.names=F)
# creates the file CHBgenetic_map_GRCh37_chr2LCT.map
#------------

selscan --ihs --vcf CHB.chr2.LCT.vcf --map CHBgenetic_map_chrom2.map --out CHB.LCT.ihs

# Create the following files:
# GBR.EDAR.ihs.ihs.out
# CHB.EDAR.ihs.ihs.out
# GBR.LCT.ihs.ihs.out
# CHB.LCT.ihs.ihs.out

# Read *ihs.out files with iHSpractical.sh

