# This is the script for summary level data preprocessing.

### Throughout all datasets A1 is the effect allele

## BMI ##
# Note we are using an older version of BMI summary statistics from Yengo et al, Hum. Mol. Genet. 27, 3641–3649 (2018). Accessed on 5/3/2018
# The summary statistics were updated by the authors on 6/25/2018.
# Link: https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz
rm(list = ls())
library(data.table)
library(dplyr)

# Load bim file for HapMap3 SNPs extracted from 1000 Genomes European sample with 1000G MAF>0.05
bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")
# Read summary statistics
data = fread("Meta-analysis_Locke_et_al+UKBiobank_2018.txt")
commonsnps = intersect(bim$SNP, data$SNP)
# Removes SNPs with 1000G MAF<0.05, or whose alleles don't match 1000G
data = data[match(commonsnps,data$SNP),] %>% mutate(z = BETA/SE, pval = 2*pnorm(-abs(z))) %>% rename(A1 = Tested_Allele, A2 = Other_Allele)
bim = bim[match(commonsnps,bim$SNP),]
allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2
# Sample size filtering
data = data[data$N>0.67*quantile(data$N,0.9),]
data = data %>% select(chr = CHR, bp = POS, rsid = SNP, A1, A2, N, z, pval) %>%
    filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6)) # Remove MHC region and SNPs with very large effects
fwrite(data, file = "data_BMI.txt", sep = "\t")


## Height ##
# Link: https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz
rm(list = ls())
library(data.table)
library(dplyr)

bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")

data = fread("Meta-analysis_Wood_et_al+UKBiobank_2018.txt")
commonsnps = intersect(bim$SNP, data$SNP)
# Removes SNPs with 1000G MAF<0.05, or whose alleles don't match 1000G
data = data[match(commonsnps,data$SNP),] %>% mutate(z = BETA/SE, pval = 2*pnorm(-abs(z))) %>% rename(A1 = Tested_Allele, A2 = Other_Allele)
bim = bim[match(commonsnps,bim$SNP),]
allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2
# Sample size filtering
data = data[data$N>0.67*quantile(data$N,0.9),]
data = data %>% select(chr = CHR, bp = POS, rsid = SNP, A1, A2, N, z, pval) %>%
    filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6))
fwrite(data, file = "data_Height.txt", sep = "\t")


## Lipids ##
# Link: http://csg.sph.umich.edu/willer/public/lipids2013/
rm(list=ls())
rm(list = ls())
library(data.table)
library(dplyr)

traitvec = c('LDL', 'HDL', 'TG', 'TC')
for (trait in traitvec){
    data = fread(paste0('jointGwasMc_', trait, '.txt')) %>%
        filter(rsid!='.') %>%
        mutate(position = strsplit(SNP_hg19, split = ':')) %>%
        mutate(chr = as.integer(gsub('chr', '',sapply(position, function(x) x[1]))), bp = as.integer(sapply(position, function(x) x[2]))) %>%
        mutate(z = beta/se, pval = 2*pnorm(-abs(z)) , A1 = toupper(A1), A2 = toupper(A2)) %>%
        select(chr, bp, rsid, A1, A2, N, z, pval, freqA1 = Freq.A1.1000G.EUR) %>% filter(!(chr==6 & bp>26e6 & bp<34e6) & (z^2<=80) & (N>=0.67*quantile(N,0.9)) & (freqA1>=0.05 & freqA1<=0.95))
    # Merge to HapMap3 common alleles
    bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
    colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")
    commonsnps = intersect(data$rsid, bim$SNP)
    data = data[match(commonsnps,data$rsid),]
    bim = bim[match(commonsnps,bim$SNP),]
    allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
    data = data[allele_ind,]
    bim = bim[allele_ind,]
    # Align alleles to 1000G EUR
    data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
    data$A1 = bim$A1
    data$A2 = bim$A2

    fwrite(data, file = paste0("data_",trait,".txt"), sep = "\t")
}


## Blood pressure ##
# Link: http://www.nealelab.is/blog/2017/7/19/rapid-gwas-of-thousands-of-phenotypes-for-337000-samples-in-the-uk-biobank
# diastolic
rm(list=ls())
library(data.table)
library(dplyr)

# bim for 1000G EUR HapMap3 with 1000G MAF>0.05
bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")

data = fread("4079.assoc.tsv")
commonsnps = intersect(bim$SNP, data$rsid)
# Removes SNPs with 1000G MAF<0.05, or whose alleles don't match 1000G
data = data[match(commonsnps,data$rsid),]
bim = bim[match(commonsnps,bim$SNP),]
temp = strsplit(data$variant, split = ":")
data$chr = sapply(temp, function(x) x[1])
data$bp = sapply(temp, function(x) x[2])
data$A1 = sapply(temp, function(x) x[4])
data$A2 = sapply(temp, function(x) x[3])

data = data %>% mutate(z = beta/se) %>% select(chr, bp, rsid, A1, A2, N = nCompleteSamples, z, pval)

allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2
# Sample size filtering
data = data[data$N>0.67*quantile(data$N,0.9),]
data = data %>% filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6))
fwrite(data, file = "data_DBP.txt", sep = "\t")

# systolic
rm(list=ls())
library(data.table)
library(dplyr)

# bim for 1000G EUR HapMap3 with 1000G MAF>0.05
bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")

data = fread("4080.assoc.tsv")
commonsnps = intersect(bim$SNP, data$rsid)
# Removes SNPs with 1000G MAF<0.05, or whose alleles don't match 1000G
data = data[match(commonsnps,data$rsid),]
bim = bim[match(commonsnps,bim$SNP),]
temp = strsplit(data$variant, split = ":")
data$chr = sapply(temp, function(x) x[1])
data$bp = sapply(temp, function(x) x[2])
data$A1 = sapply(temp, function(x) x[4])
data$A2 = sapply(temp, function(x) x[3])

data = data %>% mutate(z = beta/se) %>% select(chr, bp, rsid, A1, A2, N = nCompleteSamples, z, pval)

allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2
# Sample size filtering
data = data[data$N>0.67*quantile(data$N,0.9),]
data = data %>% filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6))
fwrite(data, file = "data_SBP.txt", sep = "\t")


## Age at menarche 2017 ##
# Link: https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip
rm(list=ls())
library(data.table)
library(dplyr)

# bim for 1000G EUR HapMap3 with 1000G MAF>0.05
bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")
data = fread("Menarche_1KG_NatGen2017_WebsiteUpload.txt")

bim$chr_bp = paste0("chr",bim$CHR,":",bim$BP)
commonvar = intersect(data$Markername,bim$chr_bp)
bim = bim[match(commonvar,bim$chr_bp),]
data = data[match(commonvar,data$Markername),]

data$rsid = bim$SNP
data$chr = bim$CHR
data$bp = bim$BP
data = data %>% mutate(z = sign(Effect)*abs(qnorm(Pvalue/2)), N = 252514, A1 = toupper(Allele1), A2 = toupper(Allele2)) %>%
    select(chr, bp, rsid, A1, A2, N, z, pval = Pvalue)
allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2

data = data %>% filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6))
fwrite(data, file = "data_menarche2017.txt", sep = "\t")


## EduYears 2018 ##
# Link: https://www.dropbox.com/s/ho58e9jmytmpaf8/GWAS_EA_excl23andMe.txt?dl=0
rm(list=ls())
library(data.table)
library(dplyr)

# bim for 1000G EUR HapMap3 with 1000G MAF>0.05
bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")
data = fread("GWAS_EA_excl23andMe.txt")

data = data %>% mutate(N = 766345, z = Beta/SE, pval = 2*pnorm(-abs(z))) %>%
    select(chr = CHR, bp = POS, rsid = MarkerName, A1, A2, N, z, pval)

commonsnps = intersect(bim$SNP, data$rsid)
# Removes SNPs with 1000G MAF<0.05, or whose alleles don't match 1000G
data = data[match(commonsnps,data$rsid),]
bim = bim[match(commonsnps,bim$SNP),]
allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2
# Sample size filtering
data = data[data$N>0.67*quantile(data$N,0.9),]
data = data %>% filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6))
fwrite(data, file = "data_EduYears2018.txt", sep = "\t")


################
### Outcomes ###
################

## Coronary artery disease ##
# http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cad.additive.Oct2015.pub.zip
rm(list=ls())
library(data.table)
library(dplyr)

# bim for 1000G EUR HapMap3 with 1000G MAF>0.05
bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")
data = fread("gwas_ichip_meta_release/cad.add.160614.website.txt")

data = data %>% mutate(N = 60801*123504/(60801+123504), z = beta/se_dgc, pval = 2*pnorm(-abs(z))) %>%
    filter(model == "FIXED") %>%
    select(chr, bp = bp_hg19, rsid = markername, A1 = effect_allele, A2 = noneffect_allele, N, z, pval)
commonsnps = intersect(bim$SNP, data$rsid)
# Removes SNPs with 1000G MAF<0.05, or whose alleles don't match 1000G
data = data[match(commonsnps,data$rsid),]
bim = bim[match(commonsnps,bim$SNP),]

allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2

# Sample size filtering
data = data[data$N>0.67*quantile(data$N,0.9),]
data = data %>% filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6))
fwrite(data, file = "data_CAD.txt", sep = "\t")


## Major depressive disorder ##
# https://www.med.unc.edu/pgc/results-and-downloads/data-use-agreement-forms/copy2_of_mdd2018_ex23andme%20_data_download_agreement
rm(list=ls())
library(data.table)
library(dplyr)

bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")
# Error appears if reading the entire data at once
data1 = fread("MDD2018_ex23andMe", nrow = 9533408)
data3 = fread("MDD2018_ex23andMe", skip = 9533471, header = FALSE)
colnames(data3) = colnames(data1)
data = rbind(data1, data3)
rm(data1, data3)

commonsnps = intersect(bim$SNP, data$SNP)
# Removes SNPs with 1000G MAF<0.05, or whose alleles don't match 1000G
data = data[match(commonsnps,data$SNP),]
bim = bim[match(commonsnps,bim$SNP),]
data = data %>% mutate(N = Nca/(Nca+Nco)*Nco, z = log(OR)/SE, pval = 2*pnorm(-abs(z))) %>%
    select(chr = CHR, bp = BP, rsid = SNP, A1, A2, N, z, pval)

allele_ind = (data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2

# Sample size filtering
data = data[data$N>0.67*quantile(data$N,0.9),]
data = data %>% filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6))
fwrite(data, file = "data_MDD.txt", sep = "\t")


## Breast cancer ##
# http://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/gwas-icogs-and-oncoarray-summary-results/
# The fourth element in rs_id is the effect allele
rm(list = ls())
library(data.table)
library(dplyr)

# bim for 1000G EUR HapMap3 with 1000G MAF>0.05
bim = fread("1000G.EUR.QC.HM3.MERGE.maf05.bim")
colnames(bim) = c("CHR","SNP","cM","BP","A1","A2")

load("case_control_meta_result_final.rdata")
data = case_control_meta_result_final
rm(case_control_meta_result_final)
data = data[complete.cases(data),]
temp = strsplit(data$rs_id, split = ":")

data$rsid = sapply(temp, function(x) x[1])
data$A1 = sapply(temp, function(x) x[4])
data$A2 = sapply(temp, function(x) x[3])
rm(temp)

commonsnps = intersect(data$rsid, bim$SNP)

data = data[match(commonsnps,data$rsid),]
data = data %>% mutate(z = logodds_meta/sd_meta, pval = 2*pnorm(-abs(z)), N = 106571*95762/(106571+95762))
bim = bim[match(commonsnps,bim$SNP),]

data = data %>% select(chr = CHR, bp = position, rsid, A1, A2, N, z, pval)

allele_ind = ((data$A1==bim$A1 & data$A2==bim$A2) | (data$A2==bim$A1 & data$A1==bim$A2)) & (data$A1%in%c("A","T","G","C")) & (data$A2%in%c("A","T","G","C"))
data = data[allele_ind,]
bim = bim[allele_ind,]
# Align alleles to 1000G EUR
data$z[data$A1==bim$A2] = - data$z[data$A1==bim$A2]
data$A1 = bim$A1
data$A2 = bim$A2

# Sample size filtering
data = data[data$N>0.67*quantile(data$N,0.9),]
data = data %>% filter((z^2<=80) & !(chr==6 & bp>26e6 & bp<34e6))
fwrite(data, file = "data_BC.txt", sep = "\t")
