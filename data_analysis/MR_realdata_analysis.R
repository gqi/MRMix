# This is the script for conducting MR analysis used in Qi and Chatterjee, Nature Communications (2018)
rm(list=ls())
library(data.table)
library(dplyr)
library(MendelianRandomization)
library(MRMix)

# Change trait_x and trait_y to analyze other traits
trait_x = "BMI"
trait_y = "CAD"
pthr = 5e-8

print(paste(trait_x, trait_y, pthr))

data_x = read.table(paste0("data_",trait_x,".txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
data_y = read.table(paste0("data_",trait_y,".txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t")

commonsnps = intersect(data_x$rsid, data_y$rsid)
data_x = data_x[match(commonsnps,data_x$rsid),]
data_y = data_y[match(commonsnps,data_y$rsid),]

# Self-computed LD score from 1000 Genomes European sample using HapMap3 SNPs
load("ldscore_thr_r0.1.rda")
ldscore = ldscore[match(commonsnps,names(ldscore))]

z2x = data_x$z^2
ldscx = ldscore*data_x$N
h2x = coef(lm(z2x~ldscx))[2]
zxy = data_x$z*data_y$z
ldscx = ldscore*sqrt(data_x$N)*sqrt(data_y$N)
rhoxy = coef(lm(zxy~ldscx))[2]
ldsc_est = rhoxy/h2x

source("ld_clump_MRsim.R")
load("LDinfo_with_R2.rda")
data_x.temp = data_x %>% select(chr, rsid, bp, z, pval)
data_x.pruned = ld_clump_MRsim(data_x.temp, LDinfo_with_R2, pthr = pthr, r2thr = 0.1)
rm(data_x.temp)

ind = match(data_x.pruned$SNP, data_x$rsid)
data_x = data_x[ind,]
data_y = data_y[ind,]
numIV = length(ind)

mr.obj = mr_input(bx = data_x$z/sqrt(data_x$N), bxse = 1/sqrt(data_x$N),
                  by = data_y$z/sqrt(data_y$N), byse = 1/sqrt(data_y$N),
                  snps = data_x.pruned$SNP)
IVW = mr_ivw(mr.obj)
med = mr_median(mr.obj)
mod = mr_mbe(mr.obj, iterations = 100)
egger = mr_egger(mr.obj)

betahat_x = data_x$z/sqrt(data_x$N)
betahat_y = data_y$z/sqrt(data_y$N)
nx = data_x$N
ny = data_y$N

MRres = list()
MRres$ldsc_est = ldsc_est
MRres$IVW_est = IVW$Estimate
MRres$IVW_sd = IVW$StdError
MRres$median_est = med$Estimate
MRres$median_sd = med$StdError
MRres$mode_est = mod$Estimate
MRres$mode_sd = mod$StdError
MRres$egger_est = egger$Estimate
MRres$egger_sd = egger$StdError.Est

theta_temp_vec = seq(-0.99,1,by=0.01) # 200 values in total

res = MRMix(betahat_x, betahat_y, 1/sqrt(nx), 1/sqrt(ny), theta_temp_vec, pi_init = 0.6, sigma_init = 1e-5)

MRres$MRmix_est = res$theta
MRres$MRmix_se_analyt = res$SE_theta
MRres$pi0_est = res$pi0
MRres$sigma2_est = res$sigma2
MRres$numIV = numIV

save(MRres, file = paste0("MRres_",trait_x,"_",trait_y,"_pthr",pthr,".rda"))
