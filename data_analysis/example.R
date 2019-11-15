# Analysis
rm(list=ls())
library(MRMix)
library(dplyr)
data("BMI15_MDD18", package = "MRMix")

# load("BMI15_MDD18.RData")
data = BMI15 %>% inner_join(MDD18,by="SNP") %>%
    mutate(beta.y = log(OR)) %>%
    select(SNP, chr = Chr, bp, EA.x = effect_allele, NEA.x = other_allele,
           nx = N, beta.x = beta, se.x = se, EAF.x = EAF,
           EA.y = A1, NEA.y = A2, beta.y, se.y = SE, EAF.y = FRQ_U_113154)
# Check if the alleles are the same for X and Y
data %>% filter((EA.x==EA.y&NEA.x==NEA.y) | (EA.x==NEA.y&NEA.x==EA.y)) %>% nrow()
# 0 - Alleles are the same.

# Check the allele frequency of palindromic SNPs.
# Results show that the are coded with respect to the same strand.
data %>% filter((EA.x=="A"&NEA.x=="T") | (EA.x=="T"&NEA.x=="A") | (EA.x=="G"&NEA.x=="C") | (EA.x=="C"&NEA.x=="G")) %>%
    select(SNP, EA.x, NEA.x, EAF.x, EA.y, NEA.y, EAF.y)

# Harmonize data
data = data %>% mutate(beta.y = ifelse(EA.x==EA.y,beta.y,-beta.y), MAF = pmin(EAF.x,1-EAF.x))
#
data_std = with(data, standardize(beta.x,beta.y,se.x,se.y,xtype = "continuous", ytype = "binary", nx = nx, ny = NULL, MAF = MAF))

# MRMix analysis
res = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std, profile = TRUE)
str(res)

# List of 7
# $ theta       : num 0.33
# $ pi0         : num 0.462
# $ sigma2      : num 0.000109
# $ SE_theta    : num 0.133
# $ zstat_theta : num 2.48
# $ pvalue_theta: num 0.013
# $ profile_mat : num [1:201, 1:3] -1 -0.99 -0.98 -0.97 -0.96 -0.95 -0.94 -0.93 -0.92 -0.91 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:3] "theta" "pi0" "sigma2"

# Profile plot
plot(res$profile_mat[,"theta"], res$profile_mat[,"pi0"], xlab = "theta", ylab = "pi0", type = "l")
abline(v = res$theta, lty = 2, col = "red")
