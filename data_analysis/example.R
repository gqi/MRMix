# Prepare data
rm(list=ls())
library(data.table)
library(dplyr)
BMI15 = fread("/dcl01/chatterj/data/gqi/multitrait_effsize/data/raw_sumstats/BMI2015_97tophits.csv")
BMI15 = BMI15 %>% mutate(bp = as.integer(gsub(",","",bp)), N = as.numeric(gsub(",","",N)))

MDD18 = fread("/dcl01/chatterj/data/Rawdata-summarydata/MDD_2018/MDD2018_ex23andMe", fill=TRUE) %>%
    filter(SNP %in% BMI15$SNP)
save(BMI15, MDD18, file = "BMI15_MDD18.RData")

# Analysis
rm(list=ls())
library(MRMix)
load("BMI15_MDD18.RData")
data = BMI15 %>% inner_join(MDD18,by="SNP") %>%
    mutate(beta.y = log(OR), MAF = pmin(EAF,1-EAF)) %>%
    select(SNP, chr = Chr, bp, EA.x = effect_allele, NEA.x = other_allele, nx = N, beta.x = beta, se.x = se,
           EA.y = A1, NEA.y = A2, beta.y, se.y = SE, MAF)

data = data %>% filter((EA.x==EA.y&NEA.x==NEA.y) | (EA.x==NEA.y&NEA.x==EA.y)) %>%
    mutate(beta.y = ifelse(EA.x==EA.y,beta.y,-beta.y))

data_std = with(data, standardize(beta.x,beta.y,se.x,se.y,contbin_x = "continuous", contbin_y = "binary", nx = nx, ny = NULL, MAF = MAF))
res = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std)
str(res)


