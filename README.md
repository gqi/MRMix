# MRMix

R package for two-sample This function conducts Mendelian randomization analysis using an underlying mixture model incorporating a fraction of the genetic instruments to have direct effect on the outcome (horizontal pleiotropy). MRMix takes GWAS summary statistics as inputs to estimate causal effects of one trait on another. For stability of the method, we recommend using summary statistics in the standardized scale: 1) For both binary and continuous traits, summary-statistics should be standardized by genotypic variance; 2) In addition, for continuous phenotype, summary-statistics should be standardized by phenotypic variance. If the data are not in the standardized scale, users may use the standardize function to standardize their data. See Example or type `?MRMix` for more information.

### System requirements

MRMix can be used on any operating system. R needs to be installed. Package `devtools` is required for the installation.

### Installation
```
devtools::install_github("gqi/MRMix")
```

### Example: body mass index (BMI) and major depressive disorder (MDD)

##### Data source
* Exposure: 97 SNPs from GIANT Consortium 2015 GWAS on BMI ([Locke et al, 2015](https://www.nature.com/articles/nature14177)).
* Outcome: Psychiatric Genomics Consortium 2018 GWAS on MDD ([Wray et al, 2018](https://www.nature.com/articles/s41588-018-0090-3)). To reduce file size, we only included the SNPs that appear in the exposure dataset.

#### Step 1. Load package and data
```
library(MRMix)
library(dplyr)
data("BMI15", package = "MRMix")
data("MDD18", package = "MRMix")
```
#### Step 2. Merge data
Merge datasets for the exposure and the outcome, compute `beta.y` from odds ratio (`OR`). Keep the following variables: SNP ID (`SNP`), chromosome (`chr`), base pair (`bp`), effect and non-effect alleles (`EA.x, NEA.x, EA,y, NEA.y`), sample size of the study associated with X (`nx`), estimate of genetic effect (`beta.x, beta.y`) and standard errors (`se.x, se.y`), effect allele frequency (`EAF.x`, `EAF.y`). Remove the SNPs of which the alleles do not match between the exposure and outcome datasets. 
```
data = BMI15 %>% inner_join(MDD18,by="SNP") %>%
    mutate(beta.y = log(OR)) %>%
    select(SNP, chr = Chr, bp, EA.x = effect_allele, NEA.x = other_allele,
           nx = N, beta.x = beta, se.x = se, EAF.x = EAF,
           EA.y = A1, NEA.y = A2, beta.y, se.y = SE, EAF.y = FRQ_U_113154) 
# Check if the alleles are the same for X and Y
data %>% filter(!((EA.x==EA.y&NEA.x==NEA.y) | (EA.x==NEA.y&NEA.x==EA.y))) %>% nrow()
# 0 - Alleles are the same.
```

#### Step 3. Harmonize data
Check the allele frequency of palindromic SNPs. Results show that the are coded with respect to the same strand.
```
data %>% filter((EA.x=="A"&NEA.x=="T") | (EA.x=="T"&NEA.x=="A") | (EA.x=="G"&NEA.x=="C") | (EA.x=="C"&NEA.x=="G")) %>%
    select(SNP, EA.x, NEA.x, EAF.x, EA.y, NEA.y, EAF.y)
    
#          SNP EA.x NEA.x EAF.x EA.y NEA.y EAF.y
# 1  rs1558902    A     T 0.415    A     T 0.410
# 2  rs4256980    G     C 0.646    C     G 0.349
# 3  rs9641123    C     G 0.429    C     G 0.410
# 4 rs17001654    G     C 0.153    C     G 0.841
# 5  rs9914578    G     C 0.211    C     G 0.795
```

Flip the sign of `beta.y` if the effect allele in the study associated with X is the non-effect allele in the study associated with Y.
```
data = data %>% mutate(beta.y = ifelse(EA.x==EA.y,beta.y,-beta.y), MAF = pmin(EAF.x,1-EAF.x))
```

#### Step 4. Standardize data
Since the exposure (BMI) is continuous, the summary statistics are standardized as z-statistics rescaled by the sample size. Since the outcome (MDD) is binary, the summary statistics are standardized by genotypic variance calculated as `2*MAF*(1-MAF)` under Hardy-Weinberg equilibrium.
```
data_std = with(data, standardize(beta.x,beta.y,se.x,se.y,xtype = "continuous", ytype = "binary", nx = nx, ny = NULL, MAF = MAF))
```
`standardize()` can be used if the summary statistics are estimates from linear or logistic regression. For other types of phenotypes or other analytic methods, the users need to standardize the data independently. **If your data have been standardized, skip this step and proceed to Step 5.**

#### Step 5. MRMix analysis
```
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
```
The results show that higher BMI increases the risk of MDD.

#### Step 6. Profile plot
Plot `pi0` against `theta`. A smooth curve with a clear maximum indicates stable performance.
```
plot(res$profile_mat[,"theta"], res$profile_mat[,"pi0"], xlab = "theta", ylab = "pi0", type = "l")
abline(v = res$theta, lty = 2, col = "red")
```



### More information 
Authors: Guanghao Qi (gqi1@jhu.edu) and Nilanjan Chatterjee (nchatte2@jhu.edu)

References: 

1. Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.
2. Qi, Guanghao, and Nilanjan Chatterjee. "A Comprehensive Evaluation of Methods for Mendelian Randomization Using Realistic Simulations of Genome-wide Association Studies." bioRxiv (2019): 702787.

Scripts for simulations in Reference 1 (scenarios A, B and C) are available [here](https://github.com/gqi/MRMix/tree/master/simulations). 

Scripts for data analysis in Reference 1 are available [here](https://github.com/gqi/MRMix/tree/master/data_analysis). A brief introduction is available [here](https://github.com/gqi/MRMix/wiki).

Scripts for simulation studies in Reference 2 are available [here](https://github.com/gqi/MR_comparison_simulations).

The software has been tested on MAC OS 10.11.5 with 2.8 GHz Intel Core i7 and R version 3.5.1.
