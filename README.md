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

##### Step 1. Load package and data
```
library(MRMix)
data("BMI15_MDD18", package = "MRMix")
```
##### Step 2. Merge data
Merge datasets for the exposure and the outcome, compute `beta.y` from odds ratio (`OR`), compute minor allele frequency (`MAF`) from effect allele frequency (`EAF`). Keep the following variables: SNP ID (`SNP`), chromosome (`chr`), base pair (`bp`), alleles (`EA.x, NEA.x, EA,y, NEA.y`), sample size of the study associated with X (`nx`), estimate of genetic effect (`beta.x, beta.y`) and standard errors (`se.x, se.y`), MAF (`MAF`).
```
data = BMI15 %>% inner_join(MDD18,by="SNP") %>%
    mutate(beta.y = log(OR), MAF = pmin(EAF,1-EAF)) %>%
    select(SNP, chr = Chr, bp, EA.x = effect_allele, NEA.x = other_allele, 
    nx = N, beta.x = beta, se.x = se,
    EA.y = A1, NEA.y = A2, beta.y, se.y = SE, MAF)
```

##### Step 3. Harmonize data
Remove the SNPs of which the alleles do not match between the exposure and outcome datasets. Flip the sign of `beta.y` if the effect allele in the study associated with X is the non-effect allele in the study associated with Y.
```
data = data %>% filter((EA.x==EA.y&NEA.x==NEA.y) | (EA.x==NEA.y&NEA.x==EA.y)) %>%
    mutate(beta.y = ifelse(EA.x==EA.y,beta.y,-beta.y))
```

##### Step 4. Standardize data
Since the exposure (BMI) is continuous, the summary statistics are standardized as z-statistics rescaled by the sample size. Since the outcome (MDD) is binary, the summary statistics are standardized by genotypic variance calculated as `2*MAF*(1-MAF)` under Hardy-Weinberg equilibrium.
```
data_std = with(data, standardize(beta.x,beta.y,se.x,se.y,contbin_x = "continuous", contbin_y = "binary", nx = nx, ny = NULL, MAF = MAF))
```
`standardize()` can be used if the summary statistics are estimates from linear or logistic regression. For other types of phenotypes or other analytic methods, the users need to standardize the data independently. **If your data have been standardize, skip this step and proceed to Step 5.**

##### Step 5. MRMix analysis
```
res = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std)
str(res)
```
The results show that higher BMI increases the risk of MDD.

### More information 
Authors: Guanghao Qi (gqi1@jhu.edu) and Nilanjan Chatterjee (nchatte2@jhu.edu)

References: 

1. Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.
2. Qi, Guanghao, and Nilanjan Chatterjee. "A Comprehensive Evaluation of Methods for Mendelian Randomization Using Realistic Simulations of Genome-wide Association Studies." bioRxiv (2019): 702787.

Scripts for simulations in Reference 1 (scenarios A, B and C) are available [here](https://github.com/gqi/MRMix/tree/master/simulations). 

Scripts for data analysis in Reference 1 are available [here](https://github.com/gqi/MRMix/tree/master/data_analysis). A brief introduction is available [here](https://github.com/gqi/MRMix/wiki).

Scripts for simulation studies in Reference 2 are available [here](https://github.com/gqi/MR_comparison_simulations).

The software has been tested on MAC OS 10.11.5 with 2.8 GHz Intel Core i7 and R version 3.5.1.
