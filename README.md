# MRMix

R package for two-sample Mendelian randomization analysis using mixture models. **For numerical stability, we recommend using summary statistics in the standardized scale, i.e. 1) for continuous phenotype, the beta and SE when both the genotypes and phenotypes are standardized to have variance 1; 2) for binary phenotype, the beta and SE when the genotypes are standardized to have variance 1**. Causal estimates are interpreted as standard deviation (SD) unit increase (continuous outcome) or log-OR (binary outcome) of Y per SD unit increase in X (assumed to be continous). If standardized-scale summary statistics are not available, users may use z-statistics and effective sample size to transform their summary statistics to the standardized scale. See Example or type `?MRMix` for more information.

### System requirements

MRMix can be used on any operating system. R needs to be installed. Package `devtools` is required for the installation.

### Installation
```
devtools::install_github("gqi/MRMix")
```

### Example
```
library(MRMix)
data("sumstats", package = "MRMix")
# Convert summary statistics to standardized scale
# beta_standardized=beta/se/sqrt(N); se_standardized=1/sqrt(N).
betahat_x = sumstats$betahat_x/sumstats$sx/sqrt(sumstats$nx)
betahat_y = sumstats$betahat_y/sumstats$sy/sqrt(sumstats$ny)
sx = 1/sqrt(sumstats$nx)
sy = 1/sqrt(sumstats$ny)
# MRMix analysis
est = MRMix(betahat_x, betahat_y, sx, sy)
data.frame(est) # True causal effect is 0.2.

#   theta       pi0       sigma2   SE_theta zstat_theta pvalue_theta
# 1  0.21 0.4602256 8.998972e-05 0.02449794    8.572151 1.015702e-17
```

If the phenotype is continuous and analyzed with **linear regression**, or binary and analyzed with **logistic regression**, the standardized-scale standard error is 1/sqrt(N) and beta is z-statistics/sqrt(N). For continuous phenotypes, N is the total sample size and for binary phenotypes N is the effective sample size, calculated as Ncase*Ncontrol/(Ncase+Ncontrol). This relationship may not hold for non-standard analysis methods and standardization would require the phenotypic variance. Type `?MRMix` in `R` for more details. The software has been tested on MAC OS 10.11.5 with 2.8 GHz Intel Core i7 and R version 3.5.1. Installation and the Example complete within seconds on this platform.



### More information 
Authors: Guanghao Qi (gqi1@jhu.edu) and Nilanjan Chatterjee (nchatte2@jhu.edu)

Reference: Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.

Scripts for simulations in this paper (scenarios A, B and C) are available [here](https://github.com/gqi/MRMix/tree/master/simulations). 

Scripts for data analysis in this paper are available [here](https://github.com/gqi/MRMix/tree/master/data_analysis). A brief introduction is available [here](https://github.com/gqi/MRMix/wiki).
