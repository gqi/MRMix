# MRMix

R package for two-sample Mendelian randomization analysis using mixture models. For numerical stability, we recommend using the function in the standardized scale, i.e. both the genotypes and phenotypes are standardized to have mean 0 and variance 1. See Example for how to obtain summary statistics in the standardized scale

### System requirements

MRMix can be used on any operating system. R needs to be installed. Package `devtools` is required for the installation.

### Installation
```
library(devtools)
install_github("gqi/MRMix")
```

### Example
```
library(MRMix)
data("sumstats", package = "MRMix")
# Convert summary statistics to standardized scale
betahat_x = sumstats$betahat_x/sumstats$sx/sqrt(sumstats$nx)
betahat_y = sumstats$betahat_y/sumstats$sy/sqrt(sumstats$ny)
sx = 1/sqrt(sumstats$nx)
sy = 1/sqrt(sumstats$ny)
# MRMix analysis
est = MRMix(betahat_x, betahat_y, sx, sy)
data.frame(est)

#   theta       pi0       sigma2   SE_theta zstat_theta pvalue_theta
# 1  0.21 0.4602256 8.998972e-05 0.02449794    8.572151 1.015702e-17
```
Type `?MRMix` and `?MRMix_se` in `R` for more details. The software has been tested on MAC OS 10.11.5 with R version 3.5.1. Installation and demo should run within seconds.



### More information 
Authors: Guanghao Qi (gqi1@jhu.edu) and Nilanjan Chatterjee (nchatte2@jhu.edu)

Reference: Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian Randomization Analysis Using Mixture Models (MRMix) for Genetic Effect-Size-Distribution Leads to Robust Estimation of Causal Effects." bioRxiv (2018): 367821.
