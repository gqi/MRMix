# MRMix

R package for two-sample Mendelian randomization analysis using mixture models. MRMix takes GWAS summary statistics as inputs to estimate causal effects of one trait on another. **We recommend using summary statistics in the standardized scale: 1) for continuous phenotypes, the data should be standardized w.r.t. genotypic and phenotypic variance; 2) for binary phenotypes, the data should be standardized w.r.t. genotypic variance**. Causal estimates are interpreted as standard deviation (SD) unit increase in the mean (continuous outcome) or log-OR (binary outcome) of Y per SD unit increase in X (if continous). If standardized-scale summary statistics are not available, users may use the `standardize` function to standardize their data. See Example or type `?MRMix` for more information.

### System requirements

MRMix can be used on any operating system. R needs to be installed. Package `devtools` is required for the installation.

### Installation
```
devtools::install_github("gqi/MRMix")
```

### Example

If the summary statistics have been standardized
```
library(MRMix)
data("sumstats_std", package = "MRMix") # sumstats_std has been standardized

est = MRMix(sumstats_std$betahat_x_std, sumstats_std$betahat_y_std, sumstats_std$sx_std, sumstats_std$sy_std)
str(est)
```

If the data have not been standardized
```
library(MRMix)
data("sumstats", package = "MRMix")
# Convert summary statistics to standardized scale
data_std = standardize(sumstats$betahat_x, sumstats$betahat_y, sumstats$sx, sumstats$sy, contbin_x = "continuous", contbin_y = "continuous", sumstats$nx, sumstats$ny, MAF = NULL)
# MRMix analysis
est = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std)
str(est) # True causal effect is 0.2.
# Include profile matrix
est = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std, profile = TRUE)
str(est)
```

`standardize()` can be used if the phenotype is continuous and analyzed with **linear regression**, or binary and analyzed with **logistic regression**. For other types of phenotypes or other analytic methods, the users need to standardize the data independently before using `MRMix`. Type `?MRMix` in `R` for more details. The software has been tested on MAC OS 10.11.5 with 2.8 GHz Intel Core i7 and R version 3.5.1. Installation and the Example complete within seconds on this platform.



### More information 
Authors: Guanghao Qi (gqi1@jhu.edu) and Nilanjan Chatterjee (nchatte2@jhu.edu)

Reference: Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.

Scripts for simulations in this paper (scenarios A, B and C) are available [here](https://github.com/gqi/MRMix/tree/master/simulations). 

Scripts for data analysis in this paper are available [here](https://github.com/gqi/MRMix/tree/master/data_analysis). A brief introduction is available [here](https://github.com/gqi/MRMix/wiki).
