# MRMix

R package for two-sample Mendelian randomization analysis using mixture models. **For numerical stability, we recommend using summary statistics in the standardized scale, i.e. 1) for continuous phenotype, the beta and SE when both the genotypes and phenotypes are standardized to have variance 1; 2) for binary phenotype, the beta and SE when the genotypes are standardized to have variance 1**. Causal estimates are interpreted as standard deviation (SD) unit increase (continuous outcome) or log-OR (binary outcome) of Y per SD unit increase in X (assumed to be continous). If standardized-scale summary statistics are not available, users may use the `standardize` function to transform their summary statistics to the standardized scale. See Example or type `?MRMix` for more information.

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
data_std = standardize(sumstats$betahat_x, sumstats$betahat_y, sumstats$sx, sumstats$sy, contbin_x = "continuous", contbin_y = "continuous", sumstats$nx, sumstats$ny, MAF = NULL)
# MRMix analysis
est = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std)
str(est) # True causal effect is 0.2.
# Include profile matrix
est = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std, profile = TRUE)
str(est)
```

`standardize()` can be used if the phenotype is continuous and analyzed with **linear regression**, or binary and analyzed with **logistic regression**. For other types of phenotypes or other analytic methods, the user needs to standardize the data independently before using `MRMix` function. Type `?MRMix` in `R` for more details. The software has been tested on MAC OS 10.11.5 with 2.8 GHz Intel Core i7 and R version 3.5.1. Installation and the Example complete within seconds on this platform.



### More information 
Authors: Guanghao Qi (gqi1@jhu.edu) and Nilanjan Chatterjee (nchatte2@jhu.edu)

Reference: Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.

Scripts for simulations in this paper (scenarios A, B and C) are available [here](https://github.com/gqi/MRMix/tree/master/simulations). 

Scripts for data analysis in this paper are available [here](https://github.com/gqi/MRMix/tree/master/data_analysis). A brief introduction is available [here](https://github.com/gqi/MRMix/wiki).
