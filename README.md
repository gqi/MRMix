# MRMix

R package for two-sample Mendelian randomization analysis using mixture models.

### System requirements

MRMix can be used on any operating system. R needs to be installed.

### Installation
```
library(devtools)
install_github("gqi/MRMix")
```

### Example
```
data("sumstats", package = "MRMix")
est = MRMix(sumstats$betahat_x, sumstats$betahat_y, sumstats$sx2, sumstats$sy2)
se = MRMix_se(sumstats$betahat_x, sumstats$betahat_y, sumstats$sx2, sumstats$sy2, est$theta, est$pi0, est$sigma2) # Standard error
print(est)
print(se)

## Expected output
# $theta
# [1] 0.22

# $pi0
# [1] 0.509395

# $sigma2
# [1] 9.04497e-05

## Standard error
# [1] 0.02084057
```
Type `?MRMix` and `?MRMix_se` in `R` for more details. Installation and demo should run within seconds.


### More information 
Authors: Guanghao Qi (gqi1@jhu.edu) and Nilanjan Chatterjee (nchatte2@jhu.edu)

Reference: Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian Randomization Analysis Using Mixture Models (MRMix) for Genetic Effect-Size-Distribution Leads to Robust Estimation of Causal Effects." bioRxiv (2018): 367821.
