# MRMix

Two-sample Mendelian randomization analysis using mixture models

### Installation
```
library(devtools)
install_github("gqi/MRMix")
```

### Example
```
data("sumstats", package = "MRMix")
est = MRMix(sumstats$betahat_x, sumstats$betahat_y, sumstats$sx2, sumstats$sy2)
se = MRMix_se(sumstats$betahat_x, sumstats$betahat_y, sumstats$sx2, sumstats$sy2, est$theta, est$pi0, est$sigma2)
```
Type `?MRMix` in `R` for more details.


### More information 
Authors: Guanghao Qi (gqi1@jhu.edu) and Nilanjan Chatterjee (nchatte2@jhu.edu)

Reference: Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian Randomization Analysis Using Mixture Models (MRMix) for Genetic Effect-Size-Distribution Leads to Robust Estimation of Causal Effects." bioRxiv (2018): 367821.
