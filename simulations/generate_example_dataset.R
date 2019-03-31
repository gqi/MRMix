rm(list = ls())
library(MASS)

Nvec = c(5e4, 1e5, 2e5, 5e5, 1e6)
pi1vec = c(0.005, 0.01)
pthrvec = c(0.005, 5e-4, 5e-6, 5e-8)
thetavec = c(-0.2, 0, 0.2)
NxNy_ratio_vec = 1:5


N = 2e5 # Sample size of exposure X
pi1 = 0.01 # Proportion of valid instruments in the genome
pthr = 5e-8 # P-value threshold for instrument selection
theta = 0.2 # Causal effect of X on Y
NxNy_ratio = 2 # Ratio of sample size for X and Y: nx/ny

M = 2e5 # Number of independent SNPs
# Simlation model parameters
pi2 = 0.02-pi1; pi3 = 0.01; pi4 = 1-pi1-pi2-pi3
sigma2x = 5e-5; sigma2y = 5e-5; rho = 0.5

print(paste("N", N, "pthr", pthr, "pi1", pi1, "theta", theta, "NxNy_ratio", NxNy_ratio))
nx = N; ny = N/NxNy_ratio # Sample size of X and Y

repind=1
set.seed(1245*repind)

# Generate SNP indices of the 4 components
ind1 = sample(M, round(M*pi1))
ind2 = sample(setdiff(1:M,ind1), round(M*pi2))
ind3 = sample(setdiff(1:M,c(ind1,ind2)), round(M*pi3))

# Generate true effect size
beta = matrix(rep(0,2*M), ncol = 2)
beta[ind1,1] = rnorm(length(ind1), mean = 0, sd = sqrt(sigma2x))
beta[ind1,2] = theta*beta[ind1,1]
beta[ind2,] = mvrnorm(length(ind2), mu = c(0,0), Sigma = matrix(c(sigma2x,rho*sqrt(sigma2x*sigma2y),rho*sqrt(sigma2x*sigma2y),sigma2y),nrow=2))
beta[ind2,2] = beta[ind2,2] + theta*beta[ind2,1]
beta[ind3,2] = rnorm(length(ind3), mean = 0, sd = sqrt(sigma2y))

# Generate observed GWAS summary statistics
betahat.x = beta[,1] + rnorm(M, mean = 0, sd = sqrt(1/nx))
betahat.y = beta[,2] + rnorm(M, mean = 0, sd = sqrt(1/ny))
betahat.x.rev = betahat.y
betahat.y.rev = betahat.x

# Instrument selection
ind_filter = which(2*(1-pnorm(sqrt(nx)*abs(betahat.x)))<pthr)
betahat.x = betahat.x[ind_filter]
betahat.y = betahat.y[ind_filter]

sumstats = data.frame(betahat_x = 3*betahat.x, betahat_y = 2*betahat.y, sx = 3/sqrt(nx), sy = 2/sqrt(ny), nx = nx, ny = ny)
save(sumstats, file = "data/sumstats.RData")

data("sumstats", package = "MRMix")
sumstats$nx = 2e5
sumstats$ny = 1e5
sumstats$betahat_x = 3*sumstats$betahat_x
sumstats$sx2 = 3^2*sumstats$sx2
sumstats$betahat_y = 2*sumstats$betahat_y
sumstats$sy2 = 2^2*sumstats$sy2
