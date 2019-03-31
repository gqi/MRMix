#' Two-sample Mendelian randomization analysis using mixture models
#'
#' @description This function conducts Mendelian randomization analysis by fitting a mixture normal model to the bivariate effect size distribution of the exposure and the outcome. MRMix uses GWAS summary level data to identify causal effects. For numerical stability, we recommend using the function in the standardized scale, i.e. both the genotypes and phenotypes are standardized to have mean 0 and variance 1. See Examples for how to standardize summary statistics.
#'
#' @param betahat_x GWAS effect estimates of the exposure. Vector of length \code{K}, where \code{K} is the number of instruments.
#' @param betahat_y GWAS effect estimates of the outcome. Vector of length \code{K}.
#' @param sx Standard error of \code{betahat_x}. Vector of length \code{K}.
#' @param sy Standard error of \code{betahat_y}. Vector of length \code{K}.
#' @param theta_temp_vec A vector of the grid search values for the causal effect \code{theta}.
#' @param pi_init Initial value of the probability mass at the null component. See Details.
#' @param sigma_init Initial value of the variance of the non-null component. See Details.
#' @param profile Logical. Default to be \code{FALSE} and return a list of 6 containing the estimate of causal effect (\code{theta}), corresponding \code{pi0} and \code{sigma2}, standard error, z statistic and p-value. If TRUE, returns a matrix of 3 columns. The first column is \code{theta_temp_vec}. The second and third columns are the corresponding \code{pi0} and \code{sigma2} values.

#' @details The algorithm searches over a grid of possible values of the causal effect \code{theta}. For each fixed \code{theta}, it fits mixture model \code{pi0*N(0,sy2+theta^2 sx2)+(1-pi0)*N(0,sigma2)} on the residual \code{betahat_y-theta*betahat_x}. It then chooses the value of \code{theta} that have the maximum \code{pi0} as the estimate of causal effect. Under the standardized scale (genotypes and phenotypes have mean 0 and variance 1), the estimated causal effect \code{theta} is the SD unit increase in \code{Y} per SD unit increase in \code{X}. Summary statistics can be standardized by \code{beta_standardized=beta/se/sqrt(N); se_standardized=1/sqrt(N)}. See Examples for how to standardize summary statistics.
#'
#' @return A list that contains
#' \item{theta}{Estimate of causal effect.}
#' \item{pi0}{The probability mass of the null component corresponding to the estimated \code{theta}.}
#' \item{sigma2}{The variance of the non-null component corresponding to the estimated \code{theta}.}
#' \item{SE_theta}{Standard error of causal effect estimate.}
#' \item{zstat_theta}{Z-statistic for test of the causal effect estimate.}
#' \item{pvalue_theta}{P-value from the z test for the causal effect.}
#' @references
#' Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian Randomization Analysis Using Mixture Models (MRMix) for Genetic Effect-Size-Distribution Leads to Robust Estimation of Causal Effects." bioRxiv (2018): 367821.
#' @export
#' @examples
#' data("sumstats", package = "MRMix")
#' # Convert summary statistics to standardized scale
#' betahat_x = sumstats$betahat_x/sumstats$sx/sqrt(sumstats$nx)
#' betahat_y = sumstats$betahat_y/sumstats$sy/sqrt(sumstats$ny)
#' sx = 1/sqrt(sumstats$nx)
#' sy = 1/sqrt(sumstats$ny)
#' # MRMix analysis
#' est = MRMix(betahat_x, betahat_y, sx, sy)
#' data.frame(est) # True causal effect is 0.2.
MRMix = function(betahat_x, betahat_y, sx, sy, theta_temp_vec = seq(-0.49,0.5,by=0.01), pi_init = 0.6, sigma_init = 1e-5, profile = FALSE){
    sx2 = sx^2; sy2 = sy^2
    EM_res = matrix(nrow = length(theta_temp_vec), ncol = 3)
    colnames(EM_res) = c("theta", "pi0", "sigma2")

    for (i in 1:length(theta_temp_vec)){
        theta_temp = theta_temp_vec[i]
        # Initial values
        pi0 = pi_init
        sigma0 = sigma_init
        f0 = 1/sqrt(2*pi*(sy2+theta_temp^2*sx2))*exp(-0.5/(sy2+theta_temp^2*sx2)*(betahat_y-theta_temp*betahat_x)^2)
        f0[f0<1e-300] = 1e-300
        # Convergence criterion
        pi0_ck = pi0
        sigma0_ck = sigma0
        iter_ck = 0
        # Fit mixture models via EM algorithm
        for (iter in 1:50000){
            f1 = 1/sqrt(2*pi*sigma0)*exp(-0.5/sigma0*(betahat_y-theta_temp*betahat_x)^2)
            f1[f1<1e-300] = 1e-300
            loglkl = sum(log(pi0*f0+(1-pi0)*f1))
            pt = pi0*f0/(pi0*f0+(1-pi0)*f1)
            pt[pt>0.9999999] = 0.9999999
            pi0 = mean(pt)
            sigma0 = sum((1-pt)*(betahat_y-theta_temp*betahat_x)^2)/(length(betahat_x)-sum(pt))

            if (pi0<0.0001) pi0 = 0.0001
            if (pi0>0.9999) pi0 = 0.9999
            if (sigma0>0.01) sigma0=0.01
            if (sigma0<1e-7) sigma0=1e-7

            if (iter%%1000==0){
                if ((abs(pi0-pi0_ck)/pi0<1e-4) & (abs(sigma0-sigma0_ck)/sigma0<1e-4)){
                    break
                } else{
                    iter_ck = iter
                    pi0_ck = pi0
                    sigma0_ck = sigma0
                }
            }
        }
        EM_res[i,] = c(theta_temp, pi0, sigma0)
    }

    if (profile==TRUE){
        return(EM_res)
    } else{
        ind = EM_res[,2]==max(EM_res[,2])
        est = list(theta = mean(EM_res[ind,1]),
                   pi0 = max(EM_res[,2]),
                   sigma2 = mean(EM_res[ind,3]))
        est$SE_theta = MRMix_se(betahat_x, betahat_y, sx, sy, est$theta, est$pi0, est$sigma2) # Standard error
        est$zstat_theta = est$theta/est$SE_theta
        est$pvalue_theta = 2*pnorm(-abs(est$zstat_theta))
        return(est)
    }
}
