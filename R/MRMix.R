#' Two-sample Mendelian randomization analysis using mixture models
#'
#' @description This function conducts Mendelian randomization analysis using the mixture model approach. MRMix takes GWAS summary level data (beta and standard error (SE)) as inputs to estimate causal effects of one trait on another. For numerical stability, we recommend using summary statistics in the standardized scale, i.e. 1) for continuous phenotype, the beta and SE when both the genotypes and phenotypes are standardized to have variance 1; 2) for binary phenotype, the beta and SE when the genotypes are standardized to have variance 1. Causal estimates are interpreted as standard deviation (SD) unit increase (continuous outcome) or log-OR (binary outcome) of Y per SD unit increase in X (assumed to be continous). If standardized-scale summary statistics are not available, users may use the \code{standardize()} function to transform their summary statistics to the standardized scale. See Details and Examples for more information.
#'
#' @param betahat_x GWAS effect estimates of the exposure, recommended to be in standardized scale. Vector of length \code{K}, where \code{K} is the number of instruments.
#' @param betahat_y GWAS effect estimates of the outcome, recommended to be in standardized scale. Vector of length \code{K}.
#' @param sx Standard error of \code{betahat_x}, recommended to be in standardized scale. Vector of length \code{K}.
#' @param sy Standard error of \code{betahat_y}, recommended to be in standardized scale. Vector of length \code{K}.
#' @param theta_temp_vec A vector of the grid search values for the causal effect \code{theta}. Default to be \code{seq(-1,1,by=0.01)}. Users may adjust the grid if larger effects are possible.
#' @param pi_init Initial value of the probability mass at the null component. Default to be 0.6. See Details.
#' @param sigma_init Initial value of the variance of the non-null component. Default to be 1e-5. See Details.
#' @param profile Whether to include the profile matrix. Default to be \code{FALSE}. If \code{TRUE}, include the profile matrix in the output. See Value \code{profile_mat} for details.
#'
#' @details The algorithm searches over a grid of possible values of the causal effect \code{theta}. For each fixed \code{theta}, it fits mixture model \code{pi0*N(0,sy^2+theta^2*sx^2)+(1-pi0)*N(0,sigma2)} on the residual \code{betahat_y-theta*betahat_x}. It then chooses the value of \code{theta} that leads to the maximum \code{pi0} as the estimate of causal effect. Under the standardized scale, the estimated causal effect \code{theta} is the SD unit increase or logOR of \code{Y} per SD unit increase in \code{X}. Summary statistics can be standardized using the \code{standardize()} function if phenotypes are continuous and analyzed with linear regression, or binary and analyzed with logistic regression.
#'
#' @return A list that contains
#' \item{theta}{Estimate of causal effect.}
#' \item{pi0}{The probability mass of the null component corresponding to the estimated \code{theta}.}
#' \item{sigma2}{The variance of the non-null component corresponding to the estimated \code{theta}.}
#' \item{SE_theta}{Standard error of causal effect estimate.}
#' \item{zstat_theta}{Z-statistic for test of the causal effect estimate.}
#' \item{pvalue_theta}{P-value from the z test for the causal effect.}
#' \item{profile_mat}{A matrix of 3 columns containing details of the grid search. The first column is \code{theta_temp_vec}. The second and third columns are the corresponding \code{pi0} and \code{sigma2} values. Only returned if \code{profile=TRUE}.}
#' @references
#' Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.
#' @export
#' @examples
#' data("sumstats", package = "MRMix")
#' # Convert summary statistics to standardized scale
#' data_std = standardize(sumstats$betahat_x, sumstats$betahat_y, sumstats$sx, sumstats$sy, contbin_x = "continuous", contbin_y = "continuous", sumstats$nx, sumstats$ny, MAF = NULL)
#' # MRMix analysis
#' est = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std)
#' str(est) # True causal effect is 0.2.
#' # Include profile matrix
#' est = MRMix(data_std$betahat_x_std, data_std$betahat_y_std, data_std$sx_std, data_std$sy_std, profile = TRUE)
#' str(est)

MRMix = function(betahat_x, betahat_y, sx, sy, theta_temp_vec = seq(-1,1,by=0.01), pi_init = 0.6, sigma_init = 1e-5, profile = FALSE){
    if (length(betahat_x)!=length(betahat_y)) stop("betahat_x and betahat_y should have the same length.")
    if (length(sx)!=length(betahat_x) & length(sx)!=1) stop("sx should have the same length as betahat_x or be a single number.")
    if (length(sy)!=length(betahat_y) & length(sy)!=1) stop("sy should have the same length as betahat_y or be a single number.")

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

    ind = EM_res[,2]==max(EM_res[,2])
    est = list(theta = mean(EM_res[ind,1]),
               pi0 = max(EM_res[,2]),
               sigma2 = mean(EM_res[ind,3]))
    est$SE_theta = MRMix_se(betahat_x, betahat_y, sx, sy, est$theta, est$pi0, est$sigma2) # Standard error
    est$zstat_theta = est$theta/est$SE_theta
    est$pvalue_theta = 2*pnorm(-abs(est$zstat_theta))

    # Diagnosis
    s0_med = median(sy)^2+EM_res[,1]^2*median(sx)^2
    if (sum(EM_res[,2]>0.99 & EM_res[,3]<s0_med)>0){
        warning("Algorithm may not have converged.")
    }

    # Output results
    if (profile==TRUE){
        est$profile_mat = EM_res
    }
    return(est)
}
