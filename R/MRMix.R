#' Mendelian randomization analysis using mixture models
#'
#' @description Conducts Mendelian randomization analysis by fitting a mixture normal model for the bivariate effect size distribution of the exposure and the outcome.
#'
#' @param betahat_x GWAS effect estimates of the exposure. Vector of length \code{K}, where \code{K} is the number of instruments.
#' @param betahat_y GWAS effect estimates of the outcome. Vector of length \code{K}.
#' @param sx2 Variance of \code{betahat_x}. Vector of length \code{K}.
#' @param sy2 Variance of \code{betahat_y}. Vector of length \code{K}.
#' @param theta_temp_vec A vector of the grid search values of the causal effect \code{theta}.
#' @param pi_init Initial value of the probability mass at the null component. See details.
#' @param sigma_init Initial value of the variance of the non-null component. See details.
#' @param profile Logical. Default to be \code{FALSE} and only a list of 3 containing the estimate of causal effect (\code{theta}) and corresponding \code{pi0} and \code{sigma2}. If TRUE, returns a matrix of 3 columns. The first column is \code{theta_temp_vec}. The second and third column are the corresponding \code{pi0} and \code{sigma2} values.

#' @details The algorithm searches over a grid of possible values of the causal effect $\theta$. For each fixed $\theta$, fit mixture model \code{pi0N(0,sy2+theta^2 sx2)+(1-\pi0)N(0,sigma2)} on the residual \code{betahat_y-theta*betahat_x}. Choose the value of \code{theta} that have the maximum \code{pi0} as the estimate of causal effect.
#'
#' @return A list that contains
#' \item{theta}{Estimate of causal effect.}
#' \item{pi0}{The probability mass of the null component corresponding to the estimated \code{theta}.}
#' \item{sigma2}{The variance of the non-null component corresponding to the estimated \code{theta}.}
#'
#' @references
#'
#' @export
MRMix = function(betahat_x, betahat_y, sx2, sy2, theta_temp_vec = seq(-0.49,0.5,by=0.01), pi_init = 0.6, sigma_init = 1e-5, profile = FALSE){
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
        for (iter in 1:50000){
            f1 = 1/sqrt(2*pi*sigma0)*exp(-0.5/sigma0*(betahat_y-theta_temp*betahat_x)^2)
            f1[f1<1e-300] = 1e-300
            loglkl = sum(log(pi0*f0+(1-pi0)*f1))
            # if (iter%%1000==0) print(c(iter,pi0,sigma0,loglkl))
            # print(c(iter, pi0,sigma0,loglkl))
            pt = pi0*f0/(pi0*f0+(1-pi0)*f1)
            pt[pt>0.9999999] = 0.9999999
            pi0 = mean(pt)
            sigma0 = sum((1-pt)*(betahat_y-theta_temp*betahat_x)^2)/(length(betahat_x)-sum(pt))

            # if (pi0<1e-10) pi0 = 1e-10
            # if (pi0>0.9999999) pi0 = 0.9999999
            # if (sigma0>0.01) sigma0=0.01
            
            if (pi0<0.0001) pi0 = 0.0001
            if (pi0>0.9999) pi0 = 0.9999
            if (sigma0>0.01) sigma0=0.01
            if (sigma0<1e-7) sigma0=1e-7

            if (iter%%1000==0){
                # print(paste("theta_temp", theta_temp, "pi0", pi0, "loglkl", loglkl))
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
        # print(paste("theta_temp", theta_temp, "pi0", pi0, "loglkl", loglkl))
    }

    if (profile==TRUE){
        return(EM_res)
    } else{
        ind = EM_res[,2]==max(EM_res[,2])
        est = list(theta = mean(EM_res[ind,1]),
                   pi0 = max(EM_res[,2]),
                   sigma2 = mean(EM_res[ind,3]))
        return(est)
    }
}
