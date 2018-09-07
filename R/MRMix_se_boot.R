#' Standard error of the MRMix estimator
#'
#' @description Calculate the standard error of the MRMix estimator using asymptotic theory.
#'
#' @param betahat_x GWAS effect estimates of the exposure. Vector of length \code{K}, where \code{K} is the number of instruments.
#' @param betahat_y GWAS effect estimates of the outcome. Vector of length \code{K}.
#' @param sx2 Variance of \code{betahat_x}. Vector of length \code{K}.
#' @param sy2 Variance of \code{betahat_y}. Vector of length \code{K}.
#' @param theta Estimate of causal effect. First output of MRMix
#' @param pi0 The probability mass of the null component corresponding to the estimated \code{theta}. Second output of MRMix
#' @param sigma2 The variance of the non-null component corresponding to the estimated \code{theta}. Third output of MRMix.
#'
#' @return The standard error of MRMix estimator.
#'
#' @references
#'
#' @export
MRMix_se_boot = function(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2){
    ### Main variance formula
    term1 = rep(NA,100)
    for (i in 1:100){
        ind = sample(length(betahat_x), replace=TRUE)
        betahat_x_boot = betahat_x[ind]
        betahat_y_boot = betahat_y[ind]
        if (length(sx2)==1){
            sx2_boot = sx2
        } else{
            sx2_boot = sx2[ind]
        }
        
        if (length(sy2)==1){
            sy2_boot = sy2
        } else{
            sy2_boot = sy2[ind]
        }
        
        term1[i] = sum(p2l_p2sigma2(betahat_x_boot, betahat_y_boot, sx2_boot, sy2_boot, theta, pi0, sigma2))*
            sum(p2l_ppi0_ptheta(betahat_x_boot, betahat_y_boot, sx2_boot, sy2_boot, theta, pi0, sigma2))-
            sum(p2l_ppi0_psigma2(betahat_x_boot, betahat_y_boot, sx2_boot, sy2_boot, theta, pi0, sigma2))*
            sum(p2l_psigma2_ptheta(betahat_x_boot, betahat_y_boot, sx2_boot, sy2_boot, theta, pi0, sigma2))
    }
    sqrt(1/(sum(p3l_p2sigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
           sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_ppi0_p2theta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
           sum(p3l_ppi0_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
           sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_psigma2_p2theta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)))^2*
        (var(term1)+sum((p3l_p2sigma2_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
                  sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_ppi0_ptheta_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
                  p3l_ppi0_psigma2_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
                  sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_psigma2_ptheta_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2*sx2)+
             sum((p3l_p2sigma2_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
                      sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_ppi0_ptheta_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
                      p3l_ppi0_psigma2_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
                      sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_psigma2_ptheta_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2*sy2)))
}
