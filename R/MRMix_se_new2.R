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
MRMix_se = function(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2){
    ### Main variance formula
    sqrt(1/(sum(p3l_p2sigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
           sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_ppi0_p2theta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
           sum(p3l_ppi0_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
           sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_psigma2_p2theta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
               (sum(p3l_p3sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
               sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_ppi0_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
               sum(p3l_ppi0_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
               sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_p2sigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)))*
               (sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
                    sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)))/
               (sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2))^2*
        (sum((p3l_p2sigma2_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
                  sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_ppi0_ptheta_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
                  p3l_ppi0_psigma2_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
                  sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_psigma2_ptheta_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2*sx2)+
             sum((p3l_p2sigma2_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
                      sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_ppi0_ptheta_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
                      p3l_ppi0_psigma2_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
                      sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_psigma2_ptheta_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2*sy2)))
}
