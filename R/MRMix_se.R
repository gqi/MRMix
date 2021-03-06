#' Standard error of the MRMix estimator
#'
#' @description This function calculates the standard error of the MRMix estimator using asymptotic theory.
#'
#' @param betahat_x GWAS effect estimates of the exposure, recommended to be in standardized scale. Vector of length \code{K}, where \code{K} is the number of instruments.
#' @param betahat_y GWAS effect estimates of the outcome, recommended to be in standardized scale. Vector of length \code{K}.
#' @param sx Standard error of \code{betahat_x}, recommended to be in standardized scale. Vector of length \code{K}.
#' @param sy Standard error of \code{betahat_y}, recommended to be in standardized scale. Vector of length \code{K}.
#' @param theta Estimate of causal effect.
#' @param pi0 The probability mass of the null component corresponding to the estimated \code{theta}.
#' @param sigma2 The variance of the non-null component corresponding to the estimated \code{theta}.
#'
#' @return The standard error of MRMix estimator.
#'
#' @references
#' Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.
#' @export
MRMix_se = function(betahat_x, betahat_y, sx, sy, theta, pi0, sigma2){
    sx2 = sx^2; sy2 = sy^2
    ### Main variance formula
    pS_ptheta = sum(p3l_p2sigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
        sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_ppi0_p2theta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
        sum(p3l_ppi0_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
        sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_psigma2_p2theta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))
    pS_psigma2 = sum(p3l_p3sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
        sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_ppi0_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
        sum(p3l_ppi0_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
        sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_p2sigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))
    psigma2_ptheta = (sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
                          sum(p2l_p2pi0(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)))/
        (sum(p2l_p2pi0(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
             sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2)
    pS_ppi0 = sum(p3l_ppi0_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
        sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_p2pi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
        sum(p3l_p2pi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
        sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p3l_ppi0_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))
    ppi0_pbetax = (sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p2l_psigma2_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
                       sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p2l_ppi0_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))/
        (sum(p2l_p2pi0(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
             sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2)
    ppi0_pbetay = (sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p2l_psigma2_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
                       sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p2l_ppi0_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))/
        (sum(p2l_p2pi0(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
             sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2)
    psigma2_pbetax = (sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p2l_ppi0_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
                       sum(p2l_p2pi0(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p2l_psigma2_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))/
        (sum(p2l_p2pi0(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
             sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2)
    psigma2_pbetay = (sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p2l_ppi0_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
                          sum(p2l_p2pi0(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p2l_psigma2_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))/
        (sum(p2l_p2pi0(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
             sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))^2)
    pS_pbetax = p3l_p2sigma2_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
        sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_ppi0_ptheta_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
        p3l_ppi0_psigma2_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
        sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_psigma2_ptheta_pbetax(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)
    pS_pbetay = p3l_p2sigma2_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_ppi0_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))+
        sum(p2l_p2sigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_ppi0_ptheta_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)-
        p3l_ppi0_psigma2_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)*sum(p2l_psigma2_ptheta(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))-
        sum(p2l_ppi0_psigma2(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2))*p3l_psigma2_ptheta_pbetay(betahat_x, betahat_y, sx2, sy2, theta, pi0, sigma2)

    sqrt(1/(pS_ptheta+pS_psigma2*psigma2_ptheta)^2*
        (sum((pS_pbetax+pS_ppi0*ppi0_pbetax+pS_psigma2*psigma2_pbetax)^2*sx2)+
             sum((pS_pbetay+pS_ppi0*ppi0_pbetay+pS_psigma2*psigma2_pbetay)^2*sy2)))
}
