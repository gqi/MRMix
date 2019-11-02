#' Standardize summary statistics for MRMix analysis
#'
#' @description This function standardizes GWAS summary statistics (beta and standard error (SE)), i.e. 1) for continuous phenotype, transform them to the beta and SE when both the genotypes and phenotypes are standardized to have variance 1; 2) for binary phenotype, transform them to the beta and SE when the genotypes are standardized to have variance 1. This function applies to continuous phenotypes analyzed with linear regression and binary phenotypes analyzed with logistic regression. Do not use for other types of phenotypes or other analytic methods.
#'
#' @param betahat_x GWAS effect estimates of the exposure. Vector of length \code{K}, where \code{K} is the number of instruments.
#' @param betahat_y GWAS effect estimates of the outcome. Vector of length \code{K}.
#' @param sx Standard error of \code{betahat_x}. Vector of length \code{K}.
#' @param sy Standard error of \code{betahat_y}. Vector of length \code{K}.
#' @param contbin_x Is the exposure a continuous or binary trait? Set to \code{contbin_x="continuous"} or \code{contbin_x="binary"}. Or set to \code{contbin_x="n"} if standardization is not needed.
#' @param contbin_y Is the outcome a continuous or binary trait? Set to \code{contbin_y="continuous"} or \code{contbin_y="binary"}. Or set to \code{contbin_y="n"} if standardization is not needed.
#' @param nx SNP-specific sample size of the study associated with the exposure. Vector of length \code{K} or a single number. A single number indicates all the SNPs have the same sample size. Set to \code{NULL} if trait is binary.
#' @param ny SNP-specific sample size of the study associated with the outcome. Vector of length \code{K} or a single number. A single number indicates all the SNPs have the same sample size. Set to \code{NULL} if trait is binary.
#' @param MAF Minor allele frequency. Vector of length \code{K}. Set to \code{NULL} if both traits are continuous.
#'
#' @details For continuous phenotypes analyzed with linear regression, data are standardized by \code{beta_standardized=beta/(se*sqrt(N)); se_standardized=1/sqrt(N)}; for binary phenotypes analyzed with logistic regression, data are standardized by \code{beta_standardized=beta*sqrt(2*MAF*(1-MAF)); se_standardized=se*sqrt(2*MAF*(1-MAF))}.
#'
#' @return A list that contains
#' \item{betahat_x_std}{Standardized \code{betahat} for the exposure.}
#' \item{betahat_y_std}{Standardized \code{betahat} for the outcome.}
#' \item{sx_std}{Standard error of \code{betahat_x}.}
#' \item{sy_std}{Standard error of \code{betahat_y}.}
#'
#' @references
#' Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.
#' @export
standardize = function(betahat_x, betahat_y, sx, sy, contbin_x, contbin_y, nx, ny, MAF){
    if (length(betahat_x)!=length(betahat_y)) stop("betahat_x and betahat_y should have the same length.")
    if (length(sx)!=length(betahat_x) & length(sx)!=1) stop("sx should have the same length as betahat_x or be a single number.")
    if (length(sy)!=length(betahat_y) & length(sy)!=1) stop("sy should have the same length as betahat_y or be a single number.")
    if (length(nx)!=length(betahat_x) & length(nx)!=1 & !is.null(nx)) stop("nx should have the same length as betahat_x or be a single number.")
    if (length(ny)!=length(betahat_y) & length(ny)!=1 & !is.null(nx)) stop("ny should have the same length as betahat_y or be a single number.")
    if (length(MAF)!=length(betahat_y) & !is.null(MAF)) stop("MAF should have the same length as betahat_y or be NULL.")


    if (contbin_x=="continuous"){
        betahat_x_std = betahat_x/sx/sqrt(nx)
        sx_std = 1/sqrt(nx)
        if (length(sx_std)==1) sx_std = rep(sx_std,length(betahat_x_std))
    } else if (contbin_x=="binary"){
        betahat_x_std = betahat_x*sqrt(2*MAF*(1-MAF))
        sx_std = sx*sqrt(2*MAF*(1-MAF))
    } else if (contbin_x=="n"){
        betahat_x_std = betahat_x
        sx_std = sx
    } else{
        stop("contbin_x out of range.")
    }

    if (contbin_y=="continuous"){
        betahat_y_std = betahat_y/sy/sqrt(ny)
        sy_std = 1/sqrt(ny)
        if (length(sy_std)==1) sy_std = rep(sy_std,length(betahat_y_std))
    } else if (contbin_y=="binary"){
        betahat_y_std = betahat_y*sqrt(2*MAF*(1-MAF))
        sy_std = sy*sqrt(2*MAF*(1-MAF))
    } else if (contbin_y=="n"){
        betahat_y_std = betahat_y
        sy_std = sy
    } else{
        stop("contbin_y out of range.")
    }

    return(list(betahat_x_std = betahat_x_std, betahat_y_std = betahat_y_std,
                sx_std = sx_std, sy_std = sy_std))
}
