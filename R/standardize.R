#' Standardize summary statistics for MRMix analysis
#'
#' @description 1) For both binary and continuous traits, this function standardizes GWAS summary statistics by genotypic variance; 2) In addition, for continuous phenotype, this function standardizes summary statistics by phenotypic variance. This function is designed for GWAS estimates from linear or logistic regression. Do not use for other models.
#'
#' @param betahat_x GWAS effect estimates of the exposure. Vector of length \code{K}, where \code{K} is the number of instruments (SNPs).
#' @param betahat_y GWAS effect estimates of the outcome. Vector of length \code{K}.
#' @param sx Standard error of \code{betahat_x}. Vector of length \code{K}.
#' @param sy Standard error of \code{betahat_y}. Vector of length \code{K}.
#' @param xtype Is the exposure a continuous or binary trait? Set to \code{xtype="continuous"} or \code{xtype="binary"}. Or set to \code{xtype="n"} if exposure summary statistics do not need to be standardized.
#' @param ytype Is the outcome a continuous or binary trait? Set to \code{ytype="continuous"} or \code{ytype="binary"}. Or set to \code{ytype="n"} if outcome summary statistics do not need to be standardized.
#' @param nx SNP-specific sample size (recommended) or total sample size of the study associated with the exposure. Vector of length \code{K} or a single number. Set to \code{NULL} if trait is binary. Summary statistics for binary traits are standardized by the genotypic variance which can be calculated using the minor allele frequency (MAF) under Hardy-Weinberg equilibrium. Hence sample size is not needed for binary traits.
#' @param ny SNP-specific sample size (recommended) or total sample size of the study associated with the outcome. Vector of length \code{K} or a single number. Set to \code{NULL} if trait is binary for the same reason as for \code{nx}.
#' @param MAF Minor allele frequency. Vector of length \code{K}. Set to \code{NULL} if both traits are continuous. Summary statistics for continuous traits are standardized as z statistics rescaled by sample size, hence MAF is not needed.
#'
#' @return A list that contains
#' \item{betahat_x_std}{Standardized \code{betahat} for the exposure.}
#' \item{betahat_y_std}{Standardized \code{betahat} for the outcome.}
#' \item{sx_std}{Standard error of \code{betahat_x_std}.}
#' \item{sy_std}{Standard error of \code{betahat_y_std}.}
#'
#' @details Using the exposure \code{X} as an example: 1) For continuous phenotypes analyzed with linear regression, data are standardized by \code{betahat_x_std=betahat_x/(sx*sqrt(nx)); sx_std=1/sqrt(nx)}. Note that this standardization assumes that GWAS was conducted without covariate adjustment or the covariates do not have strong effects on \code{Y}. If the covariates have strong effects on \code{Y}, set \code{nx} equal to the effective sample size, which can be approximated by \code{N/(1-R2)}, where \code{N} is the sample size associated with the study for \code{X} and \code{R2} is the R-squared for the covariates. 2) For binary phenotypes analyzed with logistic regression, data are standardized by \code{betahat_x_std=betahat_x*sqrt(2*MAF*(1-MAF)); sx_std=sx*sqrt(2*MAF*(1-MAF))}. Same formulas apply to the outcome Y.
#'
#' @references
#' 1. Qi, Guanghao, and Nilanjan Chatterjee. "Mendelian randomization analysis using mixture models for robust and efficient estimation of causal effects." Nature Communications 10.1 (2019): 1941.
#'
#' 2. Qi, Guanghao, and Nilanjan Chatterjee. "A Comprehensive Evaluation of Methods for Mendelian Randomization Using Realistic Simulations of Genome-wide Association Studies." bioRxiv (2019): 702787.
#' @export
standardize = function(betahat_x, betahat_y, sx, sy, xtype, ytype, nx, ny, MAF){
    if (length(betahat_x)!=length(betahat_y)) stop("betahat_x and betahat_y should have the same length.")
    if (length(sx)!=length(betahat_x) & length(sx)!=1) stop("sx should have the same length as betahat_x or be a single number.")
    if (length(sy)!=length(betahat_y) & length(sy)!=1) stop("sy should have the same length as betahat_y or be a single number.")
    if (length(nx)!=length(betahat_x) & length(nx)!=1 & !is.null(nx)) stop("nx should have the same length as betahat_x or be a single number.")
    if (length(ny)!=length(betahat_y) & length(ny)!=1 & !is.null(ny)) stop("ny should have the same length as betahat_y or be a single number.")
    if (length(MAF)!=length(betahat_y) & !is.null(MAF)) stop("MAF should have the same length as betahat_y or be NULL.")


    if (xtype=="continuous"){
        betahat_x_std = betahat_x/sx/sqrt(nx)
        sx_std = 1/sqrt(nx)
        if (length(sx_std)==1) sx_std = rep(sx_std,length(betahat_x_std))
    } else if (xtype=="binary"){
        betahat_x_std = betahat_x*sqrt(2*MAF*(1-MAF))
        sx_std = sx*sqrt(2*MAF*(1-MAF))
    } else if (xtype=="n"){
        betahat_x_std = betahat_x
        sx_std = sx
    } else{
        stop("xtype out of range.")
    }

    if (ytype=="continuous"){
        betahat_y_std = betahat_y/sy/sqrt(ny)
        sy_std = 1/sqrt(ny)
        if (length(sy_std)==1) sy_std = rep(sy_std,length(betahat_y_std))
    } else if (ytype=="binary"){
        betahat_y_std = betahat_y*sqrt(2*MAF*(1-MAF))
        sy_std = sy*sqrt(2*MAF*(1-MAF))
    } else if (ytype=="n"){
        betahat_y_std = betahat_y
        sy_std = sy
    } else{
        stop("ytype out of range.")
    }

    return(list(betahat_x_std = betahat_x_std, betahat_y_std = betahat_y_std,
                sx_std = sx_std, sy_std = sy_std))
}
