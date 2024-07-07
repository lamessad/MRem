#' Mendelian Randomization analysis
#'
#' @description This function conducts Mendelian randomization analysis based on latent phenotype of the outcome that explicitly exclude vertical pleiotropy effects and iteratively refines causal relationship  through an Expectation and Maximisation (EM) algorithm.MRvpi takes GWAS summary statistics as inputs to estimate causal effects of one trait on another. We recommend using summary statistics in the standardized scale
#' @param betaY beta of outcome, recommended to be in standardized scale.
#' @param betaX beta of exposure, recommended to be in standardized scale.
#' @param betaYse standard error of outcome, recommended to be in standardized scale.
#' @param betaXse standard error of exposure, recommended to be in standardized scale.
#' @param ny sample size of outcome gwas.
#' @param gwas_p p value of association between variants and the latent phenotype (horizontal pleiotropy test)
#' @param gwas_p2 p value of exposure gwas (SignifThreshold).
#' @param permutn number of permutations.
#' @details None.
#' @keywords Mendelian randomization.
#' @importFrom stats cor ecdf lm pchisq pt quantile
#' @importFrom utils capture.output
#' @export
#' @return A list that contains
#' \item{CausEst}{Estimate of causal effect.}
#' \item{CausEstSE}{Standard error of causal effect estimate.}
#' \item{CausEstP}{P-value from the z test for the causal effect.}
#' \item{SNPP}{P-value from the for the direct causal effect of variant on the the outcome.}
#' \item{VPI}{index of valid variants.}
mr_em <- function(betaY, betaX, betaYse, betaXse, ny, gwas_p = 5e-2, gwas_p2 = 5e-8, permutn = 0) {

  mr_vpi <- function(betaY, betaX, betaYse, betaXse, ny, gwas_p, gwas_p2) {
    sv1 <- seq_along(betaY)
    svv <- integer(0)
    ivw <- IVW(betaY[svv], betaX[svv], betaYse[svv])

    yi <- 0
    tau_t <- -99999
    while (tau_t != ivw[1] & yi <= 10) {
      yi <- yi + 1
      tau_t <- ivw[1]
      gwas5 <- matrix(0, nrow = length(betaY), ncol = 4)
      pcor <- cor(betaY, betaX)
      gwas5[, 1] <- betaY - tau_t * betaX
      gwas5[, 2] <- sqrt(betaYse^2 + (tau_t^2 / ny) - (2 * tau_t * pcor / ny))
      gwas5[, 3] <- gwas5[, 1] / gwas5[, 2]
      gwas5[, 4] <- pchisq(gwas5[, 3]^2, df = 1, lower.tail = FALSE)

      a3_p <- pchisq((betaX / betaXse)^2, df = 1, lower.tail = FALSE)
      svv <- sv1[gwas5[, 4] > gwas_p & a3_p < gwas_p2]

      ivw <- IVW(betaY[svv], betaX[svv], betaYse[svv])
    }

    ivw <- IVW(betaY[svv], betaX[svv], betaYse[svv])
    result <- list(CausEst = ivw[1], CausEstSE = ivw[2], CausEstP = ivw[4], SNPP = gwas5[, 4], VPI = svv)
    return(result)
  }


  MRem_result <- mr_vpi(betaY, betaX, betaYse, betaXse, ny, gwas_p, gwas_p2)

  if (permutn > 0) {

    permutp <- numeric(permutn)

    for (zi in seq_len(permutn)) {
      sv2 <- sample(seq_along(betaX))
      vpi_result <- mr_vpi(betaY[sv2], betaX, betaYse[sv2], betaXse, ny, gwas_p, gwas_p2)
      permutp[zi] <- vpi_result$CausEstP
    }

    permutt <- quantile(permutp, probs = 0.05, na.rm = TRUE)
    permutt2 <- ecdf(permutp)(MRem_result$CausEstP)
    permut_result <- list(sig_v = permutt, corrected_p = permutt2)


    result <- list(MRem = MRem_result, permutP = permut_result)
  } else {

    result <- MRem_result
  }

  return(result)
}
