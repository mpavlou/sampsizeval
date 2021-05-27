#' Precision-based Sample size Calculations for Validation of Risk Models
#' for Binary Outcomes

#' @description
#' This function calculates the sample size required in the validation dataset
#' to estimate the C-statistic (C), the calibration Slope (CS) and the
#' Calibration in the Large (CL) with sufficient precision. It takes as
#' arguments the anticipated values of the C-statistic and the outcome
#' prevalence (obtained, for example, from a previous study) and the required
#' standard error for C, CS and CL.

#' @details
#' The sample size calculations are valid under the assumption of marginal
#' normality for the distribution of the linear predictor.The default sample
#' size calculation based on C uses the closed-form expression in equation (9)
#' as proposed by *Pavlou et al. (2021)*. This is quick to run and accurate for
#' all values of anticipated C and p.The default sample size calculations based
#' on CS and CL use the formulae (12) and (13) that require the use of numerical
#' integration. The parameters of the assumed Normal
#' distribution used in the latter two expressions are obtained using equations
#' (7) and (8) and are fine-tuned for values of anticipated C>0.8.
#'
#' Sample size calculations from the estimator based on C that uses numerical
#' integration can also be obtained.
#'
#' @param p (numeric) The anticipated outcome prevalence, a real number between
#' 0 and 1
#' @param c (numeric) The anticipated C-statistic, a real number between
#' 0.5 and 1
#' @param se_c (numeric) The required standard error of the estimated
#' C-Statistic
#' @param se_cs (numeric) The required standard error of the estimated
#' Calibration Slope
#' @param se_cl (numeric) The required standard error of the estimated
#' Calibration in the Large
#' @param c_ni (logical) Numerical integration is used for the calculations for
#'             C-statistic (TRUE) or the closed-form expression (FALSE). Default
#'             value is 'FALSE'
#'
#' @return size_c: the sample size based on the C-statistic
#' @return size_cs: the sample size based on the Calibration Slope
#' @return size_cl: the sample size based on the Calibration in the Large
#' @return size_recommended: the final sample size recommendation (the largest
#' of the three above)
#' @export
#'
#' @examples
#' # Calculate the sample size of the validation data to estimate the
#' # C-statistic, the Calibration slope and the Calibration in the Large with
#' # sufficient precision. It is assumed that the anticipated prevalence is 0.1
#' # and the C-statistic is 0.75. The required SE for the C statistic is 0.025
#' # (corresponding to a confidence interval of width approximately 0.1) and the
#' # required SE for the calibration slope and calibration in the large is 0.1
#' # (corresponding to a confidence interval of width approximately 0.4).
#'
#' sampsizeval(p=0.1, c=0.75, se_c=0.025, se_cs =0.1, se_cl = 0.1)
#'
#' @references
#' Pavlou M, Chen Q, Omar ZR, Seaman RS, Steyerberg WE, White RI, Ambler G.
#' Estimation of required sample size for external validation of risk models
#' for binary outcomes, SMMR (2021). doi:10.1177/09622802211007522
#'
sampsizeval <- function(p, c, se_c, se_cs, se_cl, c_ni=FALSE) {


if (missing(p) == TRUE) stop("Please enter the anticipated prevalence")

if (missing(c) == TRUE) stop("Please enter the anticipated C-statistic")

if (missing(se_c) == TRUE) stop("Please enter required stadard error (SE) for the estimated C-statistic.
  For example, for a confidence interval of width 0.1, se_c=0.025")

if (missing(se_cs) == TRUE) stop("Please enter required stadard error (SE) for the estimated Calibration Slope.
  For example, for a confidence interval of width 0.1, se_cs=0.1")

if (missing(se_cl) == TRUE) stop("Please enter required stadard error (SE) for the estimated Calibration in the large.
  For example, for a confidence interval of width 0.1, se_cl=0.1")

if (p<=0 | p>=1) stop("Please enter a value of the anticipated outcome prevalence in (0,1)")

if (c<=0.5 | c>=1) stop("Please enter a value of the anticipated C-statistic in (0.5,1)")


if (p > 0.5) p <- 1-p

# Fine tune mu and sigma if anticipated c>0.8
if (c <= 0.8) {
  fc      <- 1.00
  sigma_c <- sqrt(2) * stats::qnorm(c) * fc
  mu      <- 0.5 * (2 * p - 1) * (sigma_c ^ 2) + log(p / (1 - p))
  sigma   <- sqrt((sigma_c^2) * (1 + p * (1 - p) * (sigma_c^2)))
} else {

  correct <- tune_mu_sigma(p, c)
  mu      <- correct$mu
  sigma   <- correct$sigma
}


if (c_ni==TRUE)  size_c  <- size_c_ni(mu = mu, sigma = sigma, se_c = se_c) else
                 size_c  <- size_c_app(p = p, c = c, se_c = se_c)
size_cs <- size_cs_ni(mu = mu, sigma = sigma, se_cs = se_cs)

size_cl <- size_cl_ni(mu = mu, sigma = sigma, se_cl = se_cl)


size_final <- max(size_c, size_cs, size_cl)

out <- list("size_c_statistic" = size_c,
            "size_calibration_slope" = size_cs,
            "size_calibration_large" = size_cl,
            "size_recommended" = size_final)

return(out)
}



