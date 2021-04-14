#' Precision-based sample size calculations for the validation of risk models
#' for binary outcomes
#'
#' This function calculates the sample size required in the validation dataset
#' to estimate the C-statistic (C), the calibration Slope (CS) and the
#' Calibration in the Large (CL) with sufficient precision. It takes as
#' arguments the anticipated values of the C-statistic and the outcome
#' prevalence(obtained, for example, from a previous study) and the required
#' standard error for C, CS and CL.
#'
#' The sample size calculations are valid under the assumption of marginal
#' normality for the distribution of the linear predictor.The default sample
#' size calculation based on C uses the closed-form expression in equation (9).
#' This is quick to run and accurate for all values of C and p.
#' The default sample size calculation based on CS and CL
#' use the formulae (12) and (13) that require the use of numerical integration.
#' These are also quick to run.
#'
#' Results from the estimator based on C that uses numerical integration  and
#' the closed-form estimators based on CS and CL that make the use of
#' approximations can be obtained from functions  'size_c_ni', 'size_cs_app'
#' and 'size_cl_app'.
#'
#' @param p (numeric) The outcome prevalence, a real number between 0 and 0.5
#' @param c (numeric) The C-statistic, a real number between 0.5 and 1
#' @param se_c (numeric) The required standard error of the estimated
#' C-Statistic
#' @param se_cs (numeric) The required standard error of the estimated
#' Calibration Slope
#' @param se_cl (numeric) The required standard error of the estimated
#' Calibration in the large
#'
#' @return size_c: the sample size based on the C-statistic,
#' @return size_cs: the sample size based on the Calibration Slope,
#' @return size_cl: the sample size based on the Calibration in the Large and the final recommended sample size
#' @return size_recommended: the final sample size recommendation
#' @export
#'
#' @examples
#' size_val(p=0.1, c=0.75, se_c=0.025, se_cs =0.1, se_cl = 0.1)
size_val <- function(p, c, se_c = 0.025, se_cs = 0.1, se_cl = 0.1) {

size_c  <- size_c_app(p = p, c = c, se_c = se_c)
size_cs <- size_cs_ni(p = p, c = c, se_cs = se_cs)
size_cl <- size_cil_ni(p = p, c= c, se_cl = se_cl)

size_final <- max(size_c, size_cs, size_cl)

out <- list("size_c_statistic" = size_c,
            "size_calibration_slope" = size_cs,
            "size_calibration_large" = size_cl,
            "size_recommended" = size_final)

return(out)
  }
