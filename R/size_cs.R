
#' Precision-based sample size calculation based on the Calibration Slope (CS)
#'
#' The purpose of this function is to receive the anticipated outcome prevalence
#' , C-statistic and calculate the sample size to achieve the required precision
#' based on CS.
#'
#'
#' @param p (numeric) The outcome prevalence, a real number between 0 and 0.5
#' @param c (numeric) The C-statistic, a real number between 0.5 and 1
#' @param varcs (numeric) The required variance of the Calibration Slope
#'
#' @return n The required sample size
#' @export
#'
#' @examples
#' size_cs(0.057, 0.77, 0.15^2)
#'
size_cs  <- function(p, c, varcs) {
a        <- 2 * p * (1 - p) * stats::qnorm(c)^2
n        <- 1 / varcs * (1 / a + 2)
n        <- ceiling(n)
events   <- ceiling(n * p)
return(n)
}
