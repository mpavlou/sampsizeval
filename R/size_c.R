
#' Precision-based sample size calculation based on the C-statistic (C)
#'
#' The purpose of this function is to receive the anticipated outcome prevalence
#' , C-statistic and calculate the sample size to achieve the required precision
#' based on the C-statistic.
#'
#'
#' @param p (numeric) The outcome prevalence, a real number between 0 and 0.5
#' @param c (numeric) The C-statistic,  a real number between 0.5 and 1
#' @param varc (numeric) The required variance of the C- statistic
#'
#' @return n The required sample size
#' @export
#'
#' @examples
#' size_c(0.057, 0.77, 0.025^2)
#'
size_c <- function(p, c, varc) {
  n            <- (c - 2 * sn::T.Owen(-stats::qnorm(c), 1 / sqrt(3)) - c^2) /
                  ((p - p^2) * varc)
  n            <- ceiling(n)
  events       <- ceiling(n * p)
  return(n)
}
