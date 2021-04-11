
#' Precision-based sample size calculation based on the Calibration in the Large
#' (CiL) using Numerical Integration
#'
#' The purpose of this function is to receive the anticipated outcome prevalence
#' based on CiL. It uses Numerical integration and it assumes marginal normality
#'
#'
#'
#' @param p (numeric) The outcome prevalence, a real number between 0 and 0.5
#' @param c (numeric) The C-statistic a real number between 0.5 and 1
#' @param varcl (numeric) The required variance of the Calibration in the Large
#'
#' @return n The required sample size
#' @export
#'
#' @examples
#' size_cil_ni(0.57, 0.77, 0.15^2)

size_cil_ni <- function(p, c, varcl) {

  fc      <- 1.00
  sigma_c <- sqrt(2) * stats::qnorm(c) * fc
  mu      <- 0.5 * (2 * p - 1) * (sigma_c^2) + log(p / (1 - p))
  sigma   <- sqrt((sigma_c^2) * (1 + p * (1 - p) * (sigma_c^2)))

  # Function for numerical integration
  f_omega <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 /
          (2 * s**2)) * ((1 + exp(- x))^ (-1)) * (1 - (1 + exp(- x))^ (-1))

  # Numerical integration
  e_omega <- stats::integrate(f_omega, - Inf, Inf, mu = mu, s = sigma)$value

  n            <- 1 / varcl * (1 / e_omega)
  n            <- ceiling(n)
  events       <- ceiling(n * p)

  return(n)
}
