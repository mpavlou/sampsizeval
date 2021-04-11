
#' Precision-based sample size calculation based on the Calibration Slope (CS)
#' using Numerical Integration
#'
#' The purpose of this function is to receive the anticipated outcome prevalence
#' , -statistic and calculate the sample size to achieve the required precision
#' based on CS. It uses Numerical integration and it assumes marginal normality.
#'
#'
#' @param p (numeric) The outcome prevalence, a real number between 0 and 0.5
#' @param c (numeric) The C-statistic, a real number between 0.5 and 1
#' @param varcs (numeric) The required variance of the Calibration slope
#'
#' @return n The required sample size
#' @export
#'
#' @examples
#' size_cs_ni(0.057, 0.77, 0.15^2)
#'
size_cs_ni <- function(p, c, varcs) {

  fc      <- 1.00
  sigma_c <- sqrt(2) * stats::qnorm(c) * fc
  mu      <- 0.5 * (2 * p - 1) * (sigma_c^2) + log(p / (1 - p))
  sigma   <- sqrt((sigma_c^2) * (1 + p * (1 - p) * (sigma_c^2)))


  #Functions for numerical integration

  f_omega <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 /
            (2 * s**2)) * ((1 + exp(- x))^ (-1)) * (1 - (1 + exp(- x))^ (-1))

  f_omega_eta   <- function(x, mu, s) 1 / (2 * s**2 * pi)^0.5 * exp(- (x - mu)^2 /
               (2 * s**2)) * ((1 + exp(- x))^ (- 1)) * (1 - (1 + exp(- x))^ (- 1)) * x

  f_omega_etasq <- function(x, mu, s) 1 / (2 * s**2 * pi)^0.5 * exp(- (x  - mu)^2 /
               (2 * s**2)) * ((1 + exp(- x))^ (- 1)) * (1 - (1 + exp(- x))^ (- 1)) * x^2

  #Numerical integration

  e_omega        <- stats::integrate(f_omega, - Inf, Inf, mu = mu, s = sigma)$value
  e_omega_eta_sq <- stats::integrate(f_omega_etasq, - Inf, Inf, mu = mu, s = sigma)$value
  e_omega_eta    <- stats::integrate(f_omega_eta, - Inf, Inf, mu = mu, s = sigma)$value

  numer        <- e_omega
  denom        <- e_omega * e_omega_eta_sq - e_omega_eta ^ 2
  n            <- 1 / varcs * (numer / denom)
  n            <- ceiling(n)
  events       <- ceiling(n * p)

  return(n)
}
