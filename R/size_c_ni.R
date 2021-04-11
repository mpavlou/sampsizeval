
#' Precision-based sample size calculation based on the C-statistic (C)
#' using Numerical Integration
#'
#' The purpose of this function is to receive the anticipated outcome prevalence
#', C-statistic and calculate the sample size to achieve the required precision
#' based on the C-Statistic. It uses Numerical integration and it assumes
#' marginal normality.
#'
#'
#' @param p (numeric) The outcome prevalence, a real number between 0 and 0.5
#' @param c (numeric) The C-statistic, a real number between 0.5 and 1
#' @param varc (numeric) The required variance of the C-Statistic
#'
#' @return n         The required sample size
#' @export
#'
#' @examples
#' size_c_ni(0.057, 0.77, 0.025^2)
#'
size_c_ni <- function(p, c, varc) {
  # For the numerical integration method a distribution for the linear predictor
  # needs to be assumed. In this calculation we assume marginal normality for
  # the linear predictor

  # f1 and f2 are functions that contain the integrad in equation (5) and (6)
  # as given by Gail and Pfeiffer (2005)

  f1 <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 / (2 * s**2)) * (1  + exp(-x)) ^ (-1)

  f0 <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 / (2 * s**2)) * (1 - (1 + exp(-x)) ^ (-1))

  # The following functions are used to calculate P(eta_0<eta_1)
  # and P(eta_1<eta_0)

  intf0  <- function(upper) {
    prob <- NULL
    for (i in seq(upper)) {
      prob[i] <- stats::integrate(f0, -Inf, upper = upper[i], mu = mu,
                     s = sigma)$value /
                     stats::integrate(f0, -Inf, Inf, mu = mu, s = sigma)$value
    }
    prob
    }

  intf1 <- function(upper) {
    prob <- NULL
    for (i in seq(upper)) {
      prob[i] <- stats::integrate(f1, -Inf, upper = upper[i], mu = mu,
                      s = sigma)$value /
                      stats::integrate(f1, -Inf, Inf, mu = mu, s = sigma)$value
    }
    prob
    }

  # Input values for mu and sigma^2 assuming marginal normality
  # Factor fc can be used to fine-tuned this when C is too high (>0.85),
  # not used here

  fc      <- 1.00
  sigma_c <- sqrt(2) * stats::qnorm(c) * fc
  mu      <- 0.5 * (2 * p - 1) * (sigma_c^2) + log(p / (1 - p))
  sigma   <- sqrt((sigma_c^2) * (1 + p * (1 - p) * (sigma_c^2)))

  #We now get the cumulative distribution of eta_0 and eta1
  # A vector of values to cover the range of possible values of eta_0
  #and eta_0

  #Pre-processing to narrow down the integration range

  x <- seq(- 14, 14, 0.01)

  p_eta0 <- NULL ; p_eta1 <- NULL

  #Numerical integration to get P(eta_0<x) and P(eta_1<x)
  #Equations for Gail and Pfeiffer

  for (i in seq(x)) {
    p_eta1[i] <- stats::integrate(f1, - Inf, x[i], mu = mu, subdivisions = 1000L, s = sigma)$value /
      integrate(f1, - Inf, Inf, subdivisions = 1000L, mu = mu, s = sigma)$value
    p_eta0[i] <- stats::integrate(f0, - Inf, x[i], mu = mu, subdivisions = 1000L, s = sigma)$value /
      stats::integrate(f0, - Inf, Inf, subdivisions = 1000L, mu = mu, s = sigma)$value
  }

  d0 <- data.frame(cbind(x, p_eta0))
  d1 <- data.frame(cbind(x, p_eta1))

  d01 <- subset(d0, p_eta0 < (1 - 0.00001) & p_eta0 > 0.00001)
  d11 <- subset(d1, p_eta1 < (1 - 0.00001) & p_eta1 > 0.00001)

  a <- d01$x[1]
  b <- d11$x[nrow(d11)]

  ## Actual numerical integration starts here

  step <- 0.0005

  x <- seq(a, b, step)

  p_eta0 <- NULL
  p_eta1 <- NULL

  #Numerical integration to get P(eta_0<x) and P(eta_1<x)
  #Equations for Gail and Pfeiffer (2005)

  for (i in seq(x)) {
    p_eta1[i] <- stats::integrate(f1, -Inf, x[i], mu = mu, s = sigma)$value /
      stats::integrate(f1, - Inf, Inf, mu = mu, s = sigma)$value
    p_eta0[i] <- stats::integrate(f0, -Inf, x[i], mu = mu, s = sigma)$value /
      stats::integrate(f0, - Inf, Inf, mu = mu, s = sigma)$value
  }

  p_eta0[1] <- 0 ; p_eta1[1] <- 0
  p_eta0[length(x)] <- 1 ; p_eta1[length(x)] <- 1


  # Inverse CDF to sample from the distribution of eta_0 and eta_1
  nsamps   <- 400000
  u        <- sort(stats::runif(nsamps))
  eta0     <- pracma::interp1(x = sort(p_eta0), y = x, xi = u)

  u    <- sort(stats::runif(nsamps))
  eta1 <- pracma::interp1(x = sort(p_eta1), y = x, xi = u)


  #Note: Distribution of eta_0 and eta_1 is approximately Normal, as expected


  # Compute EK=E(P(eta_0<eta_1))
  prob  <-  NULL
  prob  <-  intf0(eta1)
  e_k2  <-  mean(prob^2)

  #Calculate  the true C and p for these values of mu and sigma
  c_ni      <- mean(prob)
  p_ni      <- stats::integrate(f1, - Inf, Inf, mu = mu, s = sigma)$value
  c_true    <- c_ni
  p_true    <- p_ni

  # Compute EG=E(P(eta_1<eta_0))
  prob <- NULL
  prob <- intf1(eta0)
  e_g2 <- mean((1 - prob)^2)

  #Enter in the final formula
  n         <- (1 / varc) * ((1 - p_true) * e_k2 + p_true * e_g2 - c_true^2) /
               (p_true * (1 - p_true))
  n         <- ceiling(n)
  events    <- ceiling(n * p)

  return(n)
}
