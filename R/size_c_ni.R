
#' Precision-based sample size calculation based on the C-statistic (C) -
#' Numerical Integration
#'
#' The purpose of this function is to receive the anticipated outcome prevalence
#', C-statistic and calculate the sample size to achieve the required precision
#' based on the C-Statistic. It uses Numerical integration and it assumes
#' marginal normality.
#'
#'
#' @param p (numeric) The outcome prevalence, a real number between 0 and 0.5
#' @param c (numeric) The C-statistic, a real number between 0.5 and 1
#' @param se_c (numeric) The required standard error variance of the estimated
#' C-Statistic
#'
#' @return n         The required sample size

#'
size_c_ni <- function(p, c, se_c) {

  varc <- se_c^2

  # Input values for mu and sigma^2 assuming marginal normality
  # factor fc can be used to fine-tuned this when C is too high (>0.8),
  # not used here


  if (c <= 0.8) {
    fc      <- 1.00
    sigma_c <- sqrt(2) * stats::qnorm(c) * fc
    mu      <- 0.5 * (2 * p - 1) * (sigma_c^2) + log(p / (1 - p))
    sigma   <- sqrt((sigma_c^2) * (1 + p * (1 - p) * (sigma_c^2)))
  } else {

    correct <- tune_mu_sigma(p, c)
    mu      <- correct$mu
    sigma   <- correct$sigma
  }

  # For the numerical integration method a distribution for the linear predictor
  # needs to be assumed. In this calculation we assume marginal normality for
  # the linear predictor

  # f1 and f2 are functions that contain the integrad in equation (5) and (6)
  # as given by Gail and Pfeiffer (2005)


  f1 <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 /
                               (2 * s**2)) * (1  + exp(-x)) ^ (-1)

  f0 <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 /
                               (2 * s**2)) * (1 - (1 + exp(-x)) ^ (-1))

  # The following functions are used to calculate P(eta_0<eta_1)
  # and P(eta_1<eta_0)

  intf0   <- function(upper) {
    prob  <- NULL
    denom <- stats::integrate(f0, -Inf, Inf, mu = mu, s = sigma)$value
    for (i in seq(upper)) {
      prob[i] <- stats::integrate(f0, -Inf, upper = upper[i], mu = mu,
                     s = sigma)$value / denom

    }
    prob
    }

  intf1  <- function(upper) {
    prob  <- NULL
    denom <-stats::integrate(f1, -Inf, Inf, mu = mu, s = sigma)$value
    for (i in seq(upper)) {
      prob[i] <- stats::integrate(f1, -Inf, upper = upper[i], mu = mu,
                      s = sigma)$value / denom

    }
    prob
    }

  #We now get the cumulative distribution of eta_0 and eta1
  # A vector of values to cover the range of possible values of eta_0
  #and eta_0

  #Pre-processing to narrow down the integration range

  x <- seq(- 14, 14, 0.01)

  p_eta0 <- NULL ; p_eta1 <- NULL

  #Numerical integration to get P(eta_0<x) and P(eta_1<x)
  #Equations for Gail and Pfeiffer
  denom1 <- stats::integrate(f1, - Inf, Inf, mu = mu, s = sigma,
                             subdivisions=10000L)$value
  denom0 <- stats::integrate(f0, - Inf, Inf, mu = mu, s = sigma,
                             subdivisions=1000L)$value

  for (i in seq(x)) {
    p_eta1[i] <- stats::integrate(f1, - Inf, x[i], mu = mu, s = sigma)$value /
             denom1
    p_eta0[i] <- stats::integrate(f0, - Inf, x[i], mu = mu, s = sigma)$value /
             denom0
  }

  d0 <- data.frame(cbind(x, p_eta0))
  d1 <- data.frame(cbind(x, p_eta1))

  d01 <- subset(d0, p_eta0 < (1 - 0.00001) & p_eta0 > 0.00001)
  d11 <- subset(d1, p_eta1 < (1 - 0.00001) & p_eta1 > 0.00001)

  a <- d01$x[1]
  b <- d11$x[nrow(d11)]

  ## Actual numerical integration starts here

  step <- 0.00005

  x <- seq(a, b, step)

  p_eta0 <- NULL
  p_eta1 <- NULL

  #Numerical integration to get P(eta_0<x) and P(eta_1<x)
  #Equations for Gail and Pfeiffer (2005)

  for (i in seq(x)) {
    p_eta1[i] <- stats::integrate(f1, -Inf, x[i], mu = mu, s = sigma,
                                  subdivisions=1000L)$value / denom1
    p_eta0[i] <- stats::integrate(f0, -Inf, x[i], mu = mu, s = sigma,
                                  subdivisions=1000L)$value / denom0
  }

  p_eta0[1] <- 0
  p_eta1[1] <- 0
  if (p_eta0[length(x)] !=1 ) p_eta0[length(x)] <- 1
  if (p_eta1[length(x)] !=1 ) p_eta1[length(x)] <- 1


  #mydata <- data.frame(x,p_eta0)
  #mydata[duplicated(mydata1$p_eta1), ]
  #mydata0 <- data.frame(x,p_eta0)

  #mydata1 <- mydata[duplicated(mydata1$p_eta1), ]
  #mydata0 <- mydata0[!duplicated(mydata0$p_eta0), ]

  # Inverse CDF to sample from the distribution of eta_0 and eta_1
  #nsamps   <- length(x)
  nsamps   <- 1000000
  u        <- sort(stats::runif(nsamps))
  eta0     <- pracma::interp1(x = sort(p_eta0), y = x, xi = u)

  u        <- sort(stats::runif(nsamps))
  eta1     <- pracma::interp1(x = sort(p_eta1), y = x, xi = u)


  #Note: Distribution of eta_0 and eta_1 is approximately Normal, as expected


  # Compute EK=E(P(eta_0<eta_1))
  prob  <-  NULL
  prob  <-  intf0(eta1)
  e_k2  <-  mean(prob^2)

  #Calculate  the true C and p for these values of mu and sigma
  c_ni      <- mean(prob)
  p_ni      <- denom1
  c_true    <- c_ni
  p_true    <- p_ni

  # Compute EG=E(P(eta_1<eta_0))
  prob <- NULL
  prob <- intf1(eta0)
  e_g2 <- mean((1 - prob)^2)
  #e_g2 <- e_k2

  #Enter in the final formula
  n         <- (1 / varc) * ((1 - p_true) * e_k2 + p_true * e_g2 - c_true^2) /
               (p_true * (1 - p_true))
  n         <- ceiling(n)
  events    <- ceiling(n * p)

  return(n)
}
