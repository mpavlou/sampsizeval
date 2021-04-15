
size_cs_ni <- function(p, c, se_cs) {

  varcs <- se_cs^2

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

  #Functions for numerical integration

f_om <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 /
            (2 * s**2)) * ((1 + exp(- x))^ (-1)) * (1 - (1 + exp(- x))^ (-1))

f_om_eta <- function(x, mu, s) 1 / (2 * s**2 * pi)^0.5 * exp(- (x - mu)^2 /
        (2 * s**2)) * ((1 + exp(- x))^ (- 1)) * (1 - (1 + exp(- x))^ (- 1)) * x

f_om_etasq <- function(x, mu, s) 1 / (2 * s**2 * pi)^0.5 * exp(- (x  - mu)^2 /
      (2 * s**2)) * ((1 + exp(- x))^ (- 1)) * (1 - (1 + exp(- x))^ (- 1)) * x^2

  #Numerical integration

e_om        <- stats::integrate(f_om, -Inf, Inf, mu = mu, s = sigma)$value
e_om_eta_sq <- stats::integrate(f_om_etasq, -Inf, Inf, mu = mu, s = sigma)$value
e_om_eta    <- stats::integrate(f_om_eta, -Inf, Inf, mu = mu, s = sigma)$value

  numer        <- e_om
  denom        <- e_om * e_om_eta_sq - e_om_eta ^ 2
  n            <- 1 / varcs * (numer / denom)
  n            <- ceiling(n)
  events       <- ceiling(n * p)

  return(n)
}
