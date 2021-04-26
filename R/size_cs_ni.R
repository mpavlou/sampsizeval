
size_cs_ni <- function(mu, sigma, se_cs) {

  varcs <- se_cs^2

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

  return(n)
}
