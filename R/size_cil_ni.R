

size_cil_ni <- function(p, c, se_cl) {

  varcl <- se_cl^2

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
