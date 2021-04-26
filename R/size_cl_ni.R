

size_cl_ni <- function(mu, sigma, se_cl) {

  varcl <- se_cl^2


  # Function for numerical integration
  f_omega <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 /
          (2 * s**2)) * ((1 + exp(- x))^ (-1)) * (1 - (1 + exp(- x))^ (-1))

  # Numerical integration
  e_omega <- stats::integrate(f_omega, - Inf, Inf, mu = mu, s = sigma)$value

  n            <- 1 / varcl * (1 / e_omega)
  n            <- ceiling(n)

  return(n)
}
