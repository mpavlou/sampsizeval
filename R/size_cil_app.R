
size_cil_app <- function(p, c, se_cl) {
  varcl <- se_cl^2

  fc      <- 1.00
  sigma_c <- sqrt(2) * stats::qnorm(c) * fc
  mu      <- 0.5 * (2 * p - 1) * (sigma_c^2) + log(p / (1 - p))
  sigma   <- sqrt((sigma_c^2) * (1 + p * (1 - p) * (sigma_c^2)))

# CiL - Approximation (Formula 16 of main paper)

p_tilde       <-  exp(mu) / (1 + exp(mu))
term1m        <-  p_tilde * (1 - p_tilde)
term2m        <-  (1 / 2) * (1 - 6 * p_tilde + 6 * p_tilde^2) *
                  p_tilde * (1 - p_tilde) * sigma^2
n             <-  1 / (varcl * (term1m + term2m))
n             <-  ceiling(n)
events        <-  ceiling(n * p)
return(n)
}
