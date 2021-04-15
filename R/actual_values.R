
actual_values <- function(p, c, fc = 1) {
  nevents <- 300000
  n       <- nevents / p

  sigma_c <- sqrt(2) * stats::qnorm(c) * fc
  mu      <- 0.5 * (2 * p - 1) * (sigma_c^2) + log(p / (1 - p))
  sigma   <- sqrt((sigma_c^2) * (1 + p * (1 - p) * (sigma_c^2)))

  eta   <- stats::rnorm(n, mu, sigma)
  p_est <- (1 + exp(- eta)) ^ (-1)
  y     <- stats::rbinom(n, 1, p_est)

  cstat <- pROC::roc(y, eta, quiet = TRUE, ci = FALSE)
  c_est <- as.vector(cstat$auc)

  out   <- list("p_actual" = round(mean(y), 3),
                  "c_actual" = round(c_est, 3),
                  "mu" = round(mu, 3),
                  "sigma" = round(sigma, 3))
  return(out)
}
