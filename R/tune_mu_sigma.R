

tune_mu_sigma <- function(p, c) {

  actual <- actual_values(p = p, c = c, fc = 1)
  fc     <- 1.0

  while (abs(actual$c_actual - c) > 0.002) {

    if (c>0.92) fc <- fc+0.01 else fc  <- fc + 0.005
    actual <- actual_values(p = p, c = c, fc = fc)

  }

  out <- list("p_actual" = actual$p_actual,
              "c_actual" = actual$c_actual,
              "mu" = actual$mu,
              "sigma" = actual$sigma,
              "fc" = fc)
  return(out)
}

