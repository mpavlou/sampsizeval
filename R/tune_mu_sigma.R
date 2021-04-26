

tune_mu_sigma <- function(p, c) {

  actual <- actual_values(p = p, c = c, fc = 1)
  fc     <- 1.0

  while (abs(actual$c_actual - c) > 0.0025) {

    fc <- fc+0.01
    actual <- actual_values(p = p, c = c, fc = fc)

  }

  out <- list("p_actual" = actual$p_actual,
              "c_actual" = actual$c_actual,
              "mu" = actual$mu,
              "sigma" = actual$sigma,
              "fc" = fc)
  return(out)
}

#system.time(a<-tune_mu_sigma(0.25,0.81))
#a

