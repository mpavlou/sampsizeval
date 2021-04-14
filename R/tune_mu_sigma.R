#' Find mu and sigma that correspond to antcipated C and p (when C>=0.8)
#'
#' The purpose of this function is to check whether the values of mu
#' and sigma as given by the formulae in the paper
#' give actual value of C sufficiently close to the anticipated
#' values. If not, mu and sigma are fine-tuned to provide actual values of c and
#' p sufficiently close to the anticipated ones. This is relevant for the
#' version of the formulae that require the use of numerical integration, when
#' C>=0.8.
#'
#' @param p The anticipated outcome prevalence
#' @param c The anticipated C-statistic
#'
#' @return  The actual prevalence and c-statistic, for the updated values of
#' mu and sigma based on the tuning factor fc.

tune_mu_sigma <- function(p, c) {

  actual <- actual_values(p = p , c = c, fc = 1)
  fc     <- 1.0

  while (abs(actual$c_actual - c) > 0.002) {
    fc     <- fc + 0.005
    actual <- actual_values(p = p , c = c, fc = fc)

  }

  out <- list("p_actual" = actual$p_actual,
              "c_actual" = actual$c_actual,
              "mu" = actual$mu,
              "sigma" = actual$sigma,
              "fc" = fc)
  return(out)
}
