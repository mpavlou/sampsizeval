
size_c_app <- function(p, c, se_c) {

  varc <- se_c^2
  n            <- (c - 2 * sn::T.Owen(-stats::qnorm(c), 1 / sqrt(3)) - c^2) /
                  ((p - p^2) * varc)
  n            <- ceiling(n)
  events       <- ceiling(n * p)
  return(n)
}
