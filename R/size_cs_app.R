
size_cs_app  <- function(p, c, se_cs) {
varcs    <- se_cs^2
a        <- 2 * p * (1 - p) * stats::qnorm(c)^2
n        <- 1 / varcs * (1 / a + 2)
n        <- ceiling(n)
events   <- ceiling(n * p)
return(n)
}
