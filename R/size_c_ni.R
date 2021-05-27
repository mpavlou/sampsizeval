
size_c_ni <- function(mu, sigma, se_c) {

  varc <- se_c^2


  # For the numerical integration method a distribution for the linear predictor
  # needs to be assumed. In this calculation we assume marginal normality for
  # the linear predictor

  # f1 and f2 are functions that contain the integrad in equation (5) and (6)
  # as given by Gail and Pfeiffer (2005)


  f1 <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 /
                               (2 * s**2)) * (1  + exp(-x)) ^ (-1)

  f0 <- function(x, mu, s) 1 / (2 * s**2 * pi)**0.5 * exp(- (x - mu)**2 /
                               (2 * s**2)) * (1 - (1 + exp(-x)) ^ (-1))


  # We now get the cumulative distribution of eta_0 and eta1
  # A vector of values to cover the range of possible values of eta_0
  # and eta_0

  # Pre-processing to narrow down the integration range

  x <- seq(- 14, 14, 0.1)

  p_eta0 <- NULL ; p_eta1 <- NULL

  #Numerical integration to get P(eta_0<x) and P(eta_1<x)
  #Equations for Gail and Pfeiffer
  denom1 <- stats::integrate(f1, - Inf, Inf, mu = mu, s = sigma,
                             subdivisions = 1000L)$value
  denom0 <- stats::integrate(f0, - Inf, Inf, mu = mu, s = sigma,
                             subdivisions = 1000L)$value

  p_eta1 <- sapply(x,function(u) stats::integrate(f1, lower = -Inf, mu = mu,
                        s = sigma,upper = u)$value)/denom1
  p_eta0 <- sapply(x,function(u) stats::integrate(f0, lower = -Inf, mu = mu,
                        s = sigma,upper = u)$value)/denom0

  d0 <- data.frame(cbind(x, p_eta0))
  d1 <- data.frame(cbind(x, p_eta1))

  d01 <- subset(d0, p_eta0 < (1 - 0.00001) & p_eta0 > 0.00001)
  d11 <- subset(d1, p_eta1 < (1 - 0.00001) & p_eta1 > 0.00001)

  a <- d01$x[1]
  b <- d11$x[nrow(d11)]

  ## Actual numerical integration starts here

  step <- 0.00005

  x <- seq(a, b, step)

  p_eta0 <- NULL
  p_eta1 <- NULL

  #Numerical integration to get P(eta_0<x) and P(eta_1<x)

  #Equations from Gail and Pfeiffer (2005)


  p_eta1 <- sapply(x,function(u) stats::integrate(f1, lower = -Inf, mu = mu,
                                           s = sigma,upper = u)$value)/denom1
  p_eta0 <- sapply(x,function(u) stats::integrate(f0, lower = -Inf, mu = mu,
                                           s = sigma,upper = u)$value)/denom0

  p_eta0[1] <- 0
  p_eta1[1] <- 0
  if (p_eta0[length(x)] != 1) p_eta0[length(x)] <- 1
  if (p_eta1[length(x)] != 1) p_eta1[length(x)] <- 1

  x<-round(x,5)
  d0<-data.frame(x,p_eta0)
  d1<-data.frame(x,p_eta1)


  # Inverse CDF to sample from the distribution of eta_0 and eta_1

  nsamps <- 2000000

  u        <- sort(stats::runif(nsamps))
  eta0     <- pracma::interp1(x = sort(p_eta0), y = x, xi = u)
  eta0     <- sort(plyr::round_any(eta0,step))
  eta0     <- round(eta0,5)
  eta0s<-data.frame(eta0)


  u        <- sort(stats::runif(nsamps))
  eta1     <- pracma::interp1(x = sort(p_eta1), y = x, xi = u)
  eta1     <- sort(plyr::round_any(eta1,step))
  eta1     <- round(eta1,5)
  eta1s    <- data.frame(eta1)

  #Note: Distribution of eta_0 and eta_1 is approximately Normal, as expected


  #system.quant(prob <- sapply(eta1,function(u) integrate(f0, lower = -Inf, mu = mu,
  #                                        s = sigma,upper = u)$value)/denom0)

  # Compute EK=E(P(eta_0<eta_1))

  merged <- dplyr::left_join(eta1s,d0,by=c("eta1"="x"))
  prob   <- merged$p_eta0
  e_k2  <-  mean(prob^2)

  #prob <- mean(ifelse(sample(eta0,nsamps)<sample(eta1,nsamps),1,0))

  #Calculate  the true C and p for these values of mu and sigma
  c_ni      <- mean(prob)
  p_ni      <- denom1
  c_true    <- c_ni
  p_true    <- p_ni

  # Compute EG=E(P(eta_1<eta_0))
  #system.quant(prob <- sapply(eta0,function(u) integrate(f1, lower = -Inf, mu = mu,
                                                      #  s = sigma,upper = u)$value)/denom1)

  # without additional integration

  merged <- dplyr::left_join(eta0s,d1,by=c("eta0"="x"))
  prob   <- merged$p_eta1
  e_g2   <- mean((1 - prob)^2)

  #Enter in the final formula
  n         <- (1 / varc) * ((1 - p_true) * e_k2 + p_true * e_g2 - c_true^2) /
               (p_true * (1 - p_true))
  n         <- ceiling(n)
  #events    <- ceiling(n * p)

  return(n)
}

#system.time(a<-sampsizeval(p=0.1, c=0.75, se_c=0.025, se_cs =0.1, se_cl = 0.1,c_ni=TRUE))
#a
