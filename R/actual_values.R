
actual_values <- function(p, c, fc = 1) {
  nevents <- 400000
  n       <- nevents / p

  quickcstat<-function(y,pred,seed=1){
    set.seed(seed)
    casepred=pred[y==1]
    conpred=pred[y==0]
    conpred=conpred[sample(length(conpred),length(casepred),replace=FALSE)]
    auc.true=sum(casepred>conpred)/length(casepred)
    return(auc.true)
  }

  sigma_c <- sqrt(2) * stats::qnorm(c) * fc
  mu      <- 0.5 * (2 * p - 1) * (sigma_c^2) + log(p / (1 - p))
  sigma   <- sqrt((sigma_c^2) * (1 + p * (1 - p) * (sigma_c^2)))

  eta   <- stats::rnorm(n, mu, sigma)
  p_est <- (1 + exp(- eta)) ^ (-1)
  y     <- stats::rbinom(n, 1, p_est)

  #cstat <- pROC::roc(y, eta, quiet = TRUE, ci = FALSE)
  #c_est <- as.vector(cstat$auc)

  c_est <- quickcstat(y,eta,seed=1)

  out   <- list("p_actual" = round(mean(y), 3),
                  "c_actual" = round(c_est, 3),
                  "mu" = round(mu, 3),
                  "sigma" = round(sigma, 3))
  return(out)
}
