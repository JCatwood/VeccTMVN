#' Gradient of the psiT function computed using sparse A. For this function,
#' \code{xAndBeta} should be a vector of length 2 * n, representing x_n and
#' beta_n
#' 
#' @noRd
#' @import TruncatedNormal
gradpsiT_idea5 <- function(xAndBeta, veccCondMeanVarObj, a, b, nu) {
  n <- length(a)
  x <- beta <- co <- rep(0, n)
  x[-n] <- xAndBeta[1:(n - 1)]
  r <- exp(xAndBeta[n])
  beta[-n] <- xAndBeta[(n + 1):(2 * n - 1)]
  eta <- xAndBeta[2 * n]
  a <- a / sqrt(nu)
  b <- b / sqrt(nu)

  D <- sqrt(veccCondMeanVarObj$cond_var)
  mu_c <- as.vector(veccCondMeanVarObj$A %*% (x))
  a_tilde_shift <- (a * r - mu_c) / D - beta
  b_tilde_shift <- (b * r - mu_c) / D - beta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  Psi <- pl - pu
  # compute grad ------------------------------------------------
  dpsi_dx <- as.vector(Matrix::t(veccCondMeanVarObj$A) %*%
    (beta / D + Psi / D)) - beta / D
  dpsi_dbeta <- beta - (x - mu_c) / D + Psi
  a[is.infinite(a)] <- 0
  b[is.infinite(b)] <- 0
  dpsi_dr <- (nu - 1) / r - eta + sum(b * pu - a * pl)
  dpsi_de <- eta - r + exp(stats::dnorm(eta, log = TRUE) - 
                             TruncatedNormal::lnNpr(-eta, Inf))

  c(dpsi_dx[-n], dpsi_dr, dpsi_dbeta[-n], dpsi_de)
}

#' psiT function that wraps the C function \code{psi}
#' 
#' @noRd
psiT_wrapper <- function(x, beta, nu, veccCondMeanVarObj, NN, a, b) {
  n <- length(a)
  r <- x[n] # note that this is after exponential
  eta <- beta[n]
  x[n] <- 0
  beta[n] <- 0
  a <- a / sqrt(nu) * r
  b <- b / sqrt(nu) * r

  psi(
    a, b, NN, rep(0, length(n)), veccCondMeanVarObj$cond_mean_coeff,
    sqrt(veccCondMeanVarObj$cond_var), beta, x
  ) + log(2 * pi) / 2 - lgamma(nu / 2) - (0.5 * nu - 1) * log(2) +
    0.5 * sum(eta * eta) - r * eta + (nu - 1) * log(r) +
    TruncatedNormal::lnNpr(-eta, Inf)
}
