#' Wrapper of the \code{mvtrnd} function written in C++
#' 
#' @import truncnorm
mvtrnd_wrap <- function(a, b, df, NN, veccObj, N, verbose = 0) {
  n <- length(a)
  tmp <- cbind(a, b)
  a <- tmp[, 1]
  b <- tmp[, 2]
  # find tilting parameter beta -----------------------------------
  trunc_expect <- truncnorm::etruncnorm(a, b)
  x0 <- c(trunc_expect, rep(0, n))
  x0[2 * n] <- sqrt(df)
  x0[n] <- log(x0[2 * n])
  solv_idea_5_sp <- nleqslv::nleqslv(
    x = x0, fn = gradpsiT_idea5, veccCondMeanVarObj = veccObj,
    a = a, b = b, nu = df, global = "pwldog", method = "Broyden"
  )
  soln <- solv_idea_5_sp$x
  exitflag <- solv_idea_5_sp$termcd
  if(!(exitflag %in% c(1,2)) || !all.equal(solv_idea_5_sp$fvec, rep(0, length(x0)))){
    warning('Did not find a solution to the nonlinear system in `pmvt`!')
  }
  soln[n] <- exp(soln[n])
  x <- soln[1:n]
  beta <- soln[(n + 1):(2 * n)]
  eta <- beta[n]
  # psi_star --------------------------------------
  psi_star <- psiT_wrapper(x,
    beta = beta, df,
    veccCondMeanVarObj = veccObj,
    a = a, b = b, NN = NN
  )
  # sample until N samples are collected ------------------
  X <- matrix(0, nrow = n, ncol = N)
  R <- rep(0, N)
  accept <- 0L
  iter <- 0L
  n_sim <- floor(N / 2) * 2
  ntotsim <- 0L
  while (accept < N) { # while # of accepted is less than N
    call <- mvtrnd(
      a, b, df, NN, veccObj$cond_mean_coeff,
      sqrt(veccObj$cond_var), beta, n_sim
    )
    ntotsim <- ntotsim + n_sim
    idx <- rexp(n_sim) > (psi_star - call$logpr) # acceptance tests
    m <- sum(idx)
    if (m > N - accept) {
      m <- N - accept
      idx <- which(idx)[1:m]
    }
    if (m > 0) {
      X[, (accept + 1):(accept + m)] <-
        t(call$X_trans[idx, , drop = FALSE]) # accumulate accepted
      R[(accept + 1):(accept + m)] <- call$r[idx]
    }
    accept <- accept + m # keep track of # of accepted
    iter <- iter + 1L # keep track of while loop iterations
    if (verbose) {
      cat("Iteration", iter, "n_sim", n_sim, "n_accept", m, "\n")
    }
    n_sim <- min(c(1e6, N, ceiling(n_sim / m * (N - accept))))
    n_sim <- ceiling(n_sim / 2) * 2
    if ((ntotsim > 1e4) && (accept / ntotsim < 1e-3)) { # if iterations are getting large, give warning
      warning("Acceptance probability smaller than 0.001")
    } else if (iter > 1e5) { # if iterations too large, seek approximation only
      if (accept == 0) {
        stop("Could not sample from truncated Normal - check input")
      } else if (accept > 1) {
        X <- X[, 1:accept]
        warning("Sample of size smaller than N returned.")
      }
    }
  }
  return(t(t(X) / R * sqrt(df)))
}


# TEST -------------------------------------------------------
# library(GpGp)
# library(VeccTMVN)
# library(truncnorm)
# library(TruncatedNormal)
# ## example MVN probabilities --------------------------------
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1 * n2
# locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# covparms <- c(2, 0.3, 0.01)
# cov_name <- "matern15_isotropic"
# cov_mat <- get(cov_name)(covparms, locs)
# a <- rep(-Inf, n)
# b <- rep(-2, n)
# mu <- rep(0, n)
# N <- 1e3
# ## Sample with TruncatedNormal -----------------------
# samp_TN <- TruncatedNormal::mvrandn(a, b, cov_mat, n = N, mu = mu)
# ## Vecc approx objs --------------------------------
# m <- 30 # num of nearest neighbors
# ord <- VeccTMVN::Vecc_reorder(
#   a, b, m, locs, cov_name,
#   covparms
# )$order
# a_vecc_ord <- a[ord]
# b_vecc_ord <- b[ord]
# locs_vecc_ord <- locs[ord, , drop = FALSE]
# NN <- GpGp::find_ordered_nn(locs_vecc_ord, m)
# vecc_obj <- vecc_cond_mean_var_sp(NN,
#                                   locs = locs_vecc_ord, covName = cov_name,
#                                   covParms = covparms
# )
# ## Sample with mvnrnd_wrap -----------------------
# samp_Vecc_ord <- mvnrnd_wrap(
#   a_vecc_ord, b_vecc_ord, mu,
#   NN = NN, veccObj = vecc_obj, N = N, verbose = 1
# )
# ord_rev <- integer(n)
# ord_rev[ord] <- 1 : n
# samp_Vecc <- samp_Vecc_ord[ord_rev, , drop = FALSE]
# ## Visual comparison of two TMVN samples -------------------------------
# par(mfrow = c(1, 2))
# hist(samp_Vecc)
# hist(samp_TN)
