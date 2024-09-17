#' Applying the multi-level Monte Carlo (MLMC) technique to the pmvn function
#' The function uses \code{NLevel1 = 1} for \code{m = m2} and the same
#' exponential tilting parameter as \code{m = m1} to compute one MC estimate.
#' This MC estimate is used to correct the bias from the Vecchia approximation
#'
#' @import truncnorm
#' @importFrom stats runif
#' @importFrom TruncatedNormal cholperm
#' @importFrom utils getFromNamespace
#' 
#' @param lower lower bound vector for TMVN
#' @param upper upper bound vector for TMVN
#' @param mean MVN mean
#' @param locs location (feature) matrix n X d
#' @param covName covariance function name from the `GpGp` package
#' @param covParms parameters for `covName`
#' @param m1 the smaller Vecchia conditioning set size for Level 1 MC
#' @param m2 the bigger Vecchia conditioning set size for Level 2 MC
#' @param sigma dense covariance matrix, not needed when `locs` is not null
#' @param reorder whether to reorder integration variables. `0` for no,
#' `1` for FIC-based univariate ordering, `2` for Vecchia-based univariate
#' ordering, and `3` for the reordering implemented in TruncatedNormal, 
#' which appeared faster than `2`
#' @param NLevel1 first level Monte Carlo sample size
#' @param NLevel2 second level Monte Carlo sample size
#' @param verbose verbose or not
#' @param retlog TRUE or FALSE for whether to return loglk or not
#' @param ... could be
#' m_ord for conditioning set size for reordering
#' @return estimated MVN probability and estimation error
#'
#' @export
pmvn_MLMC <- function(lower, upper, mean, locs = NULL, covName = "matern15_isotropic",
                      covParms = c(1.0, 0.1, 0.0), m1 = 30, m2 = 100, sigma = NULL, reorder = 0,
                      NLevel1 = 12, NLevel2 = 1e4, verbose = FALSE, retlog = FALSE, ...) {
  # standardize the input MVN prob -----------------------------
  lower <- lower - mean
  upper <- upper - mean
  if (is.null(sigma)) {
    n <- nrow(locs)
    use_sigma <- FALSE
    margin_sd <- sqrt(covParms[1])
    upper <- upper / margin_sd
    lower <- lower / margin_sd
    covParms[1] <- 1
  } else {
    n <- nrow(sigma)
    use_sigma <- TRUE
    margin_sd <- sqrt(diag(sigma))
    upper <- upper / margin_sd
    lower <- lower / margin_sd
    sigma <- t(t(sigma / margin_sd) / margin_sd)
  }
  if (any(upper < lower)) {
    stop("Invalid MVN probability. Truncated marginal
         probabilities have negative value(s)\n")
  }
  if (m1 > m2) {
    stop("m1 should be greater than m2 when using
         Multi-level Monte Carlo(MLMC)\n")
  }
  lower_upper <- matrix(0, n, 2)
  lower_upper[, 1] <- lower
  lower_upper[, 2] <- upper
  lower <- lower_upper[, 1]
  upper <- lower_upper[, 2]
  # reorder --------------------------------
  if (is.null(list(...)[["m_ord"]])) {
    m_ord <- m1
  } else {
    m_ord <- list(...)[["m_ord"]]
  }
  if (reorder == 1) {
    if (use_sigma) {
      ord <- FIC_reorder_univar(lower, upper, m_ord, covMat = sigma)
    } else {
      ord <- FIC_reorder_univar(
        lower, upper, m_ord, locs, covName,
        covParms
      )
    }
    lower <- lower[ord]
    upper <- upper[ord]
    if (use_sigma) {
      sigma <- sigma[ord, ord, drop = FALSE]
    } else {
      locs <- locs[ord, , drop = FALSE]
    }
  } else if (reorder == 2) {
    if (use_sigma) {
      ord <- Vecc_reorder(lower, upper, m_ord, covMat = sigma)$order
    } else {
      ord <- Vecc_reorder(
        lower, upper, m_ord, locs, covName, covParms
      )$order
    }
    lower <- lower[ord]
    upper <- upper[ord]
    if (use_sigma) {
      sigma <- sigma[ord, ord, drop = FALSE]
    } else {
      locs <- locs[ord, , drop = FALSE]
    }
  } else if (reorder == 3) {
    if (!use_sigma) {
      cov_func_GpGp <- utils::getFromNamespace(covName, "GpGp")
      sigma <- cov_func_GpGp(covParms, locs)
    }
    ord <- TruncatedNormal::cholperm(sigma, lower, upper)$perm
    lower <- lower[ord]
    upper <- upper[ord]
    if (use_sigma) {
      sigma <- sigma[ord, ord, drop = FALSE]
    } else {
      locs <- locs[ord, , drop = FALSE]
    }
  }
  # find nearest neighbors for Vecchia --------------------------------
  if (use_sigma) {
    NN_m2 <- find_nn_corr(sigma, m2)
    NN_m1 <- NN_m2[, 1:(m1 + 1)]
  } else {
    NN_m2 <- GpGp::find_ordered_nn(locs, m2)
    NN_m1 <- NN_m2[, 1:(m1 + 1)]
  }
  # find Vecchia approx object -----------------------------------
  if (use_sigma) {
    vecc_obj_m1 <- vecc_cond_mean_var_sp(NN_m1, covMat = sigma)
    vecc_obj_m2 <- vecc_cond_mean_var_sp(NN_m2, covMat = sigma)
  } else {
    vecc_obj_m1 <- vecc_cond_mean_var_sp(NN_m1,
      locs = locs, covName = covName,
      covParms = covParms
    )
    vecc_obj_m2 <- vecc_cond_mean_var_sp(NN_m2,
      locs = locs, covName = covName,
      covParms = covParms
    )
  }
  # find tilting parameter beta -----------------------------------
  trunc_expect <- truncnorm::etruncnorm(lower, upper)
  x0 <- c(trunc_expect, rep(0, n))
  solv_idea_5_sp <- stats::optim(
    x0,
    fn = function(x, ...) {
      ret <- grad_jacprod_jacsolv_idea5(x, ...,
        retJac = FALSE,
        retProd = FALSE, retSolv = FALSE
      )
      0.5 * sum((ret$grad)^2)
    },
    gr = function(x, ...) {
      ret <- grad_jacprod_jacsolv_idea5(x, ...,
        retJac = FALSE,
        retProd = TRUE, retSolv = FALSE
      )
      ret$jac_grad
    },
    method = "L-BFGS-B",
    veccCondMeanVarObj = vecc_obj_m1,
    a = lower, b = upper, verbose = verbose,
    lower = c(lower, rep(-Inf, n)), upper = c(upper, rep(Inf, n)),
    control = list(maxit = 500)
  )
  if (verbose) {
    cat(
      "Gradient norm at the optimal beta is", sqrt(2 * solv_idea_5_sp$value),
      "\n"
    )
  }
  if (any(solv_idea_5_sp$par[1:n] < lower) ||
    any(solv_idea_5_sp$par[1:n] > upper)) {
    warning("Optimal x is outside the integration region during minmax tilting\n")
  }
  beta <- solv_idea_5_sp$par[(n + 1):(2 * n)]
  # compute MVN probs and est error ---------------------------------
  seed <- round(runif(1) * 1e6)
  set.seed(seed)
  exp_psi_m1 <- sample_psi_idea5_cpp(vecc_obj_m1, lower, upper,
    beta = beta, N_level1 = NLevel1,
    N_level2 = NLevel2
  )
  set.seed(seed)
  exp_psi_m2 <- sample_psi_idea5_cpp(vecc_obj_m2, lower, upper,
    beta = beta, N_level1 = 1,
    N_level2 = NLevel2
  )
  if (retlog) {
    exponent <- min(exp_psi_m1[[2]])
    # correction using exp_psi_m2
    log_est_prob <- exponent +
      log(mean(exp_psi_m1[[1]] * exp(exp_psi_m1[[2]] - exponent)) +
        exp_psi_m2[[1]][1] * exp(exp_psi_m2[[2]][1] - exponent) -
        exp_psi_m1[[1]][1] * exp(exp_psi_m1[[2]][1] - exponent))
    return(log_est_prob)
  } else {
    exp_psi <- exp_psi_m1[[1]] * exp(exp_psi_m1[[2]])
    exp_psi_m2 <- exp_psi_m2[[1]] * exp(exp_psi_m2[[2]])
    # correction using exp_psi_m2
    est_prob <- mean(exp_psi) + (exp_psi_m2[1] - exp_psi[1])
    est_prob_err <- stats::sd(exp_psi) / sqrt(NLevel1)
    attr(est_prob, "error") <- est_prob_err
    return(est_prob)
  }
}


# TEST -------------------------------------------------------

# library(VeccTMVN)
# library(TruncatedNormal)
# library(GpGp)
#
# ## example MVN probabilities --------------------------------
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1 * n2
# locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# covparms <- c(2, 0.3, 0)
# cov_mat <- matern15_isotropic(covparms, locs)
# a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2 - 4)
# b_list <- list(rep(-2, n), rep(1, n), -runif(n) * 2)
#
# ## Compute MVN probs --------------------------------
# N_level1 <- 12 # Level 1 MC size
# N_level2 <- 1e4 # Level 2 MC size
# m <- 30 # num of nearest neighbors
# for (i in 1:length(a_list)) {
#   a <- a_list[[i]]
#   b <- b_list[[i]]
#   ### Compute MVN prob with idea V -----------------------
#   est_Vecc <- VeccTMVN::pmvn(a, b, 0, locs, covName = "matern15_isotropic",
#                              covParms = covparms, m = m, verbose = FALSE)
#   est_Vecc_MLMC <- VeccTMVN::pmvn_MLMC(a, b, 0, locs, covName = "matern15_isotropic",
#                                        covParms = covparms, m1 = m, m2 = 2*m, verbose = FALSE)
#   est_TN <- TruncatedNormal::pmvnorm(rep(0, n), cov_mat, lb = a, ub = b)
#   cat(
#     "est_Vecc", est_Vecc, "err_Vecc", attributes(est_Vecc)$error, "\n"
#   )
#   cat(
#     "est_Vecc_MLMC", est_Vecc_MLMC, "err_Vecc", attributes(est_Vecc_MLMC)$error,
#     "\n"
#   )
#   cat(
#     "est_TN", est_TN, "err_TN", attributes(est_TN)$relerr * est_TN, "\n"
#   )
# }
