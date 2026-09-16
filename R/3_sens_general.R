# --------------------------- Helper Functions --------------------------------#

#' Project a point onto the ellipsoid
#' E = {beta: (beta_med - beta)' Sigma (beta_med - beta) <= r^2}
#'
#' @param z A \eqn{d_X \times 1} vector to project.
#' @param ellipsoid A list computed from `precompute_ellipsoid()`.
#'
#' @return A \eqn{d_X \times 1} vector.
#' @noRd
project_ellipsoid <- function(z, ellipsoid) {
  beta_med <- ellipsoid$beta_med
  radius <- ellipsoid$radius
  Q <- ellipsoid$Q # matrix of eigenvectors (dX x dX matrix)
  lam <- ellipsoid$lam # eigenvalues (dX x 1 vector)
  lam_sqrt <- ellipsoid$lam_sqrt # square root of eigenvalues (dX x 1 vector)

  if (radius <= 0) {
    return(as.numeric(beta_med))
  }

  delta <- beta_med - z # dX x 1 vector
  # w = lam^{1/2} Q' delta  (in eigenbasis) (element-wise)
  w <- lam_sqrt * drop(crossprod(Q, delta)) # dX x 1 vector

  # Check if z is already feasible (unconstrained gamma = w)
  if (sum(w^2) <= radius^2) {
    return(z)
  }

  # Bisection to find lambda > 0 such that sum((w_i / (1 + lambda * lam_i))^2) = radius^2
  # The LHS is strictly decreasing from sum(w^2) > radius^2 to 0 as lambda -> Inf.
  phi <- function(lambda) {
    sum((w / (1 + lambda * lam))^2) - radius^2
  }

  # Upper bound: need 1 + lambda * lam_lb large enough
  lam_l <- 0
  lam_u <- 1
  while (phi(lam_u) > 0) lam_u <- lam_u * 2

  # Bisection
  for (i in seq_len(100)) {
    lam_mid <- (lam_l + lam_u) / 2
    if (phi(lam_mid) > 0) lam_l <- lam_mid else lam_u <- lam_mid
    if (lam_u - lam_l < 1e-10 * max(1, lam_u)) break
  }
  lambda_star <- (lam_l + lam_u) / 2

  # Recover gamma in eigenbasis, then beta
  gamma_eig <- w / (1 + lambda_star * lam) # dX x 1 vector in eigenbasis
  # beta = beta_med - Sigma^{-1/2} gamma = beta_med - Q (gamma_eig / sqrt(lam))
  beta_proj <- drop(beta_med - Q %*% (gamma_eig / lam_sqrt))
  return(beta_proj)
}

#' Sensitivity of general functions of OLS coefficients
#'
#' `proj_grad_descent()` calculates the extremal values for general functions of
#' OLS coefficients \eqn{f_\text{g}(\beta)} by projected gradient descent.
#'
#' @inheritParams sens_general
#'
#' @return A list with `beta` (solution), `f_val` (objective), `history` (data.frame).
#' @noRd
proj_grad_descent <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid = NULL,
                              f, grad_f, beta_init, step_size = 0.01,
                              max_iter = 5000, tol = 1e-8, verbose = 0) {
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }
  beta <- project_ellipsoid(beta_init, ellipsoid) # ensure feasibility
  fval <- f(beta)
  history <- data.frame(iter = 0:max_iter, f_val = fval, step = NA, change = NA)
  history$f_val[1] <- fval

  for (k in seq_len(max_iter)) {
    g <- grad_f(beta)

    #--- Step size selection ---#
    if (identical(step_size, "backtrack")) {
      eta <- 1.0
      armijo <- 0.5
      shrink <- 0.5
      f0 <- fval
      repeat {
        candidate <- project_ellipsoid(beta - eta * g, ellipsoid)
        if (f(candidate) <= f0 - armijo * sum(g * (beta - candidate))) break
        eta <- eta * shrink
        if (eta < 1e-15) break
      }
    } else {
      eta <- step_size
    }

    beta_update <- project_ellipsoid(beta - eta * g, ellipsoid)
    fval_update <- f(beta_update)
    change <- sqrt(sum((beta_update - beta)^2))

    history$f_val[k + 1] <- fval_update
    history$step[k + 1] <- eta
    history$change[k + 1] <- change

    beta <- beta_update
    fval <- fval_update

    if (verbose > 0 && k %% verbose == 0) {
      cat(sprintf(
        "iter %5d | f = %.8f | change = %.2e | step = %.2e\n",
        k, fval, change, eta
      ))
    }
    if (change < tol) {
      if (verbose > 0) cat("Converged at iteration", k, "\n")
      break
    }
  }

  keep <- seq_len(k + 1)
  history <- history[keep, ]

  list(beta = beta, f_val = fval, history = history)
}

#' Check the validity of the input of general functions
#'
#' @inheritParams sens_general
#'
#' @noRd
check_input_general <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid,
                                f, grad_f, beta_init, step_size, max_iter, tol) {
  check_input_ellipsoid(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid
  )
  dX <- ncol(as.matrix(X))

  # f and grad_f
  if (!is.function(f) || !is.function(grad_f)) {
    stop("`f` and `grad_f` must both be functions.")
  }

  # beta_init
  if (!is.numeric(beta_init) || length(beta_init) != dX || anyNA(beta_init) || any(!is.finite(beta_init))) {
    stop("`beta_init` must be a finite numeric vector of length ncol(X).")
  }

  # step_size
  if (!identical(step_size, "backtrack")) {
    if (!is.numeric(step_size) || length(step_size) != 1 || !is.finite(step_size) || step_size <= 0) {
      stop('`step_size` must be a positive finite scalar or "backtrack".')
    }
  }

  # max_iter
  if (!is.numeric(max_iter) || length(max_iter) != 1 || !is.finite(max_iter) || max_iter < 1) {
    stop("`max_iter` must be a positive integer.")
  }

  # tol
  if (!is.numeric(tol) || length(tol) != 1 || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite scalar.")
  }
}

# ----------------------------- Main Functions --------------------------------#

#' Sensitivity of general functions of OLS coefficients
#'
#' `sens_general()` calculates the extremal values for general functions of
#' OLS coefficients \eqn{f_\text{g}(\beta)} by projected gradient descent.
#'
#' @param Y A \eqn{N \times 1} vector.
#' @param X A \eqn{N \times d_X} dataframe/matrix/vector.
#' @param W1 A \eqn{N \times d1} dataframe/matrix/vector or `NULL`.
#' @param bar_rho A scalar in \eqn{[0, 1)}.
#' @param bar_R2 A scalar in \eqn{[0, 1)}.
#' @param ellipsoid A list computed from `precompute_ellipsoid()` (optional).
#' @param f Objective function \eqn{f(\beta)} -> scalar.
#' @param grad_f Gradient function \eqn{\nabla f(\beta)} -> \eqn{d_X \times 1} vector.
#' @param beta_init Initial feasible point (\eqn{d_X \times 1} vector).
#' @param step_size Fixed step size (scalar), or "backtrack" for backtracking line search.
#' @param max_iter Maximum number of iterations.
#' @param tol Convergence tolerance on \eqn{\|\beta^{t+1} - \beta^t\|}.
#' @param verbose Print progress every `verbose` iterations (0 = silent).
#'
#' @return A list of 6 objects where:
#' * `beta_lb` stores the coefficient vector beta attaining the lower bound,
#' * `beta_ub` stores the coefficient vector beta attaining the upper bound,
#' * `obj_lb` stores the lower bound of \eqn{f_{g}(\beta)},
#' * `obj_ub` stores the upper bound of \eqn{f_{g}(\beta)},
#' * `history_lb` stores the history of optimizing lower bound of \eqn{f_{g}(\beta)},
#' * `history_ub` stores the history of optimizing upper bound of \eqn{f_{g}(\beta)}.
#' @export
sens_general <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid = NULL,
                         f, grad_f, beta_init, step_size = 0.01,
                         max_iter = 5000, tol = 1e-8, verbose = 0) {
  check_input_general(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    f = f, grad_f = grad_f, beta_init = beta_init, step_size = step_size,
    max_iter = max_iter, tol = tol
  )
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }

  result_lb <- proj_grad_descent(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    f = f, grad_f = grad_f, beta_init = beta_init, step_size = step_size,
    max_iter = max_iter, tol = tol, verbose = verbose
  )
  beta_lb <- result_lb$beta
  obj_lb <- result_lb$f_val
  history_lb <- result_lb$history

  f_neg <- function(beta) -f(beta)
  grad_f_neg <- function(beta) -grad_f(beta)
  result_ub <- proj_grad_descent(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    f = f_neg, grad_f = grad_f_neg, beta_init = beta_init, step_size = step_size,
    max_iter = max_iter, tol = tol, verbose = verbose
  )

  beta_ub <- result_ub$beta
  obj_ub <- -result_ub$f_val
  history_ub <- result_ub$history
  history_ub$f_val <- -history_ub$f_val

  list(
    beta_lb = beta_lb,
    beta_ub = beta_ub,
    obj_lb = obj_lb,
    obj_ub = obj_ub,
    history_lb = history_lb,
    history_ub = history_ub
  )
}
