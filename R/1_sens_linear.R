# --------------------------- Helper Functions --------------------------------#

#' Residualize the variable by partialing out controls
#'
#' Computes the OLS residual plus intercept from regressing each
#' column of `vars` on an intercept and `controls`.
#'
#' @param vars Variables of interest \eqn{N \times J} dataframe/matrix/vector).
#' @param controls Control variables (\eqn{N \times K} dataframe/matrix/vector or `NULL`).
#'
#' @return A residualized \eqn{N \times J} matrix.
#' @noRd
residualize <- function(vars, controls) {
  # Deal with input
  vars <- as.matrix(vars)
  if (is.null(controls)) {
    vars_resid <- vars
  } else {
    controls <- as.matrix(controls)
    vars_resid <- vars - controls %*% solve(stats::var(controls), stats::cov(controls, vars))
  }
  return(as.matrix(vars_resid))
}

#' Calculate A to the power of s
#'
#' @param A A symmetric and positive definite square matrix/dataframe.
#' @param s A numeric value.
#'
#' @return A matrix of `A` to the power of `s`.
#' @noRd
mat_power <- function(A, s) {
  # Deal with input
  A <- as.matrix(A)
  A <- (A + t(A)) / 2 # enforce symmetry
  ee <- eigen(A, symmetric = TRUE)
  if (min(ee$values) <= 0) {
    stop("matrix `A` is not positive definite.")
  }
  out <- ee$vectors %*% diag(ee$values^s, nrow = length(ee$values)) %*% t(ee$vectors)
  out <- (out + t(out)) / 2 # enforce symmetry
  return(out)
}

#' Calculate the largest canonical correlation between A and B
#'
#' @param A A \eqn{N \times p} matrix/dataframe.
#' @param B A \eqn{N \times q} matrix/dataframe.
#'
#' @return A nonnegative scalar which indicates the largest canonical correlation
#' between `A` and `B`.
#' @noRd
largest_can_cor <- function(A, B) {
  # Deal with input
  A <- as.matrix(A)
  B <- as.matrix(B)
  cc <- stats::cancor(A, B)
  return(cc$cor[1])
}

#' Precompute the ellipsoid parameters
#'
#' @param Y A \eqn{N \times 1} vector.
#' @param X A \eqn{N \times dX} dataframe/matrix/vector.
#' @param W1 A \eqn{N \times d1} dataframe/matrix/vector or `NULL`.
#' @param bar_rho A scalar in \eqn{[0, 1)}.
#' @param bar_R2 A scalar in \eqn{[0, 1)}.
#'
#' @return A list of 10 objects where:
#' * `beta_med`: OLS coefficients on \eqn{X} in the medium regression (\eqn{dX \times 1}).
#' * `Sigma`: var-cov matrix of \eqn{X^{\perp W_1}}.
#' * `Sigma_inv`: the inverse of `Sigma`.
#' * `Sigma_half`: `Sigma` to the power of 1/2.
#' * `Sigma_nhalf`: `Sigma` to the power of -1/2.
#' * `sigma2_Y_perp_XW1`: variance of \eqn{Y^{\perp X, W1}}.
#' * `radius`: ellipsoid radius scalar.
#' * `Q`: eigenvector matrix of `Sigma` (\eqn{dX \times dX}).
#' * `lam`: eigenvalues of `Sigma` (\eqn{dX \times 1}).
#' * `lam_sqrt`: square roots of eigenvalues.
#' @noRd
precompute_ellipsoid <- function(Y, X, W1, bar_rho, bar_R2) {
  # Deal with input
  Y <- as.numeric(Y)
  X <- as.matrix(X)
  if (!is.null(W1)) {
    W1 <- as.matrix(W1)
  }
  dX <- ncol(X)
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(dX))
  }

  # Residualize X and Y
  X_perp_W1 <- residualize(vars = X, controls = W1)
  Y_perp_W1 <- residualize(vars = Y, controls = W1)
  Y_perp_XW1 <- residualize(vars = Y, controls = if (is.null(W1)) X else cbind(X, W1))

  # Solve for beta_med
  beta_med <- as.numeric(
    solve(stats::cov(X_perp_W1), stats::cov(X_perp_W1, Y_perp_W1))
  )
  names(beta_med) <- colnames(X)

  # Calculate the var-cov matrix of X^{\perp W1}
  Sigma <- stats::cov(X_perp_W1)
  Sigma_inv <- mat_power(Sigma, -1) # Sigma^{-1}
  Sigma_half <- mat_power(Sigma, 1 / 2) # Sigma^{1/2}
  Sigma_nhalf <- mat_power(Sigma, -1 / 2) # Sigma^{-1/2}

  # Calculate the variance of Y^{\perp X, W1}
  sigma2_Y_perp_XW1 <- as.numeric(stats::var(Y_perp_XW1))

  # Calculate the ellipsoid radius
  radius <- bar_rho / sqrt(1 - bar_rho^2) * sqrt(sigma2_Y_perp_XW1 * bar_R2)

  # eigendecomposition: Sigma = Q lam Q'
  ee <- eigen(Sigma, symmetric = TRUE)
  Q <- ee$vectors # matrix of eigenvectors (dX x dX matrix)
  lam <- ee$values # eigenvalues (dX x 1 vector)
  lam_sqrt <- sqrt(lam) # square root of eigenvalues (dX x 1 vector)

  list(
    beta_med = beta_med, Sigma = Sigma, Sigma_inv = Sigma_inv,
    Sigma_half = Sigma_half, Sigma_nhalf = Sigma_nhalf,
    sigma2_Y_perp_XW1 = sigma2_Y_perp_XW1, radius = radius,
    Q = Q, lam = lam, lam_sqrt = lam_sqrt
  )
}

#' Solve Gurobi model
#'
#' @param model model list for Gurobi.
#' @param OutputFlag indicator of controlling solver output.
#' * 0: suppress solver output (default)
#' * 1: enable solver output
#'
#' @return A list storing Gurobi optimization result.
#' @noRd
solve_gurobi_model <- function(model, OutputFlag = 0) {
  lp <- gurobi::gurobi(model, params = list(OutputFlag = OutputFlag))

  # Check return status
  if (lp$status != "OPTIMAL") {
    stop(paste("Gurobi optimization failed with status: ", lp$status))
  }
  return(lp)
}

#' Sensitivity of linear functions of OLS coefficients
#'
#' `solve_linear_gurobi()` solves the optimization problem for linear functions
#' of OLS coefficients \eqn{f_L(\beta) = c^{\intercal}\beta + c_0} by gurobi.
#'
#' @inheritParams sens_linear
#' @param sense `"min"` or `"max"`
#'
#' @return A gurobi object.
#' @noRd
solve_linear_gurobi <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid = NULL, c, c0, sense) {
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }
  beta_med <- ellipsoid$beta_med
  Sigma <- ellipsoid$Sigma
  radius <- ellipsoid$radius

  dX <- length(beta_med)

  # Build Gurobi model
  model <- list()
  model$modelsense <- sense
  model$obj <- c
  model$objcon <- c0
  model$lb <- rep(-Inf, dX)
  model$ub <- rep(Inf, dX)

  # No linear constraints
  model$A <- matrix(0, nrow = 0, ncol = dX)
  model$rhs <- numeric(0)
  model$sense <- character(0)

  # Quadratic constraint
  # beta' Sigma beta - 2 beta_med' Sigma beta <= r^2 - beta_med' Sigma beta_med
  qc <- list()
  qc$Qc <- Sigma # dX x dX matrix
  qc$q <- -2 * as.numeric(t(Sigma) %*% beta_med) # dX x 1 vector
  qc$rhs <- radius^2 - as.numeric(t(beta_med) %*% Sigma %*% beta_med) # scalar
  qc$sense <- "<"
  model$quadcon <- list(qc)

  # Solve gurobi model
  result <- solve_gurobi_model(model)

  return(result)
}

#' Check the validity of the input of ellipsoid
#'
#' @inheritParams sens_linear
#'
#' @noRd
check_input_ellipsoid <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid) {
  # Y
  if (is.null(Y) || !is.numeric(Y) || !is.vector(Y)) {
    stop("`Y` must be a numeric vector for the outcome.")
  }
  if (anyNA(Y)) {
    stop("`Y` must not contain missing values.")
  }
  n <- length(Y)

  # X
  if (!(is.matrix(X) || is.data.frame(X) || !is.vector(X))) {
    stop("`X` must be a matrix, data frame, or vector.")
  }
  X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("`X` must be numeric.")
  }
  if (nrow(X) != n) {
    stop("`X` must have the same number of rows as the length of `Y`.")
  }
  if (anyNA(X)) {
    stop("`X` must not contain missing values.")
  }
  dX <- ncol(X)

  # W1
  if (!is.null(W1)) {
    if (!(is.matrix(W1) || is.data.frame(W1) || !is.vector(W1))) {
      stop("`W1` must be a matrix, data frame, vector, or `NULL`.")
    }
    W1 <- as.matrix(W1)
    if (!is.numeric(W1)) {
      stop("`W1` must be numeric.")
    }
    if (nrow(W1) != n) {
      stop("`W1` must have the same number of rows as the length of `Y`.")
    }
    if (anyNA(W1)) {
      stop("`W1` must not contain missing values.")
    }
  }

  # bar_rho
  if (!is.numeric(bar_rho) || length(bar_rho) != 1L || is.na(bar_rho) || !is.finite(bar_rho) || bar_rho < 0 || bar_rho >= 1) {
    stop("`bar_rho` must be a finite scalar in [0, 1).")
  }

  # bar_R2
  if (!is.numeric(bar_R2) || length(bar_R2) != 1L || is.na(bar_R2) || !is.finite(bar_R2) || bar_R2 < 0 || bar_R2 >= 1) {
    stop("`bar_R2` must be a finite scalar in [0, 1).")
  }

  # ellipsoid
  if (!is.null(ellipsoid) && !is.list(ellipsoid)) {
    stop("`ellipsoid` must be the output of `precompute_ellipsoid()` or `NULL`.")
  }
}

#' Check the validity of the input
#'
#' @inheritParams sens_linear
#'
#' @noRd
check_input_linear <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid, c, c0, method) {
  check_input_ellipsoid(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid
  )
  dX <- ncol(as.matrix(X))

  # c
  if (!is.numeric(c) || !is.vector(c) || length(c) != dX || anyNA(c) || any(!is.finite(c))) {
    stop("`c` must be a finite numeric vector of length ncol(X).")
  }

  # c0
  if (!is.numeric(c0) || length(c0) != 1L || is.na(c0) || !is.finite(c0)) {
    stop("`c0` must be a finite numeric scalar.")
  }

  # method
  if (!(method %in% c("analytical", "gurobi"))) {
    stop('`method` must be either "analytical" or "gurobi".')
  }
}

# ----------------------------- Main Functions --------------------------------#

#' Sensitivity of linear functions of OLS coefficients
#'
#' `sens_linear_analytical()` calculates the closed-form sharp identified set for
#' linear functions of OLS coefficients \eqn{f_L(\beta) = c^{\intercal}\beta + c_0}.
#'
#' @inheritParams sens_linear
#'
#' @return A list of 4 objects where:
#' * `beta_lb` stores the coefficient vector beta attaining the lower bound,
#' * `beta_ub` stores the coefficient vector beta attaining the upper bound,
#' * `obj_lb` stores the sharp lower bound of \eqn{f_L(\beta)},
#' * `obj_ub` stores the sharp upper bound of \eqn{f_L(\beta)}.
#' @export
sens_linear_analytical <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid = NULL, c, c0) {
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }
  beta_med <- ellipsoid$beta_med
  Sigma_inv <- ellipsoid$Sigma_inv
  radius <- ellipsoid$radius

  if (all(c == 0)) {
    # c = 0: f_L = c0 is constant over the ellipsoid and every beta attains it
    return(list(
      beta_lb = as.numeric(beta_med),
      beta_ub = as.numeric(beta_med),
      obj_lb  = c0,
      obj_ub  = c0
    ))
  }

  scale <- as.numeric(sqrt(t(c) %*% Sigma_inv %*% c)) # scalar
  center_obj <- t(c) %*% beta_med + c0 # scalar
  direction <- (Sigma_inv %*% c) / scale # dX x 1 vector

  list(
    beta_lb = as.numeric(beta_med - radius * direction),
    beta_ub = as.numeric(beta_med + radius * direction),
    obj_lb  = as.numeric(center_obj - radius * scale),
    obj_ub  = as.numeric(center_obj + radius * scale)
  )
}

#' Sensitivity of linear functions of OLS coefficients
#'
#' `sens_linear_gurobi()` calculates the sharp identified set for linear functions
#' of OLS coefficients \eqn{f_L(\beta) = c^{\intercal}\beta + c_0} by gurobi.
#'
#' @inheritParams sens_linear
#'
#' @return A list of 4 objects where:
#' * `beta_lb` stores the coefficient vector beta attaining the lower bound,
#' * `beta_ub` stores the coefficient vector beta attaining the upper bound,
#' * `obj_lb` stores the sharp lower bound of \eqn{f_L(\beta)},
#' * `obj_ub` stores the sharp upper bound of \eqn{f_L(\beta)}.
#' @export
sens_linear_gurobi <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid = NULL, c, c0) {
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }
  result_lb <- solve_linear_gurobi(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    c = c, c0 = c0, sense = "min"
  )
  result_ub <- solve_linear_gurobi(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    c = c, c0 = c0, sense = "max"
  )

  list(
    beta_lb = result_lb$x,
    beta_ub = result_ub$x,
    obj_lb = result_lb$objval,
    obj_ub = result_ub$objval
  )
}

#' Sensitivity of linear functions of OLS coefficients
#'
#' `sens_linear()` calculates the sharp identified set for linear functions
#' of OLS coefficients \eqn{f_L(\beta) = c^{\intercal}\beta + c_0}.
#'
#' @param Y A \eqn{N \times 1} vector.
#' @param X A \eqn{N \times dX} dataframe/matrix/vector.
#' @param W1 A \eqn{N \times d1} dataframe/matrix/vector or `NULL`.
#' @param bar_rho A scalar in \eqn{[0, 1)}.
#' @param bar_R2 A scalar in \eqn{[0, 1)}.
#' @param ellipsoid A list computed from `precompute_ellipsoid()` (optional).
#' @param c A \eqn{dX \times 1} vector.
#' @param c0 A scalar.
#' @param method `"analytical"` or `"gurobi"`
#'
#' @return A list of 4 objects where:
#' * `beta_lb` stores the coefficient vector beta attaining the lower bound,
#' * `beta_ub` stores the coefficient vector beta attaining the upper bound,
#' * `obj_lb` stores the sharp lower bound of \eqn{f_L(\beta)},
#' * `obj_ub` stores the sharp upper bound of \eqn{f_L(\beta)}.
#' @export
sens_linear <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid = NULL, c, c0, method = "analytical") {
  check_input_linear(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    c = c, c0 = c0, method = method
  )
  if (method == "analytical") {
    return(sens_linear_analytical(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
      c = c, c0 = c0
    ))
  }
  if (method == "gurobi") {
    return(sens_linear_gurobi(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
      c = c, c0 = c0
    ))
  }
}
