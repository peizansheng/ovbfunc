# --------------------------- Helper Functions --------------------------------#

#' Sensitivity of linear fractional functions of OLS coefficients
#'
#' `solve_linear_frac_gurobi()` solves the optimization problem for linear fractional
#' functions of OLS coefficients
#' \eqn{f_\text{LF}(\beta) = \frac{c^{\intercal}\beta + c_0}{d^{\intercal}\beta + d_0}} by gurobi.
#'
#' @inheritParams sens_linear_frac_gurobi
#' @param sense `"min"` or `"max"`
#'
#' @return A gurobi object.
#' @noRd
solve_linear_frac_gurobi <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid = NULL,
                                     c, c0, d, d0, sense) {
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }
  beta_med <- ellipsoid$beta_med
  Sigma_half <- ellipsoid$Sigma_half
  radius <- ellipsoid$radius

  dX <- length(beta_med)

  # Check denominator sign
  sens_denom <- sens_linear_analytical(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    c = d, c0 = d0
  )
  if (sens_denom$obj_lb > 0) {
    sign <- "positive"
  } else if (sens_denom$obj_ub < 0) {
    sign <- "negative"
  } else {
    sign <- "not sure"
  }

  if (sign == "not sure"){
    stop(paste("The denominator of linear fractional objective changes sign over the feasible region."))
  }
  # Negate c, c0, d, d0 if d'beta + d0 < 0,
  # (c'beta + c0) / (d'beta + d0) = (-c'beta - c0) / (-d'beta - d0) then the denominator is positive
  if (sign == "negative"){
    c <- -c
    c0 <- -c0
    d <- -d
    d0 <- -d0
  }

  # Build Gurobi model
  model <- list()
  model$obj <- c(c, c0, rep(0, dX)) # (phi, t, w)
  model$modelsense <- sense
  model$lb <- c(rep(-Inf, dX), 0, rep(-Inf, dX))
  model$ub <- rep(Inf, 2 * dX + 1)

  # Linear constraint: Sigma^{1/2} phi - Sigma^{1/2} beta_med t + w = 0, d' phi + d0 * t = 1
  a1 <- cbind(Sigma_half, -as.numeric(Sigma_half %*% beta_med), diag(dX)) # dX x (2 * dX + 1) matrix
  a2 <- matrix(c(d, d0, rep(0, dX)), nrow = 1) # 1 x (2 * dX + 1) matrix
  model$A <- rbind(a1, a2)
  model$rhs <- c(rep(0, dX), 1)
  model$sense <- c(rep("=", dX + 1))

  # Quadratic constraint w' w <= radius^2 t^2
  Qc <- Matrix::spMatrix(
    nrow = 2 * dX + 1,
    ncol = 2 * dX + 1,
    i = c((dX + 2):(2 * dX + 1), dX + 1),
    j = c((dX + 2):(2 * dX + 1), dX + 1),
    x = c(rep(1, dX), -radius^2)
  )

  qc <- list()
  qc$Qc <- Qc
  qc$q <- rep(0, 2 * dX + 1)
  qc$rhs <- 0
  qc$sense <- "<"
  model$quadcon <- list(qc)

  # Solve gurobi model
  result <- solve_gurobi_model(model)

  return(result)
}

# ----------------------------- Main Functions --------------------------------#

#' Sensitivity of linear fractional functions of OLS coefficients
#'
#' `sens_linear_frac_gurobi()` calculates the extremal values for linear fractional
#' functions of OLS coefficients
#' \eqn{f_\text{LF}(\beta) = \frac{c^{\intercal}\beta + c_0}{d^{\intercal}\beta + d_0}} by gurobi.
#'
#' @param Y A \eqn{N \times 1} vector.
#' @param X A \eqn{N \times dX} dataframe/matrix.
#' @param W1 A \eqn{N \times d1} dataframe/matrix or `NULL`.
#' @param bar_rho A scalar in \eqn{[0, 1)}.
#' @param bar_R2 A scalar in \eqn{[0, 1)}.
#' @param ellipsoid Optional output of `precompute_ellipsoid()`.
#' @param c A \eqn{dX \times 1} vector.
#' @param c0 A scalar.
#' @param d A \eqn{dX \times 1} vector.
#' @param d0 A scalar.
#'
#' @return A list of 4 objects where:
#' * `beta_lb` stores the coefficient vector beta attaining the lower bound,
#' * `beta_ub` stores the coefficient vector beta attaining the upper bound,
#' * `obj_lb` stores the lower bound of \eqn{f_{LF}(\beta)},
#' * `obj_ub` stores the upper bound of \eqn{f_{LF}(\beta)}.
#' @export
sens_linear_frac_gurobi <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid = NULL,
                                    c, c0, d, d0) {
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }
  dX <- length(ellipsoid$beta_med)
  result_lb <- solve_linear_frac_gurobi(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    c = c, c0 = c0, d = d, d0 = d0, sense = "min"
  )
  result_ub <- solve_linear_frac_gurobi(
    Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
    c = c, c0 = c0, d = d, d0 = d0, sense = "max"
  )

  beta_lb <- result_lb$x[1:dX] / result_lb$x[dX + 1]
  beta_ub <- result_ub$x[1:dX] / result_ub$x[dX + 1]
  obj_lb <- (sum(c * beta_lb) + c0) / (sum(d * beta_lb) + d0)
  obj_ub <- (sum(c * beta_ub) + c0) / (sum(d * beta_ub) + d0)

  list(
    beta_lb = beta_lb,
    beta_ub = beta_ub,
    obj_lb = obj_lb,
    obj_ub = obj_ub
    # obj_lb = result_lb$objval,
    # obj_ub = result_ub$objval
  )
}
