# --------------------------- Helper Functions --------------------------------#

#' Set ggplot2 theme
#'
#' @return A ggplot2 theme object.
#' @noRd
base_theme <- function() {
  theme(
    legend.position = "bottom",
    text = element_text(size = 15, family = "CMU Serif"),
    plot.title = element_text(hjust = 0.5, size = 15),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

#' Save ggplot object to PDF
#'
#' @param plot_obj A ggplot object.
#' @param filename Output filename.
#'
#' @return A saved PDF file.
#' @noRd
save_plot <- function(plot_obj, filename) {
  ggsave(
    filename, plot = plot_obj, device = CairoPDF, width = 20, height = 15, units = "cm"
  )
}

# ----------------------------- Main Functions --------------------------------#

#' Draw ellipsoid with two treatments (dX = 2)
#'
#' @param Y A \eqn{N \times 1} vector.
#' @param X A \eqn{N \times d_X} dataframe/matrix/vector.
#' @param W1 A \eqn{N \times d_1} dataframe/matrix/vector or `NULL`.
#' @param bar_rho A scalar in \eqn{[0, 1)}.
#' @param bar_R2 A scalar in \eqn{[0, 1)}.
#' @param ellipsoid A list computed from `precompute_ellipsoid()` (optional).
#' @param func Function type (`"linear"` or `"linear_frac"`).
#' @param func_list A list:
#' * if `func = "linear`, `func_list` contains `c` and `c0`.
#' * if `func = "linear_frac`, `func_list` contains `c`, `c0`, `d`, and `d0`.
#' @param npoints Number of boundary points (default 100).
#' @param output_file Output filename.
#'
#' @return A ggplot object.
#' @export
plot_ellipsoid <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid,
                           func, func_list, npoints = 100, output_file) {
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }
  beta_med <- ellipsoid$beta_med
  Sigma_nhalf <- ellipsoid$Sigma_nhalf
  radius <- ellipsoid$radius

  theta <- seq(0, 2 * pi, length.out = npoints)
  U <- rbind(cos(theta), sin(theta)) # 2 x npoints matrix

  # Boundary points
  E <- matrix(beta_med, nrow = 2, ncol = npoints) + Sigma_nhalf %*% U * radius
  E <- as.data.frame(t(E))
  colnames(E) <- c("beta1", "beta2")

  p <- ggplot(E, aes(x = .data$beta1, y = .data$beta2)) +
    geom_polygon(fill = "skyblue", alpha = 0.3, color = "blue", linewidth = 1) +
    annotate("point", x = beta_med[1], y = beta_med[2], size = 3) +
    labs(x = expression(beta[1]), y = expression(beta[2])) +
    # coord_equal() +
    base_theme()

  # Add tangent lines for the linear objective
  if (func == "linear") {
    result_linear <- sens_linear(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
      c = func_list$c, c0 = func_list$c0, method = "analytical"
    )
    obj_lb <- result_linear$obj_lb
    obj_ub <- result_linear$obj_ub
    c1 <- func_list$c[1]
    c2 <- func_list$c[2]
    c0 <- func_list$c0
    if (c2 != 0) {
      p <- p + geom_abline(
        intercept = (obj_lb - c0)/c2,
        slope = -c1/c2,
        color = "red",
        linewidth = 0.8
      ) +
        geom_abline(
          intercept = (obj_ub - c0)/c2,
          slope = -c1/c2,
          color = "red",
          linewidth = 0.8
        )
    }
    if (c2 == 0) {
      p <- p + geom_vline(
        xintercept = (obj_lb - c0)/c1,
        color = "red",
        linewidth = 0.8
      ) +
        geom_vline(
          xintercept = (obj_ub - c0)/c1,
          color = "red",
          linewidth = 0.8
        )
    }
  }

  # Add tangent lines for the linear fractional objective
  if (func == "linear_frac") {
    result_linear_frac <- sens_linear_frac(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2, ellipsoid = ellipsoid,
      c = func_list$c, c0 = func_list$c0, d = func_list$d, d0 = func_list$d0
    )
    obj_lb <- result_linear_frac$obj_lb
    obj_ub <- result_linear_frac$obj_ub
    c1 <- func_list$c[1]
    c2 <- func_list$c[2]
    d1 <- func_list$d[1]
    d2 <- func_list$d[2]
    c0 <- func_list$c0
    d0 <- func_list$d0
    if (c2 - d2 * obj_lb != 0) {
      p <- p + geom_abline(
        intercept = (d0 * obj_lb - c0) / (c2 - d2 * obj_lb),
        slope = (d1 * obj_lb - c1)/(c2 - d2 * obj_lb),
        color = "red",
        linewidth = 0.8
      )
    }
    if (c2 - d2 * obj_lb == 0) {
      p <- p + geom_vline(
        xintercept = (c0 - d0 * obj_lb)/(d1 * obj_lb - c1),
        color = "red",
        linewidth = 0.8
      )
    }
    if (c2 - d2 * obj_ub != 0) {
      p <- p + geom_abline(
        intercept = (d0 * obj_ub - c0) / (c2 - d2 * obj_ub),
        slope = (d1 * obj_ub - c1)/(c2 - d2 * obj_ub),
        color = "red",
        linewidth = 0.8
      )
    }
    if (c2 - d2 * obj_ub == 0) {
      p <- p + geom_vline(
        xintercept = (c0 - d0 * obj_ub)/(d1 * obj_ub - c1),
        color = "red",
        linewidth = 0.8
      )
    }
  }

  save_plot(p, output_file)
  return(p)
}

#' Draw breakdown frontier
#'
#' @param Y A \eqn{N \times 1} vector.
#' @param X A \eqn{N \times d_X} dataframe/matrix/vector.
#' @param W1 A \eqn{N \times d_1} dataframe/matrix/vector or `NULL`.
#' @param bar_rho A scalar in \eqn{[0, 1)}.
#' @param bar_R2 A scalar in \eqn{[0, 1)}.
#' @param ellipsoid A list computed from `precompute_ellipsoid()` (optional).
#' @param j An integer in \eqn{[1, d_X]}.
#' @param npoints Number of boundary points (default 100).
#' @param output_file Output filename.
#'
#' @return A ggplot object.
#' @export
plot_breakdown_frontier <- function(Y, X, W1, bar_rho, bar_R2, ellipsoid,
                                    j, npoints = 100, output_file) {
  if (is.null(ellipsoid)) {
    ellipsoid <- precompute_ellipsoid(
      Y = Y, X = X, W1 = W1, bar_rho = bar_rho, bar_R2 = bar_R2
    )
  }
  beta_med <- ellipsoid$beta_med
  Sigma_inv <- ellipsoid$Sigma_inv
  sigma2_Y_perp_XW1 <- ellipsoid$sigma2_Y_perp_XW1

  bar_rho_grid <- seq(0, 1, length.out = npoints)
  bf_cons <- (beta_med[j]^2)/(sigma2_Y_perp_XW1 * Sigma_inv[j, j])
  bar_R2_grid <- ifelse(bar_rho_grid == 0, Inf, bf_cons * (1 - bar_rho_grid^2) / bar_rho_grid^2)
  bar_R2_grid <- pmin(pmax(bar_R2_grid, 0), 1)

  df <- data.frame(bar_rho = bar_rho_grid, bar_R2 = bar_R2_grid)
  p <- ggplot(df, aes(x = bar_rho, y = bar_R2)) +
    # robust region
    geom_ribbon(aes(ymin = 0, ymax = bar_R2), fill = "steelblue", alpha = 0.3) +
    # breakdown frontier
    geom_line(color = "blue", linewidth = 1) +
    # coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = expression(bar(rho)),
         y = expression(bar(R)[Y]^2)) +
    base_theme()

  save_plot(p, output_file)
  return(p)
}
