# --------------------------- Helper Functions --------------------------------#

#' ajry application: run fixed effects OLS regression
#' @noRd
ajry_run_feols <- function(formula, data, cluster_var) {
  # Approach 1
  model <- stats::lm(formula, data)
  vcov_cl <- sandwich::vcovCL(model, cluster = data[[cluster_var]], type = "HC1")
  se_cl <- sqrt(diag(vcov_cl))
  regressors_vars <- names(stats::coef(model))[!is.na(stats::coef(model))]

  return(list(model = model, se_cl = se_cl, regressors_vars = regressors_vars))

  # Approach 2
  # fixest::feols(fml = formula, data = data, vcov = stats::as.formula(paste0("~", cluster_var)))
}

#' ajry application: extract Y, X, W1 from ajry data
#' @noRd
ajry_extract_vars <- function(data, outcome_var, treatment_vars, regressors_vars) {
  W1_vars <- setdiff(regressors_vars, c("(Intercept)", treatment_vars))

  keep_vars <- c(outcome_var, treatment_vars, W1_vars)
  data <- stats::na.omit(data[, keep_vars])

  list(
    Y = data[[outcome_var]],
    X = data[, treatment_vars],
    W1 = data[, W1_vars]
  )
}

#' ajry application: prepare ajry data
#' @noRd
ajry_setup_data <- function(acemoglu_et_al_2008, method) {
  if (method == "pols") {
    # data for Table 2 Column 1 (pooled OLS)
    fml_pols <- stats::as.formula(
      paste("fhpolrigaug ~ L1.fhpolrigaug + L1.lrgdpch +",
            paste(paste0("yr", 1:11), collapse = " + "))
    )
    list_pols <- ajry_run_feols(
      formula = fml_pols, data = acemoglu_et_al_2008, cluster_var = "code"
    )
    vars_pols <- ajry_extract_vars(
      data = acemoglu_et_al_2008,
      outcome_var = "fhpolrigaug",
      treatment_vars = c("L1.fhpolrigaug", "L1.lrgdpch"),
      regressors_vars = list_pols$regressors_vars
    )
    return(vars_pols)
  }

  if (method == "feols") {
    # data for Table 2 Column 2 (FE OLS)
    fml_feols <- stats::as.formula(
      paste("fhpolrigaug ~ L1.fhpolrigaug + L1.lrgdpch +",
            paste(paste0("yr", 1:11), collapse = " + "), "+",
            paste(paste0("cd", 1:210), collapse = " + "))
    )
    list_feols <- ajry_run_feols(
      formula = fml_feols, data = acemoglu_et_al_2008, cluster_var = "code"
    )
    vars_feols <- ajry_extract_vars(
      data = acemoglu_et_al_2008,
      outcome_var = "fhpolrigaug",
      treatment_vars = c("L1.fhpolrigaug", "L1.lrgdpch"),
      regressors_vars = list_feols$regressors_vars
    )
    return(vars_feols)
  }
}

# ----------------------------- Main Functions --------------------------------#

#' ajry application
#' @param acemoglu_et_al_2008 A dataframe of replication data of Acemoglu et al (2008).
#' @export
ajry_ovb_func <- function(acemoglu_et_al_2008) {
  result <- list(pols = NULL, feol = NULL)
  for (method in c("pols", "feols")) {
    # Prepare ajry data
    ajry_data <- ajry_setup_data(
      acemoglu_et_al_2008 = acemoglu_et_al_2008, method = method
    )
    Y <- ajry_data$Y
    X <- ajry_data$X
    W1 <- ajry_data$W1

    result[[method]]$sens_linear <- sens_linear(
      Y = Y, X = X, W1 = W1, bar_rho = 0.3, bar_R2 = 0.3, ellipsoid = NULL,
      c = c(1, 1), c0 = 0, method = "analytical"
    )

    result[[method]]$sens_linear_frac <- sens_linear_frac(
      Y = Y, X = X, W1 = W1, bar_rho = 0.3, bar_R2 = 0.3, ellipsoid = NULL,
      c = c(0, 1), c0 = 0, d = c(-1, 0), d0 = 1
    )

    f <- function(beta) sum(beta)
    grad_f <- function(beta) c(1,1)
    result[[method]]$sens_general_linear <- sens_general(
      Y = Y, X = X, W1 = W1, bar_rho = 0.3, bar_R2 = 0.3, ellipsoid = NULL,
      f = f, grad_f = grad_f, beta_init = c(0.7, 0.1), step_size = 0.01,
      max_iter = 5000, tol = 1e-8, verbose = 0
    )

    f <- function(beta) beta[2]/(1 - beta[1])
    grad_f <- function(beta) c(beta[2]/((1 - beta[1])^2), 1/(1 - beta[1]))
    result[[method]]$sens_general_linear_frac <- sens_general(
      Y = Y, X = X, W1 = W1, bar_rho = 0.3, bar_R2 = 0.3, ellipsoid = NULL,
      f = f, grad_f = grad_f, beta_init = c(0.7, 0.1), step_size = 0.01,
      max_iter = 5000, tol = 1e-8, verbose = 0
    )
  }
  return(result)
}
