# code to prepare `acemoglu_et_al_2008` dataset
# Output: data/acemoglu_et_al_2008.rda

library(readxl)
library(tidyverse)
library(sandwich)

# Import the data
ajry <- as.data.frame(read_excel("data-raw/acemoglu_et_al_2008/data_5_year_panel.xlsx"))

# Sort the data by country and year
sort_vars <- c("code_numeric", "year_numeric")
ajry <- arrange(ajry, across(all_of(sort_vars)))

# Clean the data
max_lag <- 5
vars_lag <- c("fhpolrigaug","lrgdpch")
lag_funs <- setNames(
  lapply(seq_len(max_lag), function(k) {
    force(k)
    function(x) dplyr::lag(x, k)
  }),
  paste0("L", seq_len(max_lag))
)

# Generate lags
ajry <- ajry %>%
  arrange(.data[["code_numeric"]], .data[["year_numeric"]]) %>%
  group_by(.data[["code_numeric"]]) %>%
  mutate(
    across(
      all_of(vars_lag),
      c(lag_funs, list(D1 = ~ .x - dplyr::lag(.x, 1))),
      .names = "{.fn}.{.col}"
    )
  ) %>%
  ungroup()

# Generate differences
for (var in vars_lag) {
  ajry[[paste0("LD.", var)]] <- ajry[[paste0("L1.", var)]] - ajry[[paste0("L2.", var)]]
}

# Generate year and country dummies
ajry$code_numeric_factor <- factor(ajry$code_numeric)
cd_mat <- model.matrix(stats::as.formula("~ code_numeric_factor - 1"), data = ajry)
colnames(cd_mat) <- paste0("cd", seq_len(ncol(cd_mat)))

ajry$year_numeric_factor <- factor(ajry$year_numeric)
yr_mat <- model.matrix(stats::as.formula("~ year_numeric_factor - 1"), data = ajry)
colnames(yr_mat) <- paste0("yr", seq_len(ncol(yr_mat)))

ajry <- cbind(ajry, cd_mat, yr_mat)
ajry <- filter(ajry, sample == 1)
acemoglu_et_al_2008 <- ajry

usethis::use_data(acemoglu_et_al_2008, overwrite = TRUE, compress = "xz")
