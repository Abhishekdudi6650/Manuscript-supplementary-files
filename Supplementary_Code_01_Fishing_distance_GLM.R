# =============================================================================
# SUPPLEMENTARY CODE 01
# Fishing-distance GLM analysis, model averaging, diagnostics and robustness
# =============================================================================
#
# Reproducibility workflow
# ------------------------
# Input:
#   model_data_recoded_full.csv
#
# Run this file from top to bottom in a fresh R session.
#
# Part A:
#   - recreates the analysis dataset;
#   - compares the global Gaussian and Gamma models;
#   - evaluates all 65,536 additive subsets of the 16-predictor Gaussian model;
#   - retains models with strict Delta AICc < 2;
#   - performs full model averaging;
#   - produces GVIF, residual, assumption and influence diagnostics;
#   - performs the Cook's-distance sensitivity analysis by refitting and
#     reweighting the same originally supported model formulas.
#
# Part B:
#   - verifies the supported-model set;
#   - computes HC3 model-averaged uncertainty;
#   - performs 2,000 non-parametric pairs-bootstrap repetitions;
#   - exports percentile confidence intervals and model-weight stability.
#
# Main output folders:
#   Fishing_distance_diagnostics_DeltaAICc_lt2/
#   Fishing_distance_HC3_bootstrap_robustness/
#
# The script does not automatically remove influential observations from the
# primary analysis. Cook's-distance exclusions are used only for sensitivity
# analysis.
# =============================================================================

# PART A. MODEL SELECTION, MODEL AVERAGING AND DIAGNOSTICS
# =============================================================================
#
# This script performs the full analysis and diagnostic workflow:
#
#   1. Recreates the analysis variables and one complete-case dataset.
#   2. Fits the global Gaussian GLM with identity link.
#   3. Fits every additive subset of the global model.
#   4. Retains ONLY models with Delta AICc < 2 (strictly below 2).
#   5. Performs AICc-weighted full model averaging.
#   6. Checks:
#        - data completeness, ranges and category counts;
#        - Spearman correlations among continuous variables;
#        - model convergence and rank deficiency;
#        - GVIF / adjusted GVIF for the global and supported models;
#        - residual normality (plots and Shapiro-Wilk test);
#        - heteroscedasticity (plots and Breusch-Pagan test);
#        - functional-form misspecification (Ramsey RESET test);
#        - fitted values below zero;
#        - studentised residuals and potential outliers;
#        - leverage, Cook's distance, DFFITS and DFBETAs;
#        - residuals from the model-averaged predictions;
#        - Gaussian versus Gamma global-model AICc as an additional
#          distributional sensitivity check (when the response is positive).
#   7. Takes the union of observations with Cook's distance > 4/n in ANY
#      supported model, excludes those observations for a sensitivity check,
#      and refits/reweights the same 10 originally supported model formulas.
#   8. Compares original and Cook's-distance sensitivity estimates.
#
# Required input file:
#   model_data_recoded_full.csv
#
# Expected packages:
# install.packages(c(
#   "tidyverse", "MuMIn", "car", "lmtest", "broom"
# ))
# ============================================================


# ------------------------------------------------------------
# 0. PACKAGES AND FIXED ANALYSIS SETTINGS
# ------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "MuMIn",
  "car",
  "lmtest",
  "broom"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(MuMIn)
  library(car)
  library(lmtest)
  library(broom)
})

# Required by MuMIn::dredge so every model uses identical observations.
options(na.action = "na.fail")

# Supported-model rule used in the manuscript: strict Delta AICc < 2.
DELTA_AICC_CUTOFF <- 2

# Diagnostic thresholds.
COOKS_MULTIPLIER <- 4
STUDENTISED_RESIDUAL_THRESHOLD <- 3
LEVERAGE_MULTIPLIER <- 2
GVIF_ADJUSTED_FLAG_1 <- 1.5
GVIF_ADJUSTED_FLAG_2 <- 2

# Change to TRUE only if you want progress printed for all 65,536 models.
SHOW_DREDGE_PROGRESS <- FALSE

output_dir <- "Fishing_distance_diagnostics_DeltaAICc_lt2"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Some fitted-model objects do not store a named "converged" component.
# Return NA in that case rather than attempting unsafe $ or [[ access.
model_convergence_status <- function(model) {
  model_names <- names(model)

  if (
    is.list(model) &&
    !is.null(model_names) &&
    "converged" %in% model_names
  ) {
    return(isTRUE(model[["converged"]]))
  }

  NA
}


# ------------------------------------------------------------
# 1. IMPORT DATA
# ------------------------------------------------------------

data_file <- "model_data_recoded_full.csv"

if (!file.exists(data_file)) {
  message(
    "model_data_recoded_full.csv was not found in the working directory. ",
    "Please select it manually."
  )
  data_file <- file.choose()
}

fishers_raw <- readr::read_csv(
  data_file,
  show_col_types = FALSE
) %>%
  mutate(.observation_id = row_number())


# ------------------------------------------------------------
# 2. RECREATE BYCATCH AND GEAR VARIABLES
# ------------------------------------------------------------

fishers_clean <- fishers_raw %>%
  mutate(
    bycatch_clean = case_when(
      bycatch_raw == "Never" ~ "Never",
      bycatch_raw %in% c("Rarely", "Occasionally", "Sometimes") ~
        "Infrequent",
      bycatch_raw == "On every trip" ~ "Every trip",
      TRUE ~ NA_character_
    ),
    bycatch_clean = factor(
      bycatch_clean,
      levels = c("Never", "Infrequent", "Every trip"),
      ordered = TRUE
    ),
    Site = relevel(factor(Site), ref = "A"),
    AIS = relevel(factor(AIS), ref = "No"),
    VMS = relevel(factor(VMS), ref = "No"),
    navigation.apps = relevel(factor(navigation.apps), ref = "No"),
    target.sharks = relevel(factor(target.sharks), ref = "No"),
    subsidy = relevel(factor(subsidy), ref = "No"),
    Top.species.group = relevel(
      factor(Top.species.group),
      ref = "Cephalopods"
    ),
    communication = relevel(
      factor(communication),
      ref = "Every trip"
    ),
    reg_awareness = relevel(
      factor(reg_awareness),
      ref = "Fully aware"
    ),
    gear_type_clean = factor(gear_type_clean),
    uses_longline = relevel(
      factor(
        if_else(
          str_detect(
            gear_type_clean,
            regex("Longline", ignore_case = TRUE)
          ),
          "Yes",
          "No"
        )
      ),
      ref = "No"
    ),
    uses_gillnet = relevel(
      factor(
        if_else(
          str_detect(
            gear_type_clean,
            regex("Gillnet", ignore_case = TRUE)
          ),
          "Yes",
          "No"
        )
      ),
      ref = "No"
    ),
    uses_handline = relevel(
      factor(
        if_else(
          str_detect(
            gear_type_clean,
            regex("Handline", ignore_case = TRUE)
          ),
          "Yes",
          "No"
        )
      ),
      ref = "No"
    )
  )


# ------------------------------------------------------------
# 3. CREATE ONE COMPLETE-CASE DATASET
# ------------------------------------------------------------

analysis_variables <- c(
  ".observation_id",
  "mean_distance",
  "Site",
  "Fishing.experience",
  "vessel_length_m",
  "cost_per_trip_usd",
  "AIS",
  "VMS",
  "navigation.apps",
  "target.sharks",
  "subsidy",
  "Top.species.group",
  "communication",
  "reg_awareness",
  "uses_gillnet",
  "uses_handline",
  "uses_longline",
  "bycatch_clean"
)

analysis_data <- fishers_clean %>%
  select(all_of(analysis_variables)) %>%
  drop_na() %>%
  droplevels()

if (nrow(analysis_data) == 0) {
  stop("No complete observations remain for the analysis.")
}

if (anyDuplicated(analysis_data$.observation_id)) {
  stop("Observation identifiers are not unique.")
}

factor_level_counts <- vapply(
  analysis_data %>% select(where(is.factor)),
  nlevels,
  integer(1)
)

if (any(factor_level_counts < 2)) {
  stop(
    "These factors have fewer than two observed levels: ",
    paste(names(factor_level_counts)[factor_level_counts < 2], collapse = ", ")
  )
}

cat("Rows in imported data:", nrow(fishers_raw), "\n")
cat("Complete cases used in every original candidate model:",
    nrow(analysis_data), "\n")


# ------------------------------------------------------------
# 4. DESCRIPTIVE DATA CHECKS
# ------------------------------------------------------------

numeric_variables <- analysis_data %>%
  select(where(is.numeric)) %>%
  select(-.observation_id)

numeric_summary <- purrr::imap_dfr(
  numeric_variables,
  function(x, variable_name) {
    tibble(
      variable = variable_name,
      n = sum(!is.na(x)),
      missing = sum(is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      minimum = min(x, na.rm = TRUE),
      maximum = max(x, na.rm = TRUE),
      skewness = {
        centred <- x - mean(x, na.rm = TRUE)
        mean(centred^3, na.rm = TRUE) /
          (sd(x, na.rm = TRUE)^3)
      }
    )
  }
)

factor_counts <- analysis_data %>%
  select(where(is.factor)) %>%
  purrr::imap_dfr(
    function(x, variable_name) {
      tibble(
        variable = variable_name,
        level = as.character(x)
      ) %>%
        count(variable, level, name = "n") %>%
        mutate(proportion = n / sum(n))
    }
  )

readr::write_csv(
  numeric_summary,
  file.path(output_dir, "01_numeric_variable_summary.csv")
)

readr::write_csv(
  factor_counts,
  file.path(output_dir, "02_factor_category_counts.csv")
)

# Spearman correlations and pairwise p-values among continuous variables.
correlation_variables <- analysis_data %>%
  select(
    mean_distance,
    Fishing.experience,
    vessel_length_m,
    cost_per_trip_usd
  )

spearman_correlation <- cor(
  correlation_variables,
  method = "spearman",
  use = "pairwise.complete.obs"
)

spearman_p_values <- matrix(
  NA_real_,
  nrow = ncol(correlation_variables),
  ncol = ncol(correlation_variables),
  dimnames = list(
    names(correlation_variables),
    names(correlation_variables)
  )
)

for (i in seq_len(ncol(correlation_variables))) {
  for (j in seq_len(ncol(correlation_variables))) {
    spearman_p_values[i, j] <- suppressWarnings(
      cor.test(
        correlation_variables[[i]],
        correlation_variables[[j]],
        method = "spearman",
        exact = FALSE
      )$p.value
    )
  }
}

readr::write_csv(
  as.data.frame(spearman_correlation) %>%
    rownames_to_column("variable"),
  file.path(output_dir, "03_spearman_correlations.csv")
)

readr::write_csv(
  as.data.frame(spearman_p_values) %>%
    rownames_to_column("variable"),
  file.path(output_dir, "04_spearman_correlation_p_values.csv")
)

png(
  file.path(output_dir, "05_response_distribution.png"),
  width = 1800,
  height = 900,
  res = 180
)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
hist(
  analysis_data$mean_distance,
  breaks = "FD",
  main = "Distribution of average fishing distance",
  xlab = "Average fishing distance (km)",
  col = "grey85",
  border = "grey30"
)
plot(
  density(analysis_data$mean_distance),
  main = "Density of average fishing distance",
  xlab = "Average fishing distance (km)",
  lwd = 2
)
rug(analysis_data$mean_distance)
dev.off()


# ------------------------------------------------------------
# 5. GLOBAL GAUSSIAN MODEL
# ------------------------------------------------------------

global_formula <- mean_distance ~
  Site +
  Fishing.experience +
  vessel_length_m +
  cost_per_trip_usd +
  AIS +
  VMS +
  navigation.apps +
  target.sharks +
  subsidy +
  Top.species.group +
  communication +
  reg_awareness +
  uses_gillnet +
  uses_handline +
  uses_longline +
  bycatch_clean

global_gaussian <- glm(
  formula = global_formula,
  data = analysis_data,
  family = gaussian(link = "identity"),
  na.action = na.fail
)

if (identical(model_convergence_status(global_gaussian), FALSE)) {
  stop("The global Gaussian model did not converge.")
}

global_matrix <- model.matrix(global_gaussian)
global_rank_deficient <- qr(global_matrix)$rank < ncol(global_matrix)

if (global_rank_deficient) {
  stop(
    "The global Gaussian model matrix is rank deficient. Resolve redundant ",
    "terms or factor levels before continuing."
  )
}


# ------------------------------------------------------------
# 6. OPTIONAL GLOBAL GAUSSIAN-VERSUS-GAMMA CHECK
# ------------------------------------------------------------

family_comparison <- tibble(
  family = "Gaussian",
  link = "identity",
  converged = model_convergence_status(global_gaussian),
  logLik = as.numeric(logLik(global_gaussian)),
  df = attr(logLik(global_gaussian), "df"),
  AICc = MuMIn::AICc(global_gaussian),
  minimum_fitted_value = min(fitted(global_gaussian)),
  negative_fitted_values = sum(fitted(global_gaussian) < 0)
)

global_gamma <- NULL

if (all(analysis_data$mean_distance > 0)) {
  global_gamma <- tryCatch(
    glm(
      formula = global_formula,
      data = analysis_data,
      family = Gamma(link = "log"),
      na.action = na.fail
    ),
    error = function(e) NULL
  )

  if (!is.null(global_gamma)) {
    family_comparison <- bind_rows(
      family_comparison,
      tibble(
        family = "Gamma",
        link = "log",
        converged = model_convergence_status(global_gamma),
        logLik = as.numeric(logLik(global_gamma)),
        df = attr(logLik(global_gamma), "df"),
        AICc = MuMIn::AICc(global_gamma),
        minimum_fitted_value = min(fitted(global_gamma)),
        negative_fitted_values = sum(fitted(global_gamma) < 0)
      )
    )
  }
}

family_comparison <- family_comparison %>%
  arrange(AICc) %>%
  mutate(delta_AICc = AICc - min(AICc))

readr::write_csv(
  family_comparison,
  file.path(output_dir, "06_global_family_comparison.csv")
)


# ------------------------------------------------------------
# 7. ALL-SUBSETS AICc SELECTION: STRICT Delta AICc < 2
# ------------------------------------------------------------

candidate_table <- MuMIn::dredge(
  global_gaussian,
  rank = "AICc",
  trace = SHOW_DREDGE_PROGRESS
)

candidate_table_df <- as.data.frame(candidate_table) %>%
  rownames_to_column("model_id")

supported_table_df <- candidate_table_df %>%
  filter(delta < DELTA_AICC_CUTOFF)

if (nrow(supported_table_df) < 2) {
  stop(
    "Only ", nrow(supported_table_df),
    " model is strictly below Delta AICc 2. Multiple-model averaging cannot ",
    "be performed. Inspect the selection table before proceeding."
  )
}

supported_models <- MuMIn::get.models(
  candidate_table,
  subset = delta < DELTA_AICC_CUTOFF
)

supported_model_ids <- supported_table_df$model_id
supported_model_labels <- sprintf(
  "supported_%03d_model_%s",
  seq_along(supported_models),
  supported_model_ids
)
names(supported_models) <- supported_model_labels

averaged_model <- MuMIn::model.avg(
  supported_models,
  revised.var = TRUE
)

readr::write_csv(
  candidate_table_df,
  file.path(output_dir, "07_all_candidate_models_AICc.csv")
)

readr::write_csv(
  supported_table_df,
  file.path(output_dir, "08_supported_models_DeltaAICc_lt2.csv")
)

component_table_renormalised <- as.data.frame(
  averaged_model$msTable
) %>%
  rownames_to_column("component_model")

readr::write_csv(
  component_table_renormalised,
  file.path(output_dir, "09_supported_models_renormalised_weights.csv")
)


# ------------------------------------------------------------
# 8. FUNCTIONS FOR MODEL-AVERAGED COEFFICIENTS AND GVIF
# ------------------------------------------------------------

# A Gaussian GLM with an identity link is equivalent to an ordinary linear
# model for these diagnostics. Rebuilding an lm object avoids method-dispatch
# failures in car, influence and lmtest functions that can occur for glm
# objects in some package-version combinations.
diagnostic_lm_from_glm <- function(model) {
  if (!inherits(model, "glm")) {
    stop(
      "Diagnostic input is not a fitted glm object; received class: ",
      paste(class(model), collapse = ", ")
    )
  }

  stats::lm(
    formula = stats::formula(model),
    data = stats::model.frame(model),
    na.action = na.fail
  )
}

extract_full_average <- function(model_average) {
  average_summary <- summary(model_average)

  coefficient_table <- as.data.frame(
    average_summary$coefmat.full
  ) %>%
    rownames_to_column("term")

  estimate_column <- grep(
    "^Estimate$",
    names(coefficient_table),
    value = TRUE
  )[1]

  se_column <- intersect(
    c("Adjusted SE", "Std. Error"),
    names(coefficient_table)
  )[1]

  p_column <- intersect(
    c("Pr(>|z|)", "Pr(>|t|)"),
    names(coefficient_table)
  )[1]

  if (is.na(estimate_column) || is.na(se_column)) {
    stop("Could not identify model-averaged estimate or SE columns.")
  }

  confidence_intervals <- as.data.frame(
    confint(model_average, full = TRUE, level = 0.95)
  ) %>%
    rownames_to_column("term")

  names(confidence_intervals)[2:3] <- c("conf.low", "conf.high")

  coefficient_table %>%
    transmute(
      term = term,
      estimate = .data[[estimate_column]],
      unconditional_se = .data[[se_column]],
      p.value = if (!is.na(p_column)) .data[[p_column]] else NA_real_
    ) %>%
    left_join(confidence_intervals, by = "term")
}

extract_gvif <- function(model, model_label) {
  if (!inherits(model, "glm")) {
    stop(
      "GVIF input '", model_label,
      "' is not a fitted glm object; class is: ",
      paste(class(model), collapse = ", ")
    )
  }

  diagnostic_lm <- diagnostic_lm_from_glm(model)
  predictor_terms <- attr(terms(model), "term.labels")

  if (length(predictor_terms) < 2) {
    return(
      tibble(
        model = model_label,
        term = predictor_terms,
        GVIF_or_VIF = NA_real_,
        df = NA_real_,
        adjusted_GVIF = NA_real_,
        adjusted_GVIF_above_1_5 = NA,
        adjusted_GVIF_above_2 = NA,
        calculation_status = "Not calculated: model has fewer than two terms"
      )
    )
  }

  vif_result <- tryCatch(
    car::vif(diagnostic_lm),
    error = function(e) e
  )

  if (inherits(vif_result, "error")) {
    return(
      tibble(
        model = model_label,
        term = NA_character_,
        GVIF_or_VIF = NA_real_,
        df = NA_real_,
        adjusted_GVIF = NA_real_,
        adjusted_GVIF_above_1_5 = NA,
        adjusted_GVIF_above_2 = NA,
        calculation_status = paste("GVIF/VIF failed:", vif_result$message)
      )
    )
  }

  if (is.matrix(vif_result)) {
    result <- tibble(
      model = model_label,
      term = rownames(vif_result),
      GVIF_or_VIF = as.numeric(vif_result[, 1]),
      df = as.numeric(vif_result[, 2]),
      adjusted_GVIF = if (ncol(vif_result) >= 3) {
        as.numeric(vif_result[, 3])
      } else {
        as.numeric(vif_result[, 1])^(
          1 / (2 * as.numeric(vif_result[, 2]))
        )
      }
    )
  } else {
    # When every term has 1 df, car::vif returns ordinary VIF values.
    result <- tibble(
      model = model_label,
      term = names(vif_result),
      GVIF_or_VIF = as.numeric(vif_result),
      df = 1,
      adjusted_GVIF = sqrt(as.numeric(vif_result))
    )
  }

  result %>%
    mutate(
      adjusted_GVIF_above_1_5 = adjusted_GVIF >= GVIF_ADJUSTED_FLAG_1,
      adjusted_GVIF_above_2 = adjusted_GVIF >= GVIF_ADJUSTED_FLAG_2,
      calculation_status = "Calculated"
    )
}


# ------------------------------------------------------------
# 9. MODEL-AVERAGED COEFFICIENTS AND TERM WEIGHTS
# ------------------------------------------------------------

original_averaged_coefficients <- extract_full_average(averaged_model)

term_weight_table <- tibble(
  term = names(averaged_model$sw),
  sum_of_model_weights = as.numeric(averaged_model$sw)
)

readr::write_csv(
  original_averaged_coefficients,
  file.path(output_dir, "10_full_model_averaged_coefficients.csv")
)

readr::write_csv(
  term_weight_table,
  file.path(output_dir, "11_term_summed_model_weights.csv")
)


# ------------------------------------------------------------
# 10. GVIF FOR GLOBAL AND EVERY SUPPORTED MODEL
# ------------------------------------------------------------

# Refit fresh diagnostic copies rather than relying on a combined object from
# the live R session. Building the list element by element prevents fitted
# model objects from being simplified into atomic vectors.
diagnostic_global_gaussian <- glm(
  formula = global_formula,
  data = analysis_data,
  family = gaussian(link = "identity"),
  na.action = na.fail
)

diagnostic_supported_models <- MuMIn::get.models(
  candidate_table,
  subset = delta < DELTA_AICC_CUTOFF
)

if (length(diagnostic_supported_models) != nrow(supported_table_df)) {
  stop(
    "The freshly refitted diagnostic model count does not match the ",
    "Delta AICc < 2 selection table."
  )
}

diagnostic_supported_labels <- sprintf(
  "supported_%03d_model_%s",
  seq_along(diagnostic_supported_models),
  supported_table_df$model_id
)

models_for_diagnostics <- list()
models_for_diagnostics[["global_gaussian"]] <- diagnostic_global_gaussian

for (i in seq_along(diagnostic_supported_models)) {
  models_for_diagnostics[[diagnostic_supported_labels[i]]] <-
    diagnostic_supported_models[[i]]
}

gvif_model_check <- vapply(
  models_for_diagnostics,
  inherits,
  logical(1),
  what = "glm"
)

if (!all(gvif_model_check)) {
  stop(
    "These GVIF objects are not fitted glm models: ",
    paste(names(gvif_model_check)[!gvif_model_check], collapse = ", ")
  )
}

gvif_results <- dplyr::bind_rows(
  lapply(
    seq_along(models_for_diagnostics),
    function(i) {
      extract_gvif(
        model = models_for_diagnostics[[i]],
        model_label = names(models_for_diagnostics)[i]
      )
    }
  )
)

readr::write_csv(
  gvif_results,
  file.path(output_dir, "12_GVIF_global_and_supported_models.csv")
)


# ------------------------------------------------------------
# 11. MODEL-LEVEL RESIDUAL AND ASSUMPTION TESTS
# ------------------------------------------------------------

safe_test_p_value <- function(expression) {
  tryCatch(
    as.numeric(expression$p.value),
    error = function(e) NA_real_
  )
}

diagnose_model <- function(model, model_label) {
  if (!inherits(model, "glm")) {
    stop(
      "Element '", model_label, "' is not a fitted glm object; class is: ",
      paste(class(model), collapse = ", ")
    )
  }

  diagnostic_lm <- diagnostic_lm_from_glm(model)
  n <- stats::nobs(diagnostic_lm)
  number_of_coefficients <- length(stats::coef(diagnostic_lm))
  residuals_response <- stats::residuals(diagnostic_lm)

  studentised <- tryCatch(
    stats::rstudent(diagnostic_lm),
    error = function(e) rep(NA_real_, n)
  )

  cook <- tryCatch(
    stats::cooks.distance(diagnostic_lm),
    error = function(e) rep(NA_real_, n)
  )

  leverage <- tryCatch(
    stats::hatvalues(diagnostic_lm),
    error = function(e) rep(NA_real_, n)
  )

  model_dffits <- tryCatch(
    stats::dffits(diagnostic_lm),
    error = function(e) rep(NA_real_, n)
  )

  model_dfbetas <- tryCatch(
    stats::dfbetas(diagnostic_lm),
    error = function(e) NULL
  )

  max_abs_dfbeta <- if (is.null(model_dfbetas)) {
    rep(NA_real_, n)
  } else {
    apply(abs(model_dfbetas), 1, max, na.rm = TRUE)
  }

  cook_threshold <- COOKS_MULTIPLIER / n
  leverage_threshold <- LEVERAGE_MULTIPLIER * number_of_coefficients / n
  dffits_threshold <- 2 * sqrt(number_of_coefficients / n)
  dfbeta_threshold <- 2 / sqrt(n)

  shapiro_p <- if (n >= 3 && n <= 5000) {
    finite_studentised <- studentised[is.finite(studentised)]
    if (length(finite_studentised) >= 3) {
      safe_test_p_value(shapiro.test(finite_studentised))
    } else {
      NA_real_
    }
  } else {
    NA_real_
  }

  breusch_pagan_p <- safe_test_p_value(
    lmtest::bptest(diagnostic_lm)
  )

  reset_p <- tryCatch(
    as.numeric(
      lmtest::resettest(
        diagnostic_lm,
        power = 2:3,
        type = "fitted"
      )$p.value
    ),
    error = function(e) NA_real_
  )

  tibble(
    model = model_label,
    diagnostic_error = NA_character_,
    n = n,
    number_of_coefficients = number_of_coefficients,
    residual_df = diagnostic_lm$df.residual,
    converged = model_convergence_status(model),
    rank_deficient = diagnostic_lm$rank < ncol(stats::model.matrix(diagnostic_lm)),
    logLik = as.numeric(stats::logLik(model)),
    AICc = tryCatch(MuMIn::AICc(model), error = function(e) NA_real_),
    residual_mean = mean(residuals_response),
    residual_sd = sd(residuals_response),
    RMSE = sqrt(mean(residuals_response^2)),
    shapiro_wilk_p = shapiro_p,
    breusch_pagan_p = breusch_pagan_p,
    ramsey_RESET_p = reset_p,
    minimum_fitted_value = min(stats::fitted(model)),
    number_negative_fitted = sum(stats::fitted(model) < 0),
    maximum_absolute_studentised_residual = max(abs(studentised), na.rm = TRUE),
    maximum_Cooks_distance = max(cook, na.rm = TRUE),
    Cooks_threshold_4_over_n = cook_threshold,
    number_Cooks_flags = sum(cook > cook_threshold, na.rm = TRUE),
    maximum_leverage = max(leverage, na.rm = TRUE),
    leverage_threshold_2p_over_n = leverage_threshold,
    number_leverage_flags = sum(leverage > leverage_threshold, na.rm = TRUE),
    number_studentised_residual_flags = sum(
      abs(studentised) > STUDENTISED_RESIDUAL_THRESHOLD,
      na.rm = TRUE
    ),
    DFFITS_threshold = dffits_threshold,
    number_DFFITS_flags = sum(
      abs(model_dffits) > dffits_threshold,
      na.rm = TRUE
    ),
    DFBETA_threshold = dfbeta_threshold,
    number_DFBETA_flags = sum(
      max_abs_dfbeta > dfbeta_threshold,
      na.rm = TRUE
    )
  )
}

safe_diagnose_model <- function(model, model_label) {
  tryCatch(
    diagnose_model(model = model, model_label = model_label),
    error = function(e) {
      tibble(
        model = model_label,
        diagnostic_error = conditionMessage(e)
      )
    }
  )
}

# Validate the collection before diagnostics. This gives a clear error if a
# model object was accidentally simplified or overwritten in the R session.
model_object_check <- vapply(
  models_for_diagnostics,
  inherits,
  logical(1),
  what = "glm"
)

if (!all(model_object_check)) {
  stop(
    "These diagnostic objects are not fitted glm models: ",
    paste(names(model_object_check)[!model_object_check], collapse = ", ")
  )
}

# Use base iteration rather than imap_dfr/map2 so complete fitted-model objects
# are passed without any possibility of simplification.
model_diagnostic_list <- lapply(
  seq_along(models_for_diagnostics),
  function(i) {
    safe_diagnose_model(
      model = models_for_diagnostics[[i]],
      model_label = names(models_for_diagnostics)[i]
    )
  }
)

model_diagnostic_tests <- dplyr::bind_rows(model_diagnostic_list)

readr::write_csv(
  model_diagnostic_tests,
  file.path(output_dir, "13_model_residual_and_assumption_tests.csv")
)


# ------------------------------------------------------------
# 12. OBSERVATION-LEVEL INFLUENCE DIAGNOSTICS
# ------------------------------------------------------------

extract_influence <- function(model, model_label, source_data) {
  if (!inherits(model, "glm")) {
    stop(
      "Influence input '", model_label,
      "' is not a fitted glm object; class is: ",
      paste(class(model), collapse = ", ")
    )
  }

  diagnostic_lm <- diagnostic_lm_from_glm(model)
  n <- stats::nobs(diagnostic_lm)
  number_of_coefficients <- length(stats::coef(diagnostic_lm))

  if (nrow(source_data) != n) {
    stop(
      "The source data and model have different observation counts for ",
      model_label, "."
    )
  }

  cook <- tryCatch(
    stats::cooks.distance(diagnostic_lm),
    error = function(e) rep(NA_real_, n)
  )

  leverage <- tryCatch(
    stats::hatvalues(diagnostic_lm),
    error = function(e) rep(NA_real_, n)
  )

  studentised <- tryCatch(
    stats::rstudent(diagnostic_lm),
    error = function(e) rep(NA_real_, n)
  )

  model_dffits <- tryCatch(
    stats::dffits(diagnostic_lm),
    error = function(e) rep(NA_real_, n)
  )

  model_dfbetas <- tryCatch(
    stats::dfbetas(diagnostic_lm),
    error = function(e) NULL
  )

  max_abs_dfbeta <- if (is.null(model_dfbetas)) {
    rep(NA_real_, n)
  } else {
    apply(abs(model_dfbetas), 1, max, na.rm = TRUE)
  }

  cook_threshold <- COOKS_MULTIPLIER / n
  leverage_threshold <- LEVERAGE_MULTIPLIER * number_of_coefficients / n
  dffits_threshold <- 2 * sqrt(number_of_coefficients / n)
  dfbeta_threshold <- 2 / sqrt(n)

  tibble(
    model = model_label,
    diagnostic_error = NA_character_,
    observation_id = source_data$.observation_id,
    fitted = stats::fitted(diagnostic_lm),
    response_residual = stats::residuals(diagnostic_lm),
    studentised_residual = studentised,
    leverage = leverage,
    Cooks_distance = cook,
    DFFITS = model_dffits,
    maximum_absolute_DFBETA = max_abs_dfbeta,
    Cooks_threshold = cook_threshold,
    leverage_threshold = leverage_threshold,
    DFFITS_threshold = dffits_threshold,
    DFBETA_threshold = dfbeta_threshold,
    flag_Cooks = Cooks_distance > Cooks_threshold,
    flag_leverage = leverage > leverage_threshold,
    flag_studentised_residual =
      abs(studentised_residual) > STUDENTISED_RESIDUAL_THRESHOLD,
    flag_DFFITS = abs(DFFITS) > DFFITS_threshold,
    flag_DFBETA = maximum_absolute_DFBETA > DFBETA_threshold,
    flag_any_influence_measure =
      flag_Cooks |
      flag_leverage |
      flag_studentised_residual |
      flag_DFFITS |
      flag_DFBETA
  )
}

safe_extract_influence <- function(model, model_label, source_data) {
  tryCatch(
    extract_influence(
      model = model,
      model_label = model_label,
      source_data = source_data
    ),
    error = function(e) {
      tibble(
        model = model_label,
        diagnostic_error = conditionMessage(e),
        observation_id = NA_integer_,
        Cooks_distance = NA_real_,
        flag_Cooks = FALSE,
        flag_any_influence_measure = FALSE
      )
    }
  )
}

influence_result_list <- lapply(
  seq_along(models_for_diagnostics),
  function(i) {
    safe_extract_influence(
      model = models_for_diagnostics[[i]],
      model_label = names(models_for_diagnostics)[i],
      source_data = analysis_data
    )
  }
)

influence_results <- dplyr::bind_rows(influence_result_list)

readr::write_csv(
  influence_results,
  file.path(output_dir, "14_observation_influence_diagnostics.csv")
)

readr::write_csv(
  influence_results %>%
    filter(flag_any_influence_measure) %>%
    arrange(model, desc(Cooks_distance)),
  file.path(output_dir, "15_flagged_observations_all_measures.csv")
)


# ------------------------------------------------------------
# 13. STANDARD RESIDUAL DIAGNOSTIC PLOTS
# ------------------------------------------------------------

pdf(
  file.path(output_dir, "16_residual_diagnostic_plots_all_models.pdf"),
  width = 11,
  height = 8.5,
  onefile = TRUE
)

for (model_label in names(models_for_diagnostics)) {
  model <- models_for_diagnostics[[model_label]]
  par(mfrow = c(2, 2), oma = c(0, 0, 2.5, 0))
  plot(model, which = 1, main = "Residuals vs fitted")
  plot(model, which = 2, main = "Normal Q-Q")
  plot(model, which = 3, main = "Scale-location")
  plot(model, which = 5, main = "Residuals vs leverage")
  mtext(model_label, outer = TRUE, cex = 1.1, font = 2)
}

dev.off()

pdf(
  file.path(output_dir, "17_Cooks_distance_supported_models.pdf"),
  width = 11,
  height = 7,
  onefile = TRUE
)

for (model_label in names(supported_models)) {
  model <- supported_models[[model_label]]
  cook <- cooks.distance(model)
  threshold <- COOKS_MULTIPLIER / nobs(model)

  plot(
    cook,
    type = "h",
    lwd = 1.2,
    main = paste("Cook's distance:", model_label),
    xlab = "Observation position in analysis dataset",
    ylab = "Cook's distance"
  )
  abline(h = threshold, lty = 2, lwd = 2, col = "red")
  flagged <- which(cook > threshold)
  if (length(flagged) > 0) {
    text(
      flagged,
      cook[flagged],
      labels = analysis_data$.observation_id[flagged],
      pos = 3,
      cex = 0.7
    )
  }
}

dev.off()

# Component-plus-residual plots examine functional form for continuous terms
# in the global Gaussian model. Failure is recorded rather than stopping all
# other diagnostics.
cr_plot_status <- tryCatch(
  {
    pdf(
      file.path(output_dir, "18_component_residual_plots_global_model.pdf"),
      width = 11,
      height = 8.5
    )
    car::crPlots(global_gaussian)
    dev.off()
    "Created"
  },
  error = function(e) {
    try(dev.off(), silent = TRUE)
    paste("Not created:", e$message)
  }
)


# ------------------------------------------------------------
# 14. MODEL-AVERAGED FITTED VALUES AND RESIDUAL CHECKS
# ------------------------------------------------------------

averaged_fitted <- as.numeric(
  predict(
    averaged_model,
    newdata = analysis_data,
    type = "response",
    se.fit = FALSE
  )
)

averaged_residual <- analysis_data$mean_distance - averaged_fitted
averaged_standardised_residual <- averaged_residual / sd(averaged_residual)

averaged_residual_data <- tibble(
  observation_id = analysis_data$.observation_id,
  observed = analysis_data$mean_distance,
  model_averaged_fitted = averaged_fitted,
  model_averaged_residual = averaged_residual,
  model_averaged_standardised_residual = averaged_standardised_residual,
  negative_fitted_value = averaged_fitted < 0
)

averaged_residual_lm <- lm(
  averaged_residual ~ averaged_fitted
)

model_averaged_residual_tests <- tibble(
  n = nrow(analysis_data),
  residual_mean = mean(averaged_residual),
  residual_sd = sd(averaged_residual),
  RMSE = sqrt(mean(averaged_residual^2)),
  shapiro_wilk_p = shapiro.test(averaged_standardised_residual)$p.value,
  breusch_pagan_p = lmtest::bptest(averaged_residual_lm)$p.value,
  minimum_fitted_value = min(averaged_fitted),
  number_negative_fitted = sum(averaged_fitted < 0),
  maximum_absolute_standardised_residual =
    max(abs(averaged_standardised_residual))
)

readr::write_csv(
  averaged_residual_data,
  file.path(output_dir, "19_model_averaged_fitted_and_residuals.csv")
)

readr::write_csv(
  model_averaged_residual_tests,
  file.path(output_dir, "20_model_averaged_residual_tests.csv")
)

png(
  file.path(output_dir, "21_model_averaged_residual_plots.png"),
  width = 1800,
  height = 1600,
  res = 180
)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 1))
plot(
  averaged_fitted,
  averaged_residual,
  xlab = "Model-averaged fitted fishing distance (km)",
  ylab = "Residual (km)",
  main = "Residuals vs fitted"
)
abline(h = 0, lty = 2, col = "red")
lines(lowess(averaged_fitted, averaged_residual), lwd = 2, col = "blue")
qqnorm(
  averaged_standardised_residual,
  main = "Normal Q-Q: model-averaged residuals"
)
qqline(averaged_standardised_residual, col = "red", lwd = 2)
plot(
  averaged_fitted,
  sqrt(abs(averaged_standardised_residual)),
  xlab = "Model-averaged fitted fishing distance (km)",
  ylab = "Square root |standardised residual|",
  main = "Scale-location"
)
lines(
  lowess(
    averaged_fitted,
    sqrt(abs(averaged_standardised_residual))
  ),
  lwd = 2,
  col = "blue"
)
hist(
  averaged_residual,
  breaks = "FD",
  main = "Model-averaged residual distribution",
  xlab = "Residual (km)",
  col = "grey85",
  border = "grey30"
)
dev.off()


# ------------------------------------------------------------
# 15. COOK'S-DISTANCE UNION ACROSS SUPPORTED MODELS
# ------------------------------------------------------------

supported_influence <- influence_results %>%
  filter(model != "global_gaussian")

influential_ids <- supported_influence %>%
  filter(flag_Cooks) %>%
  distinct(observation_id) %>%
  arrange(observation_id)

influential_id_summary <- influential_ids %>%
  left_join(
    supported_influence %>%
      group_by(observation_id) %>%
      summarise(
        supported_models_flagged_by_Cooks = sum(flag_Cooks),
        maximum_Cooks_distance = max(Cooks_distance),
        .groups = "drop"
      ),
    by = "observation_id"
  )

readr::write_csv(
  influential_id_summary,
  file.path(output_dir, "22_Cooks_union_supported_models.csv")
)


# ------------------------------------------------------------
# 16. COOK'S-DISTANCE SENSITIVITY ANALYSIS
# ------------------------------------------------------------

sensitivity_status <- "Not run: no Cook's-distance flags in supported models"
sensitivity_coefficients <- tibble()
coefficient_comparison <- tibble()
sensitivity_summary <- tibble(
  status = sensitivity_status,
  original_n = nrow(analysis_data),
  influential_observations_removed = 0L,
  sensitivity_n = nrow(analysis_data),
  original_candidate_models = nrow(candidate_table_df),
  original_supported_models_Delta_lt2 = nrow(supported_table_df),
  original_supported_set_weight_before_renormalisation =
    sum(supported_table_df$weight),
  sensitivity_scope =
    "Same original supported model set; no observations required exclusion",
  sensitivity_refitted_models = NA_integer_,
  sensitivity_refitted_model_weight_sum = NA_real_
)

if (nrow(influential_ids) > 0) {
  sensitivity_data <- analysis_data %>%
    filter(!.observation_id %in% influential_ids$observation_id) %>%
    droplevels()

  sensitivity_factor_levels <- vapply(
    sensitivity_data %>% select(where(is.factor)),
    nlevels,
    integer(1)
  )

  if (any(sensitivity_factor_levels < 2)) {
    stop(
      "Cook's-distance exclusion removed all observations from at least one ",
      "factor level. A comparable sensitivity model cannot be fitted."
    )
  }

  sensitivity_global <- glm(
    formula = global_formula,
    data = sensitivity_data,
    family = gaussian(link = "identity"),
    na.action = na.fail
  )

  if (identical(model_convergence_status(sensitivity_global), FALSE)) {
    stop("The Cook's-distance sensitivity global model did not converge.")
  }

  if (qr(model.matrix(sensitivity_global))$rank <
      ncol(model.matrix(sensitivity_global))) {
    stop("The Cook's-distance sensitivity global model is rank deficient.")
  }

  # Refit the same 10 originally supported models. This directly isolates
  # sensitivity to influential observations and avoids another 65,536-model
  # dredge operation that can exhaust R's C stack.
  sensitivity_supported_models <- vector(
    "list",
    length(supported_models)
  )

  for (i in seq_along(supported_models)) {
    sensitivity_supported_models[[i]] <- glm(
      formula = stats::formula(supported_models[[i]]),
      data = sensitivity_data,
      family = gaussian(link = "identity"),
      na.action = na.fail
    )
  }

  sensitivity_model_ids <- supported_model_ids
  names(sensitivity_supported_models) <- sprintf(
    "sensitivity_supported_%03d_model_%s",
    seq_along(sensitivity_supported_models),
    sensitivity_model_ids
  )

  sensitivity_AICc <- vapply(
    sensitivity_supported_models,
    MuMIn::AICc,
    numeric(1)
  )
  sensitivity_delta <- sensitivity_AICc - min(sensitivity_AICc)
  sensitivity_weight <- exp(-0.5 * sensitivity_delta)
  sensitivity_weight <- sensitivity_weight / sum(sensitivity_weight)

  sensitivity_candidate_df <- tibble(
    original_model_id = sensitivity_model_ids,
    refitted_model = names(sensitivity_supported_models),
    df = vapply(
      sensitivity_supported_models,
      function(x) attr(stats::logLik(x), "df"),
      numeric(1)
    ),
    logLik = vapply(
      sensitivity_supported_models,
      function(x) as.numeric(stats::logLik(x)),
      numeric(1)
    ),
    AICc = sensitivity_AICc,
    delta = sensitivity_delta,
    weight = sensitivity_weight,
    sensitivity_scope =
      "Original Delta-AICc-below-2 model set refitted after Cook exclusion"
  ) %>%
    arrange(AICc)

  sensitivity_average <- MuMIn::model.avg(
    sensitivity_supported_models,
    revised.var = TRUE
  )

  sensitivity_weight_table <- as.data.frame(
    sensitivity_average$msTable
  ) %>%
    rownames_to_column("component_model")

  sensitivity_coefficients <- extract_full_average(sensitivity_average)

  coefficient_comparison <- full_join(
    original_averaged_coefficients %>%
      rename_with(~ paste0(.x, "_original"), -term),
    sensitivity_coefficients %>%
      rename_with(~ paste0(.x, "_sensitivity"), -term),
    by = "term"
  ) %>%
    mutate(
      estimate_difference = estimate_sensitivity - estimate_original,
      absolute_estimate_difference = abs(estimate_difference),
      proportional_change = if_else(
        !is.na(estimate_original) & abs(estimate_original) > 1e-12,
        estimate_difference / abs(estimate_original),
        NA_real_
      ),
      sign_changed = case_when(
        is.na(estimate_original) | is.na(estimate_sensitivity) ~ NA,
        TRUE ~ sign(estimate_original) != sign(estimate_sensitivity)
      ),
      original_CI_crosses_zero = case_when(
        is.na(conf.low_original) | is.na(conf.high_original) ~ NA,
        TRUE ~ conf.low_original <= 0 & conf.high_original >= 0
      ),
      sensitivity_CI_crosses_zero = case_when(
        is.na(conf.low_sensitivity) | is.na(conf.high_sensitivity) ~ NA,
        TRUE ~ conf.low_sensitivity <= 0 & conf.high_sensitivity >= 0
      ),
      CI_zero_conclusion_changed = case_when(
        is.na(original_CI_crosses_zero) |
          is.na(sensitivity_CI_crosses_zero) ~ NA,
        TRUE ~ original_CI_crosses_zero != sensitivity_CI_crosses_zero
      )
    )

  sensitivity_models_for_diagnostics <- list()
  sensitivity_models_for_diagnostics[["sensitivity_global"]] <-
    sensitivity_global

  for (i in seq_along(sensitivity_supported_models)) {
    sensitivity_models_for_diagnostics[[
      names(sensitivity_supported_models)[i]
    ]] <- sensitivity_supported_models[[i]]
  }

  sensitivity_object_check <- vapply(
    sensitivity_models_for_diagnostics,
    inherits,
    logical(1),
    what = "glm"
  )

  if (!all(sensitivity_object_check)) {
    stop(
      "These sensitivity objects are not fitted glm models: ",
      paste(
        names(sensitivity_object_check)[!sensitivity_object_check],
        collapse = ", "
      )
    )
  }

  sensitivity_diagnostic_tests <- dplyr::bind_rows(
    lapply(
      seq_along(sensitivity_models_for_diagnostics),
      function(i) {
        safe_diagnose_model(
          model = sensitivity_models_for_diagnostics[[i]],
          model_label = names(sensitivity_models_for_diagnostics)[i]
        )
      }
    )
  )

  sensitivity_influence <- dplyr::bind_rows(
    lapply(
      seq_along(sensitivity_models_for_diagnostics),
      function(i) {
        safe_extract_influence(
          model = sensitivity_models_for_diagnostics[[i]],
          model_label = names(sensitivity_models_for_diagnostics)[i],
          source_data = sensitivity_data
        )
      }
    )
  )

  sensitivity_gvif <- dplyr::bind_rows(
    lapply(
      seq_along(sensitivity_models_for_diagnostics),
      function(i) {
        extract_gvif(
          model = sensitivity_models_for_diagnostics[[i]],
          model_label = names(sensitivity_models_for_diagnostics)[i]
        )
      }
    )
  )

  sensitivity_summary <- tibble(
    status = "Completed",
    original_n = nrow(analysis_data),
    influential_observations_removed = nrow(influential_ids),
    sensitivity_n = nrow(sensitivity_data),
    original_candidate_models = nrow(candidate_table_df),
    original_supported_models_Delta_lt2 = nrow(supported_table_df),
    original_supported_set_weight_before_renormalisation =
      sum(supported_table_df$weight),
    sensitivity_scope = paste(
      "The same", length(supported_models),
      "original Delta-AICc-below-2 models were refitted after Cook exclusion"
    ),
    sensitivity_refitted_models = length(sensitivity_supported_models),
    sensitivity_refitted_model_weight_sum = sum(sensitivity_weight)
  )

  readr::write_csv(
    sensitivity_candidate_df,
    file.path(
      output_dir,
      "23_sensitivity_refitted_supported_models_AICc.csv"
    )
  )

  readr::write_csv(
    sensitivity_weight_table,
    file.path(output_dir, "24_sensitivity_model_averaging_weights.csv")
  )

  readr::write_csv(
    sensitivity_coefficients,
    file.path(output_dir, "25_sensitivity_full_averaged_coefficients.csv")
  )

  readr::write_csv(
    coefficient_comparison,
    file.path(output_dir, "26_original_vs_sensitivity_coefficients.csv")
  )

  readr::write_csv(
    sensitivity_diagnostic_tests,
    file.path(output_dir, "27_sensitivity_residual_and_assumption_tests.csv")
  )

  readr::write_csv(
    sensitivity_influence,
    file.path(output_dir, "28_sensitivity_observation_influence.csv")
  )

  readr::write_csv(
    sensitivity_gvif,
    file.path(output_dir, "29_sensitivity_GVIF.csv")
  )

  pdf(
    file.path(output_dir, "31_sensitivity_residual_diagnostic_plots.pdf"),
    width = 11,
    height = 8.5,
    onefile = TRUE
  )

  for (model_label in names(sensitivity_models_for_diagnostics)) {
    model <- diagnostic_lm_from_glm(
      sensitivity_models_for_diagnostics[[model_label]]
    )
    par(mfrow = c(2, 2), oma = c(0, 0, 2.5, 0))
    plot(model, which = 1, main = "Residuals vs fitted")
    plot(model, which = 2, main = "Normal Q-Q")
    plot(model, which = 3, main = "Scale-location")
    plot(model, which = 5, main = "Residuals vs leverage")
    mtext(model_label, outer = TRUE, cex = 1.1, font = 2)
  }

  dev.off()

  sensitivity_status <- "Completed"
}

# This status table is always written. If no observation exceeds Cook's 4/n
# threshold, the sensitivity reanalysis is correctly recorded as unnecessary.
readr::write_csv(
  sensitivity_summary,
  file.path(output_dir, "30_sensitivity_analysis_summary.csv")
)


# ------------------------------------------------------------
# 17. MASTER ANALYSIS SUMMARY AND SESSION INFORMATION
# ------------------------------------------------------------

master_summary <- tibble(
  imported_n = nrow(fishers_raw),
  complete_case_n = nrow(analysis_data),
  number_of_global_terms = length(attr(terms(global_gaussian), "term.labels")),
  number_of_candidate_models = nrow(candidate_table_df),
  delta_AICc_rule = "strictly less than 2",
  number_of_supported_models = nrow(supported_table_df),
  supported_set_weight_before_renormalisation =
    sum(supported_table_df$weight),
  lowest_AICc = min(candidate_table_df$AICc),
  highest_retained_delta_AICc = max(supported_table_df$delta),
  first_excluded_delta_AICc = min(
    candidate_table_df$delta[
      candidate_table_df$delta >= DELTA_AICC_CUTOFF
    ]
  ),
  observations_flagged_by_Cooks_union = nrow(influential_ids),
  Cooks_threshold = COOKS_MULTIPLIER / nrow(analysis_data),
  sensitivity_status = sensitivity_status,
  component_residual_plot_status = cr_plot_status
)

readr::write_csv(
  master_summary,
  file.path(output_dir, "32_master_analysis_summary.csv")
)

capture.output(
  summary(global_gaussian),
  file = file.path(output_dir, "33_global_Gaussian_model_summary.txt")
)

capture.output(
  summary(averaged_model),
  file = file.path(output_dir, "34_full_model_averaged_summary.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "35_R_session_information.txt")
)

diagnostic_notes <- c(
  "FISHING-DISTANCE DIAGNOSTIC OUTPUT NOTES",
  "",
  "Supported-model rule: Delta AICc is strictly less than 2.",
  paste0("Cook's-distance threshold: 4/n = ",
         signif(COOKS_MULTIPLIER / nrow(analysis_data), 6), "."),
  paste0("Studentised-residual flag: absolute value > ",
         STUDENTISED_RESIDUAL_THRESHOLD, "."),
  "Leverage flag: leverage > 2p/n, where p is the number of fitted coefficients.",
  "DFFITS flag: absolute DFFITS > 2*sqrt(p/n).",
  "DFBETA flag: maximum absolute DFBETA > 2/sqrt(n).",
  "Adjusted GVIF is GVIF^(1/(2*df)); outputs flag values >=1.5 and >=2.",
  "Shapiro-Wilk and Breusch-Pagan p-values should be interpreted together with the residual plots, not alone.",
  "The Ramsey RESET test is a general functional-form check; inspect component-residual plots for continuous predictors.",
  "The Gamma comparison is a distributional sensitivity check only. Gaussian and Gamma coefficients are not model-averaged together.",
  "Cook's sensitivity removes the union of observations exceeding 4/n in any supported original model, refits the same original supported model set, recalculates AICc weights, and repeats full model averaging.",
  "Do not automatically delete influential observations from the primary analysis. Use the sensitivity comparison to assess whether conclusions change."
)

writeLines(
  diagnostic_notes,
  con = file.path(output_dir, "36_READ_ME_diagnostic_thresholds.txt")
)

# Verify that a successful run has produced every non-conditional diagnostic
# output. Files 23--29 and 31 are conditional on at least one Cook's-distance
# flag; file 30 records whether that conditional analysis was required.
mandatory_diagnostic_files <- c(
  "13_model_residual_and_assumption_tests.csv",
  "14_observation_influence_diagnostics.csv",
  "15_flagged_observations_all_measures.csv",
  "16_residual_diagnostic_plots_all_models.pdf",
  "17_Cooks_distance_supported_models.pdf",
  "18_component_residual_plots_global_model.pdf",
  "19_model_averaged_fitted_and_residuals.csv",
  "20_model_averaged_residual_tests.csv",
  "21_model_averaged_residual_plots.png",
  "22_Cooks_union_supported_models.csv",
  "30_sensitivity_analysis_summary.csv",
  "32_master_analysis_summary.csv",
  "33_global_Gaussian_model_summary.txt",
  "34_full_model_averaged_summary.txt",
  "35_R_session_information.txt",
  "36_READ_ME_diagnostic_thresholds.txt"
)

missing_diagnostic_files <- mandatory_diagnostic_files[
  !file.exists(file.path(output_dir, mandatory_diagnostic_files))
]

if (length(missing_diagnostic_files) > 0) {
  stop(
    "The run ended without these mandatory diagnostic files: ",
    paste(missing_diagnostic_files, collapse = ", ")
  )
}

message(
  "Complete diagnostic analysis finished. Results saved in: ",
  normalizePath(output_dir)
)


# =============================================================================
# ============================================================
# PART B. HC3 + PAIRS-BOOTSTRAP ROBUSTNESS ANALYSIS
# ============================================================
#
# PURPOSE
# -------
# This section performs the final robustness analysis for the fishing-distance
# models generated in Part A. It:
#
#   1. uses all 148 complete observations as the primary dataset;
#   2. reconstructs and verifies the same 10 models selected in Part A by
#      the strict Delta AICc < 2 rule;
#   3. identifies the union of observations with Cook's distance > 4/n and
#      exports their original records for source-data verification;
#   4. DOES NOT automatically remove any influential observation;
#   5. calculates full model-averaged estimates with HC3
#      heteroscedasticity-consistent within-model variances;
#   6. performs 2,000 non-parametric pairs-bootstrap repetitions, refitting
#      and reweighting the same 10 supported models in every repetition;
#   7. produces percentile bootstrap confidence intervals and compares them
#      with the original and HC3 intervals.
#
# The robustness step deliberately does NOT repeat the 65,536-model search.
# Bootstrap inference is conditional on the verified 10-model set from Part A.
#
# INPUTS FOR PART B
# -----------------
#   model_data_recoded_full.csv
#   Fishing_distance_diagnostics_DeltaAICc_lt2/08_supported_models_DeltaAICc_lt2.csv
#   Fishing_distance_diagnostics_DeltaAICc_lt2/10_full_model_averaged_coefficients.csv
#
# Part A creates the two analysis-output CSV files above automatically.
#
# INSTALLATION
# ------------
# Run once if sandwich is not installed:
#   install.packages("sandwich")
#
# RUNNING
# -------
# Run the complete file from the project directory containing
# model_data_recoded_full.csv. Part A must finish before this section runs.
# Do not run isolated sections.
# ============================================================


# ------------------------------------------------------------
# 0. FIXED SETTINGS
# ------------------------------------------------------------

required_packages <- c("sandwich")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following package(s), restart R, and rerun the script: ",
    paste(missing_packages, collapse = ", "),
    ". Installation command: install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "),
    "))"
  )
}

options(na.action = "na.fail")

EXPECTED_COMPLETE_N <- 148L
EXPECTED_SUPPORTED_MODELS <- 10L
DELTA_AICC_CUTOFF <- 2
COOKS_MULTIPLIER <- 4
AICC_REPRODUCTION_TOLERANCE <- 0.05

# Increase this only if more precision is required. Two thousand successful
# repetitions are suitable for the requested final robustness check.
N_BOOTSTRAP <- 2000L
RANDOM_SEED <- 20260919L
MAX_BOOTSTRAP_ATTEMPTS <- N_BOOTSTRAP * 10L
PROGRESS_EVERY <- 100L

# These are the previously verified analysis-row identifiers. Keeping this
# check TRUE protects against accidentally selecting a different CSV. If the
# source data are intentionally corrected, set it to FALSE and rerun the full
# original analysis before using this robustness script again.
VERIFY_EXPECTED_COOK_UNION <- TRUE
EXPECTED_COOK_UNION <- c(22L, 95L, 102L, 103L, 115L, 121L, 144L)

OUTPUT_DIR <- "Fishing_distance_HC3_bootstrap_robustness"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# 1. LOCATE INPUT FILES SAFELY
# ------------------------------------------------------------

locate_file <- function(filename, required = TRUE) {
  search_paths <- unique(c(
    filename,
    file.path("upload", filename),
    file.path("Fishing_distance_diagnostics_DeltaAICc_lt2", filename)
  ))
  
  found <- search_paths[file.exists(search_paths)]
  
  if (length(found) == 0L) {
    if (required) {
      stop(
        "Required file not found: ", filename, ". Place it in the working ",
        "directory or in an 'upload' subfolder."
      )
    }
    return(NA_character_)
  }
  
  normalizePath(found[1L], winslash = "/", mustWork = TRUE)
}

data_file <- locate_file("model_data_recoded_full.csv")
supported_file <- locate_file("08_supported_models_DeltaAICc_lt2.csv")
original_coefficient_file <- locate_file(
  "10_full_model_averaged_coefficients.csv",
  required = FALSE
)

message("DATA FILE USED: ", data_file)
message("SUPPORTED-MODEL FILE USED: ", supported_file)
if (!is.na(original_coefficient_file)) {
  message("ORIGINAL COEFFICIENT FILE USED: ", original_coefficient_file)
}


# ------------------------------------------------------------
# 2. RECREATE THE VERIFIED 148-OBSERVATION ANALYSIS DATASET
# ------------------------------------------------------------

fishers_raw <- read.csv(
  data_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(fishers_raw) != EXPECTED_COMPLETE_N) {
  stop(
    "The selected raw CSV contains ", nrow(fishers_raw),
    " rows, but the verified dataset contains ", EXPECTED_COMPLETE_N,
    ". Check the file printed above."
  )
}

fishers_raw$.observation_id <- seq_len(nrow(fishers_raw))

required_raw_columns <- c(
  "mean_distance",
  "Site",
  "Fishing.experience",
  "vessel_length_m",
  "cost_per_trip_usd",
  "AIS",
  "VMS",
  "navigation.apps",
  "target.sharks",
  "subsidy",
  "Top.species.group",
  "communication",
  "reg_awareness",
  "bycatch_raw",
  "gear_type_clean"
)

missing_raw_columns <- setdiff(required_raw_columns, names(fishers_raw))
if (length(missing_raw_columns) > 0L) {
  stop(
    "The selected raw CSV is missing: ",
    paste(missing_raw_columns, collapse = ", ")
  )
}

fishers_raw$bycatch_clean <- ifelse(
  fishers_raw$bycatch_raw == "Never",
  "Never",
  ifelse(
    fishers_raw$bycatch_raw %in% c(
      "Rarely", "Occasionally", "Sometimes"
    ),
    "Infrequent",
    ifelse(
      fishers_raw$bycatch_raw == "On every trip",
      "Every trip",
      NA_character_
    )
  )
)

fishers_raw$bycatch_clean <- ordered(
  fishers_raw$bycatch_clean,
  levels = c("Never", "Infrequent", "Every trip")
)

make_binary_gear <- function(x, pattern) {
  answer <- ifelse(
    is.na(x),
    NA_character_,
    ifelse(grepl(pattern, x, ignore.case = TRUE), "Yes", "No")
  )
  relevel(factor(answer), ref = "No")
}

fishers_raw$Site <- relevel(factor(fishers_raw$Site), ref = "A")
fishers_raw$AIS <- relevel(factor(fishers_raw$AIS), ref = "No")
fishers_raw$VMS <- relevel(factor(fishers_raw$VMS), ref = "No")
fishers_raw$navigation.apps <- relevel(
  factor(fishers_raw$navigation.apps), ref = "No"
)
fishers_raw$target.sharks <- relevel(
  factor(fishers_raw$target.sharks), ref = "No"
)
fishers_raw$subsidy <- relevel(factor(fishers_raw$subsidy), ref = "No")
fishers_raw$Top.species.group <- relevel(
  factor(fishers_raw$Top.species.group), ref = "Cephalopods"
)
fishers_raw$communication <- relevel(
  factor(fishers_raw$communication), ref = "Every trip"
)
fishers_raw$reg_awareness <- relevel(
  factor(fishers_raw$reg_awareness), ref = "Fully aware"
)
fishers_raw$uses_gillnet <- make_binary_gear(
  fishers_raw$gear_type_clean, "Gillnet"
)
fishers_raw$uses_handline <- make_binary_gear(
  fishers_raw$gear_type_clean, "Handline"
)
fishers_raw$uses_longline <- make_binary_gear(
  fishers_raw$gear_type_clean, "Longline"
)

analysis_variables <- c(
  ".observation_id",
  "mean_distance",
  "Site",
  "Fishing.experience",
  "vessel_length_m",
  "cost_per_trip_usd",
  "AIS",
  "VMS",
  "navigation.apps",
  "target.sharks",
  "subsidy",
  "Top.species.group",
  "communication",
  "reg_awareness",
  "uses_gillnet",
  "uses_handline",
  "uses_longline",
  "bycatch_clean"
)

analysis_data <- fishers_raw[
  complete.cases(fishers_raw[, analysis_variables, drop = FALSE]),
  analysis_variables,
  drop = FALSE
]
analysis_data <- droplevels(analysis_data)

if (nrow(analysis_data) != EXPECTED_COMPLETE_N) {
  missing_counts <- colSums(
    is.na(fishers_raw[, analysis_variables, drop = FALSE])
  )
  incomplete_rows <- which(
    !complete.cases(fishers_raw[, analysis_variables, drop = FALSE])
  )
  
  write.csv(
    data.frame(
      variable = names(missing_counts),
      missing_values = as.integer(missing_counts),
      row.names = NULL
    ),
    file.path(OUTPUT_DIR, "00_missing_value_audit.csv"),
    row.names = FALSE
  )
  
  if (length(incomplete_rows) > 0L) {
    write.csv(
      fishers_raw[incomplete_rows, , drop = FALSE],
      file.path(OUTPUT_DIR, "00_incomplete_source_rows.csv"),
      row.names = FALSE
    )
  }
  
  stop(
    "Expected ", EXPECTED_COMPLETE_N,
    " complete observations, but obtained ", nrow(analysis_data),
    ". See the two file-00 audits in the output folder."
  )
}

if (anyDuplicated(analysis_data$.observation_id)) {
  stop("Observation identifiers are not unique.")
}


# ------------------------------------------------------------
# 3. RECONSTRUCT AND VERIFY THE STRICT DELTA-AICc < 2 MODEL SET
# ------------------------------------------------------------

supported_table <- read.csv(
  supported_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_supported_columns <- c("model_id", "AICc", "delta")
missing_supported_columns <- setdiff(
  required_supported_columns,
  names(supported_table)
)

if (length(missing_supported_columns) > 0L) {
  stop(
    "The supported-model table is missing: ",
    paste(missing_supported_columns, collapse = ", ")
  )
}

if (nrow(supported_table) != EXPECTED_SUPPORTED_MODELS) {
  stop(
    "Expected ", EXPECTED_SUPPORTED_MODELS,
    " supported models, but the selected table contains ",
    nrow(supported_table), "."
  )
}

if (any(!is.finite(supported_table$delta)) ||
    any(supported_table$delta >= DELTA_AICC_CUTOFF)) {
  stop("The selected model table is not a strict Delta AICc < 2 table.")
}

candidate_terms <- c(
  "AIS",
  "bycatch_clean",
  "communication",
  "cost_per_trip_usd",
  "Fishing.experience",
  "navigation.apps",
  "reg_awareness",
  "Site",
  "subsidy",
  "target.sharks",
  "Top.species.group",
  "uses_gillnet",
  "uses_handline",
  "uses_longline",
  "vessel_length_m",
  "VMS"
)

missing_term_columns <- setdiff(candidate_terms, names(supported_table))
if (length(missing_term_columns) > 0L) {
  stop(
    "The supported-model table lacks model-term columns: ",
    paste(missing_term_columns, collapse = ", ")
  )
}

supported_formulas <- vector("list", nrow(supported_table))

for (i in seq_len(nrow(supported_table))) {
  included_terms <- candidate_terms[
    vapply(
      candidate_terms,
      function(term_name) !is.na(supported_table[[term_name]][i]),
      logical(1)
    )
  ]
  
  supported_formulas[[i]] <- reformulate(
    included_terms,
    response = "mean_distance"
  )
}

model_ids <- as.character(supported_table$model_id)
model_labels <- sprintf(
  "supported_%03d_model_%s",
  seq_along(model_ids),
  model_ids
)
names(supported_formulas) <- model_labels

fit_supported_models <- function(data) {
  fitted_models <- lapply(
    supported_formulas,
    function(model_formula) {
      lm(
        model_formula,
        data = data,
        na.action = na.fail,
        model = TRUE,
        x = TRUE,
        y = TRUE
      )
    }
  )
  names(fitted_models) <- model_labels
  fitted_models
}

aicc_lm <- function(model) {
  model_loglik <- logLik(model)
  number_parameters <- attr(model_loglik, "df")
  sample_size <- nobs(model)
  
  if (sample_size <= number_parameters + 1) {
    return(Inf)
  }
  
  -2 * as.numeric(model_loglik) +
    2 * number_parameters +
    (2 * number_parameters * (number_parameters + 1)) /
    (sample_size - number_parameters - 1)
}

aicc_weights <- function(aicc_values) {
  if (any(!is.finite(aicc_values))) {
    stop("At least one AICc value is non-finite.")
  }
  
  delta <- aicc_values - min(aicc_values)
  unscaled <- exp(-0.5 * delta)
  unscaled / sum(unscaled)
}

supported_models <- fit_supported_models(analysis_data)

if (any(vapply(
  supported_models,
  function(model) anyNA(coef(model)) ||
  model$rank < ncol(model.matrix(model)),
  logical(1)
))) {
  stop("At least one original supported model is rank deficient.")
}

refitted_aicc <- vapply(supported_models, aicc_lm, numeric(1))
maximum_aicc_difference <- max(
  abs(refitted_aicc - supported_table$AICc)
)

if (!is.finite(maximum_aicc_difference) ||
    maximum_aicc_difference > AICC_REPRODUCTION_TOLERANCE) {
  stop(
    "The models do not reproduce the saved AICc table. Maximum difference = ",
    signif(maximum_aicc_difference, 8),
    ". Confirm that files 08 and the raw CSV belong to the same run."
  )
}

original_weights <- aicc_weights(refitted_aicc)

refit_verification <- data.frame(
  model = model_labels,
  model_id = model_ids,
  formula = vapply(
    supported_formulas,
    function(x) paste(deparse(x), collapse = " "),
    character(1)
  ),
  saved_AICc = supported_table$AICc,
  refitted_AICc = refitted_aicc,
  absolute_AICc_difference = abs(refitted_aicc - supported_table$AICc),
  refitted_delta = refitted_aicc - min(refitted_aicc),
  refitted_weight = original_weights,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(
  refit_verification,
  file.path(OUTPUT_DIR, "37_robustness_model_refit_verification.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. EXPORT THE COOK-FLAGGED SOURCE RECORDS FOR MANUAL CHECKING
# ------------------------------------------------------------

cook_values <- lapply(supported_models, cooks.distance)
cook_threshold <- COOKS_MULTIPLIER / nrow(analysis_data)

cook_union_ids <- sort(unique(unlist(lapply(
  cook_values,
  function(x) analysis_data$.observation_id[which(x > cook_threshold)]
))))

if (VERIFY_EXPECTED_COOK_UNION &&
    !identical(as.integer(cook_union_ids), EXPECTED_COOK_UNION)) {
  stop(
    "The Cook's-distance union does not match the verified seven rows. ",
    "Obtained: ", paste(cook_union_ids, collapse = ", "),
    ". Expected: ", paste(EXPECTED_COOK_UNION, collapse = ", "),
    ". Confirm the selected inputs."
  )
}

cook_summary <- data.frame(
  analysis_observation_id = cook_union_ids,
  supported_models_flagged = vapply(
    cook_union_ids,
    function(id) {
      position <- match(id, analysis_data$.observation_id)
      sum(vapply(
        cook_values,
        function(x) x[position] > cook_threshold,
        logical(1)
      ))
    },
    integer(1)
  ),
  maximum_Cooks_distance = vapply(
    cook_union_ids,
    function(id) {
      position <- match(id, analysis_data$.observation_id)
      max(vapply(cook_values, function(x) x[position], numeric(1)))
    },
    numeric(1)
  ),
  verification_instruction = paste(
    "Compare with the original questionnaire/database;",
    "do not exclude unless confirmed erroneous"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

source_positions <- match(
  cook_summary$analysis_observation_id,
  fishers_raw$.observation_id
)

flagged_source_records <- cbind(
  cook_summary,
  fishers_raw[
    source_positions,
    setdiff(names(fishers_raw), ".observation_id"),
    drop = FALSE
  ]
)

write.csv(
  flagged_source_records,
  file.path(OUTPUT_DIR, "38_Cook_flagged_records_for_source_verification.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. FULL MODEL AVERAGING WITH CLASSICAL AND HC3 VARIANCES
# ------------------------------------------------------------

coefficient_names <- unique(unlist(lapply(
  supported_models,
  function(model) names(coef(model))
)))

full_average_table <- function(
  models,
  weights,
  coefficient_names,
  covariance_method = c("classical", "HC3")
) {
  covariance_method <- match.arg(covariance_method)
  
  if (length(models) != length(weights)) {
    stop("The number of models and weights differs.")
  }
  
  if (abs(sum(weights) - 1) > 1e-10) {
    stop("Model weights do not sum to one.")
  }
  
  beta_matrix <- matrix(
    0,
    nrow = length(models),
    ncol = length(coefficient_names),
    dimnames = list(names(models), coefficient_names)
  )
  
  within_variance_matrix <- beta_matrix
  
  for (i in seq_along(models)) {
    model_coefficients <- coef(models[[i]])
    
    if (anyNA(model_coefficients)) {
      stop("A model contains an aliased or missing coefficient.")
    }
    
    model_covariance <- if (covariance_method == "HC3") {
      sandwich::vcovHC(models[[i]], type = "HC3")
    } else {
      vcov(models[[i]])
    }
    
    covariance_names <- rownames(model_covariance)
    covariance_diagonal <- diag(model_covariance)
    
    beta_matrix[i, names(model_coefficients)] <- model_coefficients
    within_variance_matrix[i, covariance_names] <- covariance_diagonal
  }
  
  averaged_estimate <- colSums(sweep(
    beta_matrix,
    MARGIN = 1,
    STATS = weights,
    FUN = "*"
  ))
  
  squared_between_model_difference <- sweep(
    beta_matrix,
    MARGIN = 2,
    STATS = averaged_estimate,
    FUN = "-"
  )^2
  
  unconditional_variance <- colSums(sweep(
    within_variance_matrix + squared_between_model_difference,
    MARGIN = 1,
    STATS = weights,
    FUN = "*"
  ))
  
  unconditional_se <- sqrt(unconditional_variance)
  z_value <- averaged_estimate / unconditional_se
  critical_value <- qnorm(0.975)
  
  data.frame(
    term = coefficient_names,
    estimate = as.numeric(averaged_estimate),
    unconditional_se = as.numeric(unconditional_se),
    statistic = as.numeric(z_value),
    p.value = 2 * pnorm(abs(z_value), lower.tail = FALSE),
    conf.low = as.numeric(
      averaged_estimate - critical_value * unconditional_se
    ),
    conf.high = as.numeric(
      averaged_estimate + critical_value * unconditional_se
    ),
    covariance_method = covariance_method,
    row.names = NULL,
    check.names = FALSE
  )
}

classical_recomputed <- full_average_table(
  supported_models,
  original_weights,
  coefficient_names,
  covariance_method = "classical"
)

hc3_averaged <- full_average_table(
  supported_models,
  original_weights,
  coefficient_names,
  covariance_method = "HC3"
)

write.csv(
  classical_recomputed,
  file.path(OUTPUT_DIR, "39_classical_full_model_averaged_recomputed.csv"),
  row.names = FALSE
)

write.csv(
  hc3_averaged,
  file.path(OUTPUT_DIR, "40_HC3_full_model_averaged_coefficients.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. PREPARE THE NON-PARAMETRIC PAIRS BOOTSTRAP
# ------------------------------------------------------------

# Every bootstrap sample resamples complete observations (rows), refits all
# 10 supported models, recalculates AICc weights across those 10 models, and
# then calculates full model-averaged coefficients. No observation is deleted
# merely because it was influential in the original data.

used_variables <- unique(unlist(lapply(supported_formulas, all.vars)))
used_predictors <- setdiff(used_variables, "mean_distance")
used_factor_variables <- used_predictors[vapply(
  analysis_data[, used_predictors, drop = FALSE],
  is.factor,
  logical(1)
)]

original_factor_levels <- lapply(
  analysis_data[, used_factor_variables, drop = FALSE],
  levels
)

all_required_factor_levels_present <- function(data) {
  if (length(used_factor_variables) == 0L) {
    return(TRUE)
  }
  
  all(vapply(
    used_factor_variables,
    function(variable_name) {
      observed_levels <- unique(as.character(data[[variable_name]]))
      all(original_factor_levels[[variable_name]] %in% observed_levels)
    },
    logical(1)
  ))
}

full_average_point_estimate <- function(
  models,
  weights,
  coefficient_names
) {
  beta_matrix <- matrix(
    0,
    nrow = length(models),
    ncol = length(coefficient_names),
    dimnames = list(names(models), coefficient_names)
  )
  
  for (i in seq_along(models)) {
    model_coefficients <- coef(models[[i]])
    
    if (anyNA(model_coefficients)) {
      stop("Missing bootstrap coefficient.")
    }
    
    unexpected_names <- setdiff(
      names(model_coefficients),
      coefficient_names
    )
    
    if (length(unexpected_names) > 0L) {
      stop(
        "Unexpected bootstrap coefficient(s): ",
        paste(unexpected_names, collapse = ", ")
      )
    }
    
    beta_matrix[i, names(model_coefficients)] <- model_coefficients
  }
  
  as.numeric(colSums(sweep(
    beta_matrix,
    MARGIN = 1,
    STATS = weights,
    FUN = "*"
  )))
}

bootstrap_once <- function() {
  sampled_rows <- sample.int(
    n = nrow(analysis_data),
    size = nrow(analysis_data),
    replace = TRUE
  )
  
  bootstrap_data <- analysis_data[sampled_rows, , drop = FALSE]
  
  if (!all_required_factor_levels_present(bootstrap_data)) {
    return(list(
      success = FALSE,
      reason = "at_least_one_required_factor_level_absent"
    ))
  }
  
  bootstrap_models <- tryCatch(
    fit_supported_models(bootstrap_data),
    error = function(e) e
  )
  
  if (inherits(bootstrap_models, "error")) {
    return(list(
      success = FALSE,
      reason = paste0("model_fit_error: ", conditionMessage(bootstrap_models))
    ))
  }
  
  rank_problem <- any(vapply(
    bootstrap_models,
    function(model) {
      anyNA(coef(model)) || model$rank < ncol(model.matrix(model))
    },
    logical(1)
  ))
  
  if (rank_problem) {
    return(list(
      success = FALSE,
      reason = "rank_deficient_bootstrap_model"
    ))
  }
  
  bootstrap_aicc <- tryCatch(
    vapply(bootstrap_models, aicc_lm, numeric(1)),
    error = function(e) e
  )
  
  if (inherits(bootstrap_aicc, "error") ||
      any(!is.finite(bootstrap_aicc))) {
    return(list(
      success = FALSE,
      reason = "non_finite_bootstrap_AICc"
    ))
  }
  
  bootstrap_weights <- aicc_weights(bootstrap_aicc)
  
  averaged_coefficients <- tryCatch(
    full_average_point_estimate(
      bootstrap_models,
      bootstrap_weights,
      coefficient_names
    ),
    error = function(e) e
  )
  
  if (inherits(averaged_coefficients, "error") ||
      any(!is.finite(averaged_coefficients))) {
    return(list(
      success = FALSE,
      reason = "non_finite_model_averaged_coefficient"
    ))
  }
  
  list(
    success = TRUE,
    coefficients = averaged_coefficients,
    weights = bootstrap_weights
  )
}


# ------------------------------------------------------------
# 7. RUN 2,000 SUCCESSFUL BOOTSTRAP REPETITIONS
# ------------------------------------------------------------

if (!is.numeric(N_BOOTSTRAP) || length(N_BOOTSTRAP) != 1L ||
    !is.finite(N_BOOTSTRAP) || N_BOOTSTRAP < 100L) {
  stop("N_BOOTSTRAP must be a single finite integer of at least 100.")
}

N_BOOTSTRAP <- as.integer(N_BOOTSTRAP)
set.seed(RANDOM_SEED)

bootstrap_coefficients <- matrix(
  NA_real_,
  nrow = N_BOOTSTRAP,
  ncol = length(coefficient_names),
  dimnames = list(NULL, coefficient_names)
)

bootstrap_weights_matrix <- matrix(
  NA_real_,
  nrow = N_BOOTSTRAP,
  ncol = length(supported_models),
  dimnames = list(NULL, model_labels)
)

successful_repetitions <- 0L
total_attempts <- 0L
failure_reasons <- character(0)

message(
  "Beginning ", N_BOOTSTRAP,
  " successful pairs-bootstrap repetitions."
)

while (successful_repetitions < N_BOOTSTRAP &&
       total_attempts < MAX_BOOTSTRAP_ATTEMPTS) {
  total_attempts <- total_attempts + 1L
  bootstrap_result <- bootstrap_once()
  
  if (isTRUE(bootstrap_result$success)) {
    successful_repetitions <- successful_repetitions + 1L
    bootstrap_coefficients[successful_repetitions, ] <-
      bootstrap_result$coefficients
    bootstrap_weights_matrix[successful_repetitions, ] <-
      bootstrap_result$weights
    
    if (successful_repetitions %% PROGRESS_EVERY == 0L ||
        successful_repetitions == N_BOOTSTRAP) {
      message(
        "Completed ", successful_repetitions, "/", N_BOOTSTRAP,
        " successful repetitions (", total_attempts, " total attempts)."
      )
    }
  } else {
    failure_reasons <- c(failure_reasons, bootstrap_result$reason)
  }
}

if (successful_repetitions < N_BOOTSTRAP) {
  stop(
    "Only ", successful_repetitions,
    " successful bootstrap repetitions were obtained after ",
    total_attempts, " attempts. Increase MAX_BOOTSTRAP_ATTEMPTS only after ",
    "reviewing the failure reasons."
  )
}

# Save the raw bootstrap arrays so every summary can be independently checked.
saveRDS(
  list(
    coefficients = bootstrap_coefficients,
    weights = bootstrap_weights_matrix,
    random_seed = RANDOM_SEED,
    successful_repetitions = successful_repetitions,
    total_attempts = total_attempts,
    failure_reasons = failure_reasons
  ),
  file.path(OUTPUT_DIR, "41_bootstrap_raw_results.rds")
)


# ------------------------------------------------------------
# 8. SUMMARISE BOOTSTRAP COEFFICIENTS AND MODEL WEIGHTS
# ------------------------------------------------------------

original_point_estimate <- setNames(
  classical_recomputed$estimate,
  classical_recomputed$term
)[coefficient_names]

bootstrap_two_sided_p <- function(x) {
  lower_probability <- (sum(x <= 0) + 1) / (length(x) + 1)
  upper_probability <- (sum(x >= 0) + 1) / (length(x) + 1)
  min(1, 2 * min(lower_probability, upper_probability))
}

bootstrap_coefficient_summary <- do.call(
  rbind,
  lapply(
    seq_along(coefficient_names),
    function(j) {
      values <- bootstrap_coefficients[, j]
      percentile_interval <- quantile(
        values,
        probs = c(0.025, 0.975),
        names = FALSE,
        type = 7
      )
      
      data.frame(
        term = coefficient_names[j],
        original_point_estimate = original_point_estimate[j],
        bootstrap_mean = mean(values),
        bootstrap_median = median(values),
        bootstrap_standard_error = sd(values),
        bootstrap_bias = mean(values) - original_point_estimate[j],
        bootstrap_conf_low_percentile = percentile_interval[1L],
        bootstrap_conf_high_percentile = percentile_interval[2L],
        bootstrap_p_two_sided = bootstrap_two_sided_p(values),
        successful_repetitions = length(values),
        row.names = NULL,
        check.names = FALSE
      )
    }
  )
)

top_model_index <- max.col(bootstrap_weights_matrix, ties.method = "first")

bootstrap_weight_summary <- do.call(
  rbind,
  lapply(
    seq_along(model_labels),
    function(j) {
      values <- bootstrap_weights_matrix[, j]
      intervals <- quantile(
        values,
        probs = c(0.025, 0.5, 0.975),
        names = FALSE,
        type = 7
      )
      
      data.frame(
        model = model_labels[j],
        model_id = model_ids[j],
        original_weight = original_weights[j],
        bootstrap_mean_weight = mean(values),
        bootstrap_sd_weight = sd(values),
        bootstrap_weight_2.5_percent = intervals[1L],
        bootstrap_median_weight = intervals[2L],
        bootstrap_weight_97.5_percent = intervals[3L],
        proportion_repetitions_top_weight = mean(top_model_index == j),
        row.names = NULL,
        check.names = FALSE
      )
    }
  )
)

write.csv(
  bootstrap_coefficient_summary,
  file.path(
    OUTPUT_DIR,
    "42_pairs_bootstrap_full_model_averaged_coefficients.csv"
  ),
  row.names = FALSE
)

write.csv(
  bootstrap_weight_summary,
  file.path(OUTPUT_DIR, "43_pairs_bootstrap_model_weight_summary.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. COMPARE ORIGINAL, HC3 AND BOOTSTRAP INFERENCE
# ------------------------------------------------------------

if (!is.na(original_coefficient_file)) {
  original_reported <- read.csv(
    original_coefficient_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  required_original_columns <- c(
    "term", "estimate", "unconditional_se", "p.value",
    "conf.low", "conf.high"
  )
  
  if (!all(required_original_columns %in% names(original_reported))) {
    stop(
      "The optional original coefficient file does not contain the expected ",
      "columns."
    )
  }
  
  matched_original_estimates <- original_reported$estimate[
    match(classical_recomputed$term, original_reported$term)
  ]
  
  maximum_original_estimate_difference <- max(
    abs(classical_recomputed$estimate - matched_original_estimates),
    na.rm = TRUE
  )
  
  if (!is.finite(maximum_original_estimate_difference) ||
      maximum_original_estimate_difference > 1e-5) {
    stop(
      "The recomputed full-average estimates do not reproduce file 10. ",
      "Maximum estimate difference = ",
      signif(maximum_original_estimate_difference, 8), "."
    )
  }
} else {
  original_reported <- classical_recomputed[, c(
    "term", "estimate", "unconditional_se", "p.value",
    "conf.low", "conf.high"
  )]
  maximum_original_estimate_difference <- 0
}

term_order <- original_reported$term

hc3_match <- match(term_order, hc3_averaged$term)
bootstrap_match <- match(
  term_order,
  bootstrap_coefficient_summary$term
)

comparison_table <- data.frame(
  term = term_order,
  estimate_original = original_reported$estimate,
  se_original = original_reported$unconditional_se,
  p_original = original_reported$p.value,
  conf_low_original = original_reported$conf.low,
  conf_high_original = original_reported$conf.high,
  estimate_HC3 = hc3_averaged$estimate[hc3_match],
  se_HC3 = hc3_averaged$unconditional_se[hc3_match],
  p_HC3 = hc3_averaged$p.value[hc3_match],
  conf_low_HC3 = hc3_averaged$conf.low[hc3_match],
  conf_high_HC3 = hc3_averaged$conf.high[hc3_match],
  bootstrap_mean = bootstrap_coefficient_summary$bootstrap_mean[
    bootstrap_match
  ],
  bootstrap_median = bootstrap_coefficient_summary$bootstrap_median[
    bootstrap_match
  ],
  bootstrap_se = bootstrap_coefficient_summary$bootstrap_standard_error[
    bootstrap_match
  ],
  bootstrap_conf_low =
    bootstrap_coefficient_summary$bootstrap_conf_low_percentile[
      bootstrap_match
    ],
  bootstrap_conf_high =
    bootstrap_coefficient_summary$bootstrap_conf_high_percentile[
      bootstrap_match
    ],
  bootstrap_p_two_sided =
    bootstrap_coefficient_summary$bootstrap_p_two_sided[
      bootstrap_match
    ],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

comparison_table$original_CI_crosses_zero <- with(
  comparison_table,
  conf_low_original <= 0 & conf_high_original >= 0
)
comparison_table$HC3_CI_crosses_zero <- with(
  comparison_table,
  conf_low_HC3 <= 0 & conf_high_HC3 >= 0
)
comparison_table$bootstrap_CI_crosses_zero <- with(
  comparison_table,
  bootstrap_conf_low <= 0 & bootstrap_conf_high >= 0
)
comparison_table$HC3_conclusion_changed <- with(
  comparison_table,
  original_CI_crosses_zero != HC3_CI_crosses_zero
)
comparison_table$bootstrap_conclusion_changed <- with(
  comparison_table,
  original_CI_crosses_zero != bootstrap_CI_crosses_zero
)

write.csv(
  comparison_table,
  file.path(
    OUTPUT_DIR,
    "44_original_vs_HC3_vs_bootstrap_coefficients.csv"
  ),
  row.names = FALSE
)

# A separate cost-effect table reports the coefficient in the manuscript's
# interpretable unit: change in fishing distance per USD 1,000.
cost_row <- comparison_table[
  comparison_table$term == "cost_per_trip_usd",
  ,
  drop = FALSE
]

if (nrow(cost_row) == 1L) {
  cost_effect_per_1000 <- data.frame(
    method = c("Original", "HC3", "Pairs bootstrap"),
    estimate_km_per_USD_1000 = 1000 * c(
      cost_row$estimate_original,
      cost_row$estimate_HC3,
      cost_row$bootstrap_mean
    ),
    confidence_low_km_per_USD_1000 = 1000 * c(
      cost_row$conf_low_original,
      cost_row$conf_low_HC3,
      cost_row$bootstrap_conf_low
    ),
    confidence_high_km_per_USD_1000 = 1000 * c(
      cost_row$conf_high_original,
      cost_row$conf_high_HC3,
      cost_row$bootstrap_conf_high
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  write.csv(
    cost_effect_per_1000,
    file.path(OUTPUT_DIR, "45_cost_effect_per_USD_1000.csv"),
    row.names = FALSE
  )
}


# ------------------------------------------------------------
# 10. COMPLETION SUMMARY, FAILURE LOG AND REPRODUCIBILITY FILES
# ------------------------------------------------------------

non_intercept <- comparison_table$term != "(Intercept)"

robustness_summary <- data.frame(
  status = "Completed",
  primary_dataset_n = nrow(analysis_data),
  supported_models = length(supported_models),
  model_set_rule = "Original strict Delta AICc < 2 model set",
  bootstrap_scope = paste(
    "Same 10 supported models refitted and reweighted in every",
    "pairs-bootstrap repetition"
  ),
  bootstrap_successful_repetitions = successful_repetitions,
  bootstrap_total_attempts = total_attempts,
  bootstrap_failed_attempts = total_attempts - successful_repetitions,
  random_seed = RANDOM_SEED,
  maximum_refitted_AICc_difference = maximum_aicc_difference,
  maximum_original_estimate_difference =
    maximum_original_estimate_difference,
  Cooks_threshold = cook_threshold,
  Cook_union_observations = length(cook_union_ids),
  Cook_union_ids = paste(cook_union_ids, collapse = ", "),
  HC3_nonintercept_CI_conclusions_changed = sum(
    comparison_table$HC3_conclusion_changed[non_intercept],
    na.rm = TRUE
  ),
  bootstrap_nonintercept_CI_conclusions_changed = sum(
    comparison_table$bootstrap_conclusion_changed[non_intercept],
    na.rm = TRUE
  ),
  influential_rows_automatically_removed = 0L,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(
  robustness_summary,
  file.path(OUTPUT_DIR, "46_robustness_analysis_summary.csv"),
  row.names = FALSE
)

if (length(failure_reasons) == 0L) {
  failure_summary <- data.frame(
    reason = "No failed bootstrap attempts",
    count = 0L,
    stringsAsFactors = FALSE
  )
} else {
  failure_counts <- sort(table(failure_reasons), decreasing = TRUE)
  failure_summary <- data.frame(
    reason = names(failure_counts),
    count = as.integer(failure_counts),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

write.csv(
  failure_summary,
  file.path(OUTPUT_DIR, "47_bootstrap_failure_log.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(OUTPUT_DIR, "48_R_session_information.txt")
)

writeLines(
  c(
    "FISHING-DISTANCE HC3 AND BOOTSTRAP ROBUSTNESS ANALYSIS",
    "",
    paste("Primary sample size:", nrow(analysis_data)),
    paste("Supported models:", length(supported_models)),
    paste("Successful bootstrap repetitions:", successful_repetitions),
    paste("Total bootstrap attempts:", total_attempts),
    paste("Random seed:", RANDOM_SEED),
    paste("Cook's-distance threshold:", signif(cook_threshold, 8)),
    paste("Cook-flagged analysis rows:", paste(cook_union_ids, collapse = ", ")),
    "",
    "Interpretation:",
    "- The primary analysis retains all 148 valid observations.",
    "- File 38 is for source-data verification only.",
    "- Do not remove a flagged record unless the source check confirms an error.",
    "- HC3 changes the within-model covariance estimates, not the coefficients.",
    "- The pairs bootstrap resamples complete observations and recalculates",
    "  AICc weights among the same 10 verified supported models.",
    "- Bootstrap percentile intervals are the main distribution-robust check.",
    "- This analysis is conditional on the original 10-model set; it does not",
    "  rerun the 65,536-model all-subsets search in every bootstrap sample.",
    "- Compare zero-crossing conclusions in file 44 before revising the paper."
  ),
  file.path(OUTPUT_DIR, "49_READ_ME_robustness_outputs.txt")
)

mandatory_outputs <- c(
  "37_robustness_model_refit_verification.csv",
  "38_Cook_flagged_records_for_source_verification.csv",
  "39_classical_full_model_averaged_recomputed.csv",
  "40_HC3_full_model_averaged_coefficients.csv",
  "41_bootstrap_raw_results.rds",
  "42_pairs_bootstrap_full_model_averaged_coefficients.csv",
  "43_pairs_bootstrap_model_weight_summary.csv",
  "44_original_vs_HC3_vs_bootstrap_coefficients.csv",
  "46_robustness_analysis_summary.csv",
  "47_bootstrap_failure_log.csv",
  "48_R_session_information.txt",
  "49_READ_ME_robustness_outputs.txt"
)

missing_outputs <- mandatory_outputs[
  !file.exists(file.path(OUTPUT_DIR, mandatory_outputs))
]

if (length(missing_outputs) > 0L) {
  stop(
    "The analysis ended without creating: ",
    paste(missing_outputs, collapse = ", ")
  )
}

message("")
message("ROBUSTNESS ANALYSIS COMPLETED SUCCESSFULLY.")
message("Outputs: ", normalizePath(OUTPUT_DIR, winslash = "/"))
message("")
message("First inspect:")
message("  38_Cook_flagged_records_for_source_verification.csv")
message("  44_original_vs_HC3_vs_bootstrap_coefficients.csv")
message("  46_robustness_analysis_summary.csv")




