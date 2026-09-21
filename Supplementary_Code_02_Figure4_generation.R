# =============================================================================
# SUPPLEMENTARY CODE 02
# Figure 4: model-averaged effects and adjusted fishing-distance predictions
# =============================================================================
#
# Run order
# ---------
# 1. Run Supplementary_Code_01_Fishing_distance_GLM.R first.
# 2. Keep this script in the same project directory as:
#      model_data_recoded_full.csv
#    and the output folder:
#      Fishing_distance_HC3_bootstrap_robustness/
# 3. Source this file from top to bottom.
#
# Figure panels
# -------------
# A. Full model-averaged effects for every non-intercept term represented in
#    the 10 supported models, using pairs-bootstrap 95% percentile intervals.
# B. Adjusted relationship between trip cost and average fishing distance.
# C. Adjusted average fishing distance for each departure harbour.
#
# Output:
#   Figure_4_output/Figure_4_three_panel_1000dpi.png
# plus the plotting data, bootstrap prediction object, model-selection check,
# and R session information.
# =============================================================================


# ------------------------------ USER SETTINGS -------------------------------

project_dir <- "."

n_boot <- 2000L
bootstrap_seed <- 20260919L

figure_width_in <- 9.2
figure_height_in <- 10.2
figure_dpi <- 1000

output_dir <- file.path(project_dir, "Figure_4_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------- REQUIRED PACKAGES -----------------------------

required_packages <- c(
  "dplyr",
  "ggplot2",
  "patchwork",
  "readr",
  "scales",
  "stringr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    ". Installation command: install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
})


# -------------------------------- INPUT FILES -------------------------------

locate_project_file <- function(filename, search_dirs) {
  candidates <- file.path(search_dirs, filename)
  found <- candidates[file.exists(candidates)]

  if (length(found) == 0L) {
    stop(
      "Required file not found: ", filename, "\nSearched:\n",
      paste(candidates, collapse = "\n"),
      call. = FALSE
    )
  }

  normalizePath(found[1L], winslash = "/", mustWork = TRUE)
}

data_file <- locate_project_file(
  "model_data_recoded_full.csv",
  c(project_dir)
)

coefficient_file <- locate_project_file(
  "42_pairs_bootstrap_full_model_averaged_coefficients.csv",
  c(
    project_dir,
    file.path(project_dir, "Fishing_distance_HC3_bootstrap_robustness")
  )
)

message("DATA FILE USED: ", data_file)
message("BOOTSTRAP COEFFICIENT FILE USED: ", coefficient_file)

raw_data <- read_csv(data_file, show_col_types = FALSE)
bootstrap_coefficients <- read_csv(coefficient_file, show_col_types = FALSE)

# --------------------------- RECREATE MODEL CODING --------------------------
required_source_columns <- c(
  "mean_distance",
  "Site",
  "Fishing.experience",
  "cost_per_trip_usd",
  "target.sharks",
  "communication",
  "gear_type_clean"
)

missing_columns <- setdiff(required_source_columns, names(raw_data))
if (length(missing_columns) > 0L) {
  stop(
    "The data file is missing these required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

analysis_data <- raw_data %>%
  mutate(
    Site = factor(Site, levels = c("A", "B", "C")),
    communication = factor(
      communication,
      levels = c("Every trip", "Less than every trip")
    ),
    target.sharks = factor(target.sharks, levels = c("No", "Yes")),
    uses_handline = factor(
      if_else(
        str_detect(str_to_lower(gear_type_clean), fixed("handline")),
        "Yes",
        "No"
      ),
      levels = c("No", "Yes")
    ),
    uses_longline = factor(
      if_else(
        str_detect(str_to_lower(gear_type_clean), fixed("longline")),
        "Yes",
        "No"
      ),
      levels = c("No", "Yes")
    )
  ) %>%
  filter(
    if_all(
      all_of(c(
        "mean_distance",
        "Site",
        "Fishing.experience",
        "cost_per_trip_usd",
        "target.sharks",
        "communication",
        "uses_handline",
        "uses_longline"
      )),
      ~ !is.na(.x)
    )
  )

if (nrow(analysis_data) != 148L) {
  stop(
    "Expected 148 complete observations, but obtained ",
    nrow(analysis_data),
    ". Confirm that model_data_recoded_full.csv is the original analysis file.",
    call. = FALSE
  )
}

options(contrasts = c("contr.treatment", "contr.poly"))

# These are the same 10 model specifications retained by the strict
# Delta AICc < 2 rule. Their order matches the verified analysis outputs.
supported_formulas <- list(
  M1 = mean_distance ~ cost_per_trip_usd + Site +
    uses_handline + uses_longline,
  M2 = mean_distance ~ cost_per_trip_usd + Site + uses_longline,
  M3 = mean_distance ~ cost_per_trip_usd + Site,
  M4 = mean_distance ~ cost_per_trip_usd + Site + uses_handline,
  M5 = mean_distance ~ cost_per_trip_usd + Fishing.experience +
    Site + uses_longline,
  M6 = mean_distance ~ cost_per_trip_usd + Fishing.experience +
    Site + uses_handline + uses_longline,
  M7 = mean_distance ~ communication + cost_per_trip_usd +
    Site + uses_longline,
  M8 = mean_distance ~ cost_per_trip_usd + Site +
    target.sharks + uses_longline,
  M9 = mean_distance ~ cost_per_trip_usd + Fishing.experience + Site,
  M10 = mean_distance ~ cost_per_trip_usd + Site +
    target.sharks + uses_handline + uses_longline
)

# --------------------------- MODEL-AVERAGING TOOLS --------------------------

fit_supported_models <- function(data) {
  lapply(
    supported_formulas,
    function(model_formula) {
      suppressWarnings(
        glm(
          formula = model_formula,
          data = data,
          family = gaussian(link = "identity"),
          na.action = na.fail,
          control = glm.control(maxit = 100)
        )
      )
    }
  )
}

calculate_aicc <- function(model) {
  model_loglik <- logLik(model)
  k <- attr(model_loglik, "df")
  n <- nobs(model)
  
  if (!is.finite(as.numeric(model_loglik)) || n <= (k + 1)) {
    return(Inf)
  }
  
  -2 * as.numeric(model_loglik) +
    2 * k +
    (2 * k * (k + 1)) / (n - k - 1)
}

calculate_weights <- function(models) {
  aicc <- vapply(models, calculate_aicc, numeric(1))
  delta <- aicc - min(aicc)
  raw_weight <- exp(-0.5 * delta)
  weight <- raw_weight / sum(raw_weight)
  
  list(aicc = aicc, delta = delta, weight = weight)
}

valid_model_set <- function(models) {
  all(
    vapply(
      models,
      function(model) {
        isTRUE(model$converged) &&
          model$rank == length(coef(model)) &&
          all(is.finite(coef(model))) &&
          is.finite(calculate_aicc(model))
      },
      logical(1)
    )
  )
}

# Because all supported models are additive and use an identity link, the
# adjusted cost curve can be calculated exactly from each model's average
# prediction and its cost coefficient. Other observed characteristics are
# averaged over the supplied reference data.
model_averaged_cost_prediction <- function(
  models,
  weights,
  reference_data,
  cost_values
) {
  mean_reference_cost <- mean(reference_data$cost_per_trip_usd)
  
  predictions_by_model <- vapply(
    models,
    function(model) {
      mean_prediction <- mean(
        predict(model, newdata = reference_data, type = "response")
      )
      cost_coefficient <- unname(coef(model)[["cost_per_trip_usd"]])
      
      mean_prediction +
        cost_coefficient * (cost_values - mean_reference_cost)
    },
    numeric(length(cost_values))
  )
  
  as.numeric(predictions_by_model %*% weights)
}

model_averaged_harbour_prediction <- function(
  models,
  weights,
  reference_data,
  harbour_levels = c("A", "B", "C")
) {
  predictions_by_model <- vapply(
    models,
    function(model) {
      vapply(
        harbour_levels,
        function(harbour) {
          prediction_data <- reference_data
          prediction_data$Site <- factor(
            harbour,
            levels = levels(analysis_data$Site)
          )
          mean(predict(model, newdata = prediction_data, type = "response"))
        },
        numeric(1)
      )
    },
    numeric(length(harbour_levels))
  )
  
  as.numeric(predictions_by_model %*% weights)
}

# ---------------------- FIT AND VERIFY ORIGINAL MODELS ----------------------

original_models <- fit_supported_models(analysis_data)

if (!valid_model_set(original_models)) {
  stop(
    "At least one original supported model failed to converge or was rank deficient.",
    call. = FALSE
  )
}

original_selection <- calculate_weights(original_models)

# These verified AICc values protect against accidental changes in factor
# coding, input data, or the model definitions.
expected_aicc <- c(
  2228.984434,
  2229.049674,
  2229.832813,
  2230.109028,
  2230.553155,
  2230.762744,
  2230.807756,
  2230.876547,
  2230.953918,
  2230.966333
)

maximum_aicc_difference <- max(abs(original_selection$aicc - expected_aicc))
if (maximum_aicc_difference > 0.01) {
  stop(
    "The supported models did not reproduce the verified AICc values. ",
    "Maximum absolute difference = ",
    signif(maximum_aicc_difference, 5),
    ". Check the input data and factor coding.",
    call. = FALSE
  )
}

model_selection_check <- tibble(
  model = names(supported_formulas),
  AICc = original_selection$aicc,
  delta_AICc = original_selection$delta,
  renormalised_weight = original_selection$weight
)

write_csv(
  model_selection_check,
  file.path(output_dir, "Figure_4_model_selection_check.csv")
)

# ------------------ ORIGINAL ADJUSTED PREDICTION ESTIMATES ------------------

cost_grid <- seq(
  from = min(analysis_data$cost_per_trip_usd),
  to = max(analysis_data$cost_per_trip_usd),
  length.out = 120L
)

harbour_levels <- c("A", "B", "C")

original_cost_prediction <- model_averaged_cost_prediction(
  models = original_models,
  weights = original_selection$weight,
  reference_data = analysis_data,
  cost_values = cost_grid
)

original_harbour_prediction <- model_averaged_harbour_prediction(
  models = original_models,
  weights = original_selection$weight,
  reference_data = analysis_data,
  harbour_levels = harbour_levels
)

# -------------------- PAIRS BOOTSTRAP FOR PANELS B AND C --------------------

set.seed(bootstrap_seed)

n <- nrow(analysis_data)
cost_bootstrap <- matrix(
  NA_real_,
  nrow = n_boot,
  ncol = length(cost_grid)
)
harbour_bootstrap <- matrix(
  NA_real_,
  nrow = n_boot,
  ncol = length(harbour_levels)
)

successful_repetitions <- 0L
total_attempts <- 0L
maximum_attempts <- n_boot + 500L

message("Starting ", n_boot, " pairs-bootstrap repetitions...")

while (
  successful_repetitions < n_boot &&
  total_attempts < maximum_attempts
) {
  total_attempts <- total_attempts + 1L
  sampled_rows <- sample.int(n, size = n, replace = TRUE)
  bootstrap_data <- analysis_data[sampled_rows, , drop = FALSE]
  
  bootstrap_result <- tryCatch(
    {
      bootstrap_models <- fit_supported_models(bootstrap_data)
      
      if (!valid_model_set(bootstrap_models)) {
        stop("Invalid bootstrap model set")
      }
      
      bootstrap_selection <- calculate_weights(bootstrap_models)
      
      cost_prediction <- model_averaged_cost_prediction(
        models = bootstrap_models,
        weights = bootstrap_selection$weight,
        reference_data = bootstrap_data,
        cost_values = cost_grid
      )
      
      harbour_prediction <- model_averaged_harbour_prediction(
        models = bootstrap_models,
        weights = bootstrap_selection$weight,
        reference_data = bootstrap_data,
        harbour_levels = harbour_levels
      )
      
      if (
        any(!is.finite(cost_prediction)) ||
        any(!is.finite(harbour_prediction))
      ) {
        stop("Non-finite bootstrap prediction")
      }
      
      list(
        cost = cost_prediction,
        harbour = harbour_prediction
      )
    },
    error = function(error_condition) NULL
  )
  
  if (!is.null(bootstrap_result)) {
    successful_repetitions <- successful_repetitions + 1L
    cost_bootstrap[successful_repetitions, ] <- bootstrap_result$cost
    harbour_bootstrap[successful_repetitions, ] <- bootstrap_result$harbour
    
    if (
      successful_repetitions %% 100L == 0L ||
      successful_repetitions == n_boot
    ) {
      message(
        "Completed ",
        successful_repetitions,
        " of ",
        n_boot,
        " successful repetitions"
      )
    }
  }
}

if (successful_repetitions < n_boot) {
  stop(
    "Only ",
    successful_repetitions,
    " successful bootstrap repetitions were obtained after ",
    total_attempts,
    " attempts.",
    call. = FALSE
  )
}

cost_lower <- apply(
  cost_bootstrap,
  2,
  quantile,
  probs = 0.025,
  names = FALSE,
  type = 7
)
cost_upper <- apply(
  cost_bootstrap,
  2,
  quantile,
  probs = 0.975,
  names = FALSE,
  type = 7
)

harbour_lower <- apply(
  harbour_bootstrap,
  2,
  quantile,
  probs = 0.025,
  names = FALSE,
  type = 7
)
harbour_upper <- apply(
  harbour_bootstrap,
  2,
  quantile,
  probs = 0.975,
  names = FALSE,
  type = 7
)

cost_prediction_data <- tibble(
  cost_per_trip_usd = cost_grid,
  predicted_distance_km = original_cost_prediction,
  bootstrap_conf_low = cost_lower,
  bootstrap_conf_high = cost_upper
)

harbour_prediction_data <- tibble(
  Site = factor(harbour_levels, levels = harbour_levels),
  predicted_distance_km = original_harbour_prediction,
  bootstrap_conf_low = harbour_lower,
  bootstrap_conf_high = harbour_upper
)

write_csv(
  cost_prediction_data,
  file.path(output_dir, "Figure_4_panel_B_cost_predictions.csv")
)
write_csv(
  harbour_prediction_data,
  file.path(output_dir, "Figure_4_panel_C_harbour_predictions.csv")
)

saveRDS(
  list(
    seed = bootstrap_seed,
    successful_repetitions = successful_repetitions,
    total_attempts = total_attempts,
    cost_grid = cost_grid,
    cost_predictions = cost_bootstrap,
    harbour_levels = harbour_levels,
    harbour_predictions = harbour_bootstrap
  ),
  file.path(output_dir, "Figure_4_prediction_bootstrap_results.rds")
)

# ---------------------- PANEL A: BOOTSTRAP EFFECT SIZES ----------------------

required_coefficient_columns <- c(
  "term",
  "original_point_estimate",
  "bootstrap_conf_low_percentile",
  "bootstrap_conf_high_percentile"
)

missing_coefficient_columns <- setdiff(
  required_coefficient_columns,
  names(bootstrap_coefficients)
)

if (length(missing_coefficient_columns) > 0L) {
  stop(
    "The bootstrap coefficient file is missing these columns: ",
    paste(missing_coefficient_columns, collapse = ", "),
    call. = FALSE
  )
}

effect_order <- c(
  "Harbour B (vs A)",
  "Harbour C (vs A)",
  "Cost per trip (+USD 1,000)",
  "Fishing experience (+1 year)",
  "Communication (< every trip)",
  "Target sharks (Yes vs No)",
  "Handline (Yes vs No)",
  "Longline (Yes vs No)"
)

effect_data <- bootstrap_coefficients %>%
  filter(term != "(Intercept)") %>%
  mutate(
    display_term = case_when(
      term == "SiteB" ~ "Harbour B (vs A)",
      term == "SiteC" ~ "Harbour C (vs A)",
      term == "cost_per_trip_usd" ~ "Cost per trip (+USD 1,000)",
      term == "Fishing.experience" ~ "Fishing experience (+1 year)",
      term == "communicationLess than every trip" ~
        "Communication (< every trip)",
      term == "target.sharksYes" ~ "Target sharks (Yes vs No)",
      term == "uses_handlineYes" ~ "Handline (Yes vs No)",
      term == "uses_longlineYes" ~ "Longline (Yes vs No)",
      TRUE ~ term
    ),
    scale_multiplier = if_else(
      term == "cost_per_trip_usd",
      1000,
      1
    ),
    estimate_km = original_point_estimate * scale_multiplier,
    conf_low_km = bootstrap_conf_low_percentile * scale_multiplier,
    conf_high_km = bootstrap_conf_high_percentile * scale_multiplier,
    display_term = factor(
      display_term,
      levels = rev(effect_order)
    )
  ) %>%
  arrange(display_term)

if (!setequal(as.character(effect_data$display_term), effect_order)) {
  stop(
    "Panel A does not contain exactly the expected eight non-intercept terms.",
    call. = FALSE
  )
}

write_csv(
  effect_data %>%
    select(
      term,
      display_term,
      estimate_km,
      conf_low_km,
      conf_high_km
    ),
  file.path(output_dir, "Figure_4_panel_A_effect_data.csv")
)

# ------------------------------- PLOT STYLING -------------------------------

publication_theme <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title = element_text(colour = "black", size = 11),
    axis.text = element_text(colour = "black", size = 9.5),
    axis.line = element_line(colour = "black", linewidth = 0.45),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    legend.position = "none",
    plot.margin = margin(8, 10, 8, 10),
    plot.tag = element_text(
      colour = "black",
      face = "bold",
      size = 15
    ),
    plot.tag.position = c(0, 1)
  )

# Panel A: original model-averaged point estimates with bootstrap intervals.
panel_a <- ggplot(
  effect_data,
  aes(x = estimate_km, y = display_term)
) +
  geom_vline(
    xintercept = 0,
    colour = "black",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      x = conf_low_km,
      xend = conf_high_km,
      y = display_term,
      yend = display_term
    ),
    colour = "black",
    linewidth = 0.65,
    lineend = "round"
  ) +
  geom_point(
    shape = 21,
    size = 2.8,
    stroke = 0.65,
    colour = "black",
    fill = "black"
  ) +
  scale_x_continuous(labels = label_number(big.mark = ",")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 34)) +
  labs(
    x = "Full model-averaged change in average fishing distance (km)",
    y = NULL,
    tag = "A"
  ) +
  publication_theme +
  theme(
    axis.text.y = element_text(size = 10),
    plot.margin = margin(8, 14, 10, 12)
  )

# Panels B and C share the same y-axis range so their adjusted predictions are
# directly comparable.
shared_y_range <- range(
  c(
    cost_prediction_data$bootstrap_conf_low,
    cost_prediction_data$bootstrap_conf_high,
    harbour_prediction_data$bootstrap_conf_low,
    harbour_prediction_data$bootstrap_conf_high
  ),
  finite = TRUE
)
shared_y_padding <- 0.05 * diff(shared_y_range)
shared_y_limits <- shared_y_range + c(-shared_y_padding, shared_y_padding)

panel_b <- ggplot(
  cost_prediction_data,
  aes(x = cost_per_trip_usd, y = predicted_distance_km)
) +
  geom_ribbon(
    aes(
      ymin = bootstrap_conf_low,
      ymax = bootstrap_conf_high
    ),
    fill = "grey82",
    colour = NA
  ) +
  geom_line(colour = "black", linewidth = 0.9) +
  scale_x_continuous(
    labels = label_number(big.mark = ","),
    breaks = pretty_breaks(n = 4)
  ) +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    breaks = pretty_breaks(n = 5)
  ) +
  coord_cartesian(ylim = shared_y_limits) +
  labs(
    x = "Cost per trip (USD)",
    y = "Adjusted average fishing distance (km)",
    tag = "B"
  ) +
  publication_theme

panel_c <- ggplot(
  harbour_prediction_data,
  aes(x = Site, y = predicted_distance_km)
) +
  geom_errorbar(
    aes(
      ymin = bootstrap_conf_low,
      ymax = bootstrap_conf_high
    ),
    width = 0.12,
    colour = "black",
    linewidth = 0.75
  ) +
  geom_point(
    shape = 21,
    size = 3.5,
    stroke = 0.75,
    colour = "black",
    fill = "black"
  ) +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    breaks = pretty_breaks(n = 5)
  ) +
  coord_cartesian(ylim = shared_y_limits) +
  labs(
    x = "Departure harbour",
    y = "Adjusted average fishing distance (km)",
    tag = "C"
  ) +
  publication_theme

# Landscape arrangement: Panel A occupies the left side, while Panels B and C
# are stacked on the right. The extra width given to Panel A accommodates its
# predictor labels and confidence intervals without crowding.
figure_4 <- (panel_a | (panel_b / panel_c)) +
  plot_layout(widths = c(1.55, 1))

# No title or subtitle is added anywhere in the figure.
output_png <- file.path(output_dir, "Figure_4_three_panel_1000dpi.png")

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    filename = output_png,
    plot = figure_4,
    device = ragg::agg_png,
    width = figure_width_in,
    height = figure_height_in,
    units = "in",
    dpi = figure_dpi,
    background = "white",
    limitsize = FALSE
  )
} else {
  ggsave(
    filename = output_png,
    plot = figure_4,
    device = "png",
    width = figure_width_in,
    height = figure_height_in,
    units = "in",
    dpi = figure_dpi,
    bg = "white",
    limitsize = FALSE
  )
}

message("Figure saved to: ", normalizePath(output_png, winslash = "/"))
message("Resolution: ", figure_dpi, " dpi")
message(
  "Bootstrap: ",
  successful_repetitions,
  " successful repetitions from ",
  total_attempts,
  " attempts"
)

# --------------------------- REPRODUCIBILITY RECORD --------------------------

capture.output(
  sessionInfo(),
  file = file.path(output_dir, "Figure_4_R_session_information.txt")
)

writeLines(
  c(
    "FIGURE 4 REPRODUCIBILITY NOTES",
    "",
    paste("Data file:", data_file),
    paste("Coefficient file:", coefficient_file),
    paste("Prediction bootstrap repetitions:", successful_repetitions),
    paste("Prediction bootstrap total attempts:", total_attempts),
    paste("Random seed:", bootstrap_seed),
    paste("Figure dpi:", figure_dpi)
  ),
  con = file.path(output_dir, "Figure_4_READ_ME.txt")
)
