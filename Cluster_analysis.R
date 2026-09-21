# Load required libraries
library(dplyr)
library(cluster)
library(ggplot2)
library(officer)
install.packages(c("systemfonts", "gdtools", "flextable"), type = "binary")
library(flextable)
library(ggdist)

#import csv file
fishers_raw = read.csv("model_data_recoded_full.csv")

# ==============================================================================
# PRE-PROCESSING: Aligning raw data columns with model requirements
# ==============================================================================
fishers_raw = fishers_raw %>%
  mutate(
    # Rename to correct terminology
    top_species_catch_group = Top.species.group,
    
    # Map raw bycatch to the expected clean variable
    bycatch_clean = bycatch_raw,
    
    # Generate the missing binary gear type columns from gear_type_clean
    uses_gillnet = ifelse(grepl("gillnet", gear_type_clean, ignore.case = TRUE), 1, 0),
    uses_handline = ifelse(grepl("handline", gear_type_clean, ignore.case = TRUE), 1, 0),
    uses_longline = ifelse(grepl("longline", gear_type_clean, ignore.case = TRUE), 1, 0)
  )

# ==============================================================================
# SUBSET & FACTOR CONVERSION
# ==============================================================================
fishers_clustering = fishers_raw %>%
  select(
    mean_distance,
    
    # Continuous variables
    Fishing.experience,
    vessel_length_m,
    cost_per_trip_usd,
    
    # Site
    Site,
    
    # Technology
    AIS,
    VMS,
    navigation.apps,
    
    # Fishing characteristics
    target.sharks,
    top_species_catch_group, 
    uses_gillnet,
    uses_handline,
    uses_longline,
    
    # Communication / awareness
    communication,
    reg_awareness,
    
    # Bycatch
    bycatch_clean,
    
    # Other relevant variables
    subsidy
  )

#ensure all categorical variables are factors
fishers_clustering = fishers_clustering %>%
  mutate(
    Site = as.factor(Site),
    AIS = as.factor(AIS),
    VMS = as.factor(VMS),
    navigation.apps = as.factor(navigation.apps),
    target.sharks = as.factor(target.sharks),
    top_species_catch_group = as.factor(top_species_catch_group),
    uses_gillnet = as.factor(uses_gillnet),
    uses_handline = as.factor(uses_handline),
    uses_longline = as.factor(uses_longline),
    communication = as.factor(communication),
    reg_awareness = as.factor(reg_awareness),
    bycatch_clean = as.factor(bycatch_clean),
    subsidy = as.factor(subsidy)
  )

#check data structure
str(fishers_clustering)

# ==============================================================================
# GOWER DISTANCE & SILHOUETTE ANALYSIS
# ==============================================================================
#calculate gower distance
gower_dist = daisy(
  fishers_clustering,
  metric = "gower"
)
#silhouette analysis to identify optimal no of clusters
silhouette_results = data.frame(
  k = 2:10,
  silhouette = sapply(
    2:10,
    function(k) {
      
      pam_fit <- pam(
        gower_dist,
        k = k,
        diss = TRUE
      )
      
      mean(
        silhouette(
          pam_fit$clustering,
          gower_dist
        )[, 3]
      )
    }
  )
)

silhouette_results

#plot silhouette results
ggplot(
  silhouette_results,
  aes(
    x = k,
    y = silhouette
  )
) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = 2:10
  ) +
  labs(
    x = "Number of clusters (k)",
    y = "Average silhouette width",
    title = "Silhouette analysis for selecting the number of clusters"
  ) +
  theme_classic()

#fit two cluster model
pam_2_clust = pam(
  gower_dist,
  k = 2,
  diss = TRUE
)
#add clusters to data
fishers_clustering_2clust = fishers_clustering %>%
  mutate(
    cluster_2 = factor(
      pam_2_clust$clustering,
      levels = c(1, 2),
      labels = c("Cluster 1", "Cluster 2")
    )
  )
#check cluster size
fishers_clustering_2clust %>%
  dplyr::count(cluster_2) %>%
  mutate(
    percentage = round(100 * n / sum(n),2)
  )

#fit three cluster model
pam_3_clust = pam(
  gower_dist,
  k = 3,
  diss = TRUE
)
#add clusters to data
fishers_clustering_3clust = fishers_clustering %>%
  mutate(
    cluster_3 = factor(
      pam_3_clust$clustering,
      levels = c(1, 2,3),
      labels = c("Cluster 1", "Cluster 2", "Cluster 3")
    )
  )
#check cluster size
fishers_clustering_3clust %>%
  dplyr::count(cluster_3) %>%
  mutate(
    percentage = round(100 * n / sum(n),2)
  )

#check silhouette quality of 2 cluster model
sil_2_clust = silhouette(
  pam_2_clust$clustering,
  gower_dist
)

summary(sil_2_clust[, "sil_width"])

#check silhouette quality of 3 clust
sil_3_clust = silhouette(
  pam_3_clust$clustering,
  gower_dist
)

summary(sil_3_clust[, "sil_width"])

#extract widths
sil_width_2 <- sil_2_clust[, "sil_width"]
sil_width_3 <- sil_3_clust[, "sil_width"]

#create a function to summarize silhouette
summarize_silhouette <- function(x, k) {
  
  tibble(
    clusters = paste0(k, " clusters"),
    min_silhouette = min(x),
    q1_silhouette = quantile(x, 0.25),
    median_silhouette = median(x),
    mean_silhouette = mean(x),
    q3_silhouette = quantile(x, 0.75),
    max_silhouette = max(x),
    n_negative = sum(x < 0),
    pct_negative = 100 * mean(x < 0)
  )
}

#create table comparing silhouette in 2/3 cluster
silhouette_comparison <- bind_rows(
  
  summarize_silhouette(
    sil_width_2,
    2
  ),
  
  summarize_silhouette(
    sil_width_3,
    3
  )
) %>%
  
  select(
    clusters,
    mean_silhouette,
    median_silhouette,
    n_negative,
    pct_negative
  ) %>%
  
  mutate(
    mean_silhouette = round(
      mean_silhouette,
      3
    ),
    
    median_silhouette = round(
      median_silhouette,
      3
    ),
    
    pct_negative = round(
      pct_negative,
      2
    )
  )

silhouette_comparison

# Create flextable
silhouette_word <- flextable(
  silhouette_comparison
)

# Rename headers
silhouette_word <- set_header_labels(
  silhouette_word,
  clusters = "Cluster solution",
  mean_silhouette = "Mean silhouette",
  median_silhouette = "Median silhouette",
  n_negative = "Negative silhouette (n)",
  pct_negative = "Negative silhouette (%)"
)

# Professional style
silhouette_word <- theme_booktabs(
  silhouette_word
)

# Bold header
silhouette_word <- bold(
  silhouette_word,
  part = "header"
)

# Align cluster solution left
silhouette_word <- align(
  silhouette_word,
  j = "clusters",
  align = "left",
  part = "all"
)

# Center numerical columns
silhouette_word <- align(
  silhouette_word,
  j = c(
    "mean_silhouette",
    "median_silhouette",
    "n_negative",
    "pct_negative"
  ),
  align = "center",
  part = "all"
)

# Autofit
silhouette_word <- autofit(
  silhouette_word
)

doc <- read_docx()

doc <- body_add_par(
  doc,
  "Table X. Comparison of Silhouette Quality for Two- and Three-Cluster Solutions",
  style = "heading 1"
)

doc <- body_add_flextable(
  doc,
  value = silhouette_word
)

doc <- body_add_par(
  doc,
  paste0(
    "Note. Silhouette widths range from -1 to 1, with higher values ",
    "indicating better within-cluster cohesion and between-cluster separation. ",
    "Negative silhouette widths indicate observations that may be assigned ",
    "to an inappropriate cluster based on their relative distances to clusters."
  )
)

print(
  doc,
  target = "Silhouette_Comparison_2_vs_3_Clusters.docx"
)

# ==============================================================================
# CONTINUOUS CHARACTERISTICS
# ==============================================================================
#compare continous variables in 2 clust
continous_profile_2 = fishers_clustering_2clust %>%
  group_by(cluster_2) %>%
  summarise(
    n = n(),
    mean_distance = round(mean(
      mean_distance,
      na.rm = TRUE
    ),2),
    Fishing.experience = round(mean(
      Fishing.experience,
      na.rm = TRUE
    ),2),
    vessel_length_m = round(mean(
      vessel_length_m,
      na.rm = TRUE
    ),2),
    cost_per_trip_usd = round(mean(
      cost_per_trip_usd,
      na.rm = TRUE
    ),2)
  ) %>%
  dplyr::rename(
    Cluster = cluster_2,
    `Sample size` = n,
    `Mean distance` = mean_distance,
    `Fishing experience` = Fishing.experience,
    `Vessel length (m)` = vessel_length_m,
    `Cost per trip (USD)` = cost_per_trip_usd
  )

continous_profile_2

#compare continous variables in 3 cluster
continous_profile_3 = fishers_clustering_3clust %>%
  group_by(cluster_3) %>%
  summarise(
    n = n(),
    mean_distance = round(mean(
      mean_distance,
      na.rm = TRUE
    ),2),
    Fishing.experience = round(mean(
      Fishing.experience,
      na.rm = TRUE
    ),2),
    vessel_length_m = round(mean(
      vessel_length_m,
      na.rm = TRUE
    ),2),
    cost_per_trip_usd = round(mean(
      cost_per_trip_usd,
      na.rm = TRUE
    ),2)
  ) %>%
  dplyr::rename(
    Cluster = cluster_3,
    `Sample size` = n,
    `Mean distance` = mean_distance,
    `Fishing experience` = Fishing.experience,
    `Vessel length (m)` = vessel_length_m,
    `Cost per trip (USD)` = cost_per_trip_usd
  )

continous_profile_3

#export to word
profile_cont_word_2 <- flextable(
  continous_profile_2
)

# Professional style
profile_cont_word_2 <- theme_booktabs(
  profile_cont_word_2
)

# Bold header
profile_cont_word_2 <- bold(
  profile_cont_word_2,
  part = "header"
)

# Center numerical columns
profile_cont_word_2 <- align(
  profile_cont_word_2,
  j = c(
    "Sample size",
    "Mean distance",
    "Fishing experience",
    "Vessel length (m)",
    "Cost per trip (USD)"
  ),
  align = "center",
  part = "all"
)

# Left-align cluster names
profile_cont_word_2 <- align(
  profile_cont_word_2,
  j = "Cluster",
  align = "left",
  part = "all"
)

# Autofit
profile_cont_word_2 <- autofit(
  profile_cont_word_2
)

#export to word
doc <- read_docx()

doc <- body_add_par(
  doc,
  "Table X. Continuous Characteristics of Fisher Clusters",
  style = "heading 1"
)

doc <- body_add_flextable(
  doc,
  value = profile_cont_word_2
)

doc <- body_add_par(
  doc,
  paste0(
    "Note. Values are cluster-specific means, except sample size. ",
    "Clusters were identified using the selected two-cluster PAM solution."
  )
)

print(
  doc,
  target = "2_Cluster_Profile_Continuous.docx"
)

#export to word
profile_cont_word_3 <- flextable(
  continous_profile_3
)

# Professional style
profile_cont_word_3 <- theme_booktabs(
  profile_cont_word_3
)

# Bold header
profile_cont_word_3 <- bold(
  profile_cont_word_3,
  part = "header"
)

# Center numerical columns
profile_cont_word_3 <- align(
  profile_cont_word_3,
  j = c(
    "Sample size",
    "Mean distance",
    "Fishing experience",
    "Vessel length (m)",
    "Cost per trip (USD)"
  ),
  align = "center",
  part = "all"
)

# Left-align cluster names
profile_cont_word_3 <- align(
  profile_cont_word_3,
  j = "Cluster",
  align = "left",
  part = "all"
)

# Autofit
profile_cont_word_3 <- autofit(
  profile_cont_word_3
)

#export to word
doc <- read_docx()

doc <- body_add_par(
  doc,
  "Table X. Continuous Characteristics of Fisher Clusters",
  style = "heading 1"
)

doc <- body_add_flextable(
  doc,
  value = profile_cont_word_3
)

doc <- body_add_par(
  doc,
  paste0(
    "Note. Values are cluster-specific means, except sample size. ",
    "Clusters were identified using the selected three-cluster PAM solution."
  )
)

print(
  doc,
  target = "3_Cluster_Profile_Continuous.docx"
)

# ==============================================================================
# CATEGORICAL CHARACTERISTICS
# ==============================================================================
#create single categorical profile for 2 cluster
profile_categorical_2 <- fishers_clustering_2clust %>%
  dplyr::select(
    cluster_2,
    Site,
    AIS,
    VMS,
    navigation.apps,
    top_species_catch_group,
    target.sharks,
    uses_gillnet,
    uses_handline,
    uses_longline,
    subsidy,
    communication,
    reg_awareness
  ) %>%
  
  tidyr::pivot_longer(
    cols = -cluster_2,
    names_to = "Variable",
    values_to = "Category"
  ) %>%
  
  dplyr::count(
    Variable,
    Category,
    cluster_2
  ) %>%
  
  dplyr::group_by(
    Variable,
    cluster_2
  ) %>%
  
  dplyr::mutate(
    Percentage = 100 * n / sum(n)
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    `n (%)` = sprintf(
      "%d (%.1f%%)",
      n,
      Percentage
    )
  ) %>%
  
  dplyr::select(
    Variable,
    Category,
    cluster_2,
    `n (%)`
  ) %>%
  
  tidyr::pivot_wider(
    names_from = cluster_2,
    values_from = `n (%)`
  )
#formatting
profile_cat_word_2 <- flextable(
  profile_categorical_2
)

# Professional style
profile_cat_word_2 <- theme_booktabs(
  profile_cat_word_2
)

# Bold header
profile_cat_word_2 <- bold(
  profile_cat_word_2,
  part = "header"
)

# Center cluster columns
profile_cat_word_2 <- align(
  profile_cat_word_2,
  j = c(
    "Cluster 1",
    "Cluster 2"
  ),
  align = "center",
  part = "all"
)

# Left-align variable and category
profile_cat_word_2 <- align(
  profile_cat_word_2,
  j = c(
    "Variable",
    "Category"
  ),
  align = "left",
  part = "all"
)

# Autofit
profile_cat_word_2 <- autofit(
  profile_cat_word_2
)

#export
doc <- read_docx()

doc <- body_add_par(
  doc,
  "Table X. Categorical Characteristics of Fisher Clusters",
  style = "heading 1"
)

doc <- body_add_flextable(
  doc,
  value = profile_cat_word_2
)

doc <- body_add_par(
  doc,
  paste0(
    "Note. Values are presented as frequency and percentage, n (%), ",
    "within each cluster."
  )
)

print(
  doc,
  target = "2_Cluster_Profile_Categorical.docx"
)

#create single categorical profile for 3 cluster
profile_categorical_3 <- fishers_clustering_3clust %>%
  dplyr::select(
    cluster_3,
    Site,
    AIS,
    VMS,
    navigation.apps,
    top_species_catch_group,
    target.sharks,
    uses_gillnet,
    uses_handline,
    uses_longline,
    subsidy,
    communication,
    reg_awareness
  ) %>%
  
  tidyr::pivot_longer(
    cols = -cluster_3,
    names_to = "Variable",
    values_to = "Category"
  ) %>%
  
  dplyr::count(
    Variable,
    Category,
    cluster_3
  ) %>%
  
  dplyr::group_by(
    Variable,
    cluster_3
  ) %>%
  
  dplyr::mutate(
    Percentage = 100 * n / sum(n)
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::mutate(
    `n (%)` = sprintf(
      "%d (%.1f%%)",
      n,
      Percentage
    )
  ) %>%
  
  dplyr::select(
    Variable,
    Category,
    cluster_3,
    `n (%)`
  ) %>%
  
  tidyr::pivot_wider(
    names_from = cluster_3,
    values_from = `n (%)`
  )
#formatting
profile_cat_word_3 <- flextable(
  profile_categorical_3
)

# Professional style
profile_cat_word_3 <- theme_booktabs(
  profile_cat_word_3
)

# Bold header
profile_cat_word_3 <- bold(
  profile_cat_word_3,
  part = "header"
)

# Center cluster columns
profile_cat_word_3 <- align(
  profile_cat_word_3,
  j = c(
    "Cluster 1",
    "Cluster 2",
    "Cluster 3"
  ),
  align = "center",
  part = "all"
)

# Left-align variable and category
profile_cat_word_3 <- align(
  profile_cat_word_3,
  j = c(
    "Variable",
    "Category"
  ),
  align = "left",
  part = "all"
)

# Autofit
profile_cat_word_3 <- autofit(
  profile_cat_word_3
)

#export
doc <- read_docx()

doc <- body_add_par(
  doc,
  "Table X. Categorical Characteristics of Fisher Clusters",
  style = "heading 1"
)

doc <- body_add_flextable(
  doc,
  value = profile_cat_word_3
)

doc <- body_add_par(
  doc,
  paste0(
    "Note. Values are presented as frequency and percentage, n (%), ",
    "within each cluster."
  )
)

print(
  doc,
  target = "3_Cluster_Profile_Categorical.docx"
)

# ==============================================================================
# BOOTSTRAPPING STABILITY
# ==============================================================================
#custom bootstrapping function
bootstrap_pam_stability <- function(
    data,
    original_clusters,
    k,
    B = 1000,
    seed = 2026
) {
  
  set.seed(seed)
  
  n <- nrow(data)
  
  # Store Jaccard values for each original cluster
  jaccard_results <- vector(
    "list",
    k
  )
  
  for (j in 1:k) {
    jaccard_results[[j]] <- numeric(B)
  }
  
  # Bootstrap iterations
  for (b in 1:B) {
    
    # Bootstrap sample rows
    boot_index <- sample(
      1:n,
      size = n,
      replace = TRUE
    )
    
    boot_data <- data[boot_index, , drop = FALSE]
    
    # Calculate Gower distance
    boot_gower <- daisy(
      boot_data,
      metric = "gower"
    )
    
    # PAM clustering
    boot_pam <- pam(
      boot_gower,
      k = k,
      diss = TRUE
    )
    
    boot_clusters <- boot_pam$clustering
    
    # Map bootstrap observations back to original observations
    original_for_boot <- original_clusters[boot_index]
    
    # Calculate Jaccard for each original cluster
    for (j in 1:k) {
      
      # Original cluster membership
      original_members <- (
        original_for_boot == j
      )
      
      # Bootstrap cluster membership
      for (l in 1:k) {
        
        boot_members <- (
          boot_clusters == l
        )
        
        # Jaccard similarity
        intersection <- sum(
          original_members & boot_members
        )
        
        union <- sum(
          original_members | boot_members
        )
        
        if (union > 0) {
          
          jaccard_value <- (
            intersection / union
          )
          
          # Keep the best matching bootstrap cluster
          jaccard_results[[j]][b] <- max(
            jaccard_results[[j]][b],
            jaccard_value
          )
        }
      }
    }
  }
  
  # Summarize
  stability_table <- data.frame(
    cluster = paste(
      "Cluster",
      1:k
    ),
    
    mean_jaccard = sapply(
      jaccard_results,
      mean,
      na.rm = TRUE
    ),
    
    median_jaccard = sapply(
      jaccard_results,
      median,
      na.rm = TRUE
    ),
    
    lower_95 = sapply(
      jaccard_results,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    upper_95 = sapply(
      jaccard_results,
      quantile,
      probs = 0.975,
      na.rm = TRUE
    )
  )
  
  return(
    list(
      stability = stability_table,
      raw_jaccard = jaccard_results
    )
  )
}

#bootstrap 2 cluster model using 1000 iterations
set.seed(2026)

stability_2 <- bootstrap_pam_stability(
  data = fishers_clustering,
  original_clusters = pam_2_clust$clustering,
  k = 2,
  B = 1000,
  seed = 2026
)
#show bootstrap results
stability_2$stability

#bootstrap 3 cluster model using 1000 iterations
set.seed(2026)

stability_3 <- bootstrap_pam_stability(
  data = fishers_clustering,
  original_clusters = pam_3_clust$clustering,
  k = 3,
  B = 1000,
  seed = 2026
)
#show results
stability_3$stability

#compare 2/3 cluster stability
stability_comparison <- bind_rows(
  
  stability_2$stability %>%
    mutate(
      solution = "2 clusters"
    ),
  
  stability_3$stability %>%
    mutate(
      solution = "3 clusters"
    )
  
) %>%
  
  select(
    solution,
    cluster,
    mean_jaccard,
    median_jaccard,
    lower_95,
    upper_95
  ) %>%
  mutate(
    across(
      c(
        mean_jaccard,
        median_jaccard,
        lower_95,
        upper_95
      ),
      ~ round(.x, 3)
    )
  )

stability_comparison

#export to word
# Create flextable
stability_word <- flextable(
  stability_comparison
)

# Rename headers
stability_word <- set_header_labels(
  stability_word,
  solution = "Cluster solution",
  cluster = "Cluster",
  mean_jaccard = "Mean Jaccard",
  median_jaccard = "Median Jaccard",
  lower_95 = "95% CI Lower",
  upper_95 = "95% CI Upper"
)

# Apply professional table style
stability_word <- theme_booktabs(
  stability_word
)

# Bold header
stability_word <- bold(
  stability_word,
  part = "header"
)

# Left-align solution and cluster
stability_word <- align(
  stability_word,
  j = c(
    "solution",
    "cluster"
  ),
  align = "left",
  part = "all"
)

# Center numerical columns
stability_word <- align(
  stability_word,
  j = c(
    "mean_jaccard",
    "median_jaccard",
    "lower_95",
    "upper_95"
  ),
  align = "center",
  part = "all"
)

# Set widths
stability_word <- width(
  stability_word,
  j = "solution",
  width = 1.3
)

stability_word <- width(
  stability_word,
  j = "cluster",
  width = 1.2
)

stability_word <- width(
  stability_word,
  j = c(
    "mean_jaccard",
    "median_jaccard",
    "lower_95",
    "upper_95"
  ),
  width = 1.1
)

# Autofit
stability_word <- autofit(
  stability_word
)

# Create Word document
doc <- read_docx()

# Add table title
doc <- body_add_par(
  doc,
  "Table X. Bootstrap Stability of Cluster Solutions",
  style = "heading 1"
)

# Add table
doc <- body_add_flextable(
  doc,
  value = stability_word
)

# Add explanatory note
doc <- body_add_par(
  doc,
  paste0(
    "Note. Cluster stability was assessed using bootstrap resampling and ",
    "Jaccard similarity. Mean and median Jaccard values summarize the ",
    "stability of cluster membership across bootstrap samples. The 95% ",
    "confidence intervals represent the uncertainty around the estimated ",
    "Jaccard stability."
  )
)

# Export Word document
print(
  doc,
  target = "Cluster_Stability_Comparison.docx"
)

==============================================================================
  # MANUSCRIPT FIGURE: OPERATIONAL PROFILES OF THE RETAINED TWO-CLUSTER SOLUTION
  # ==============================================================================
# This section creates:
#   Figure 5A: continuous cluster profiles (mean +/- SD)
#   Figure 5B: categorical cluster profiles (within-cluster percentages)
#   Figure 5:  combined side-by-side manuscript figure
#
# The plots are descriptive. Do not add hypothesis-test results or confidence
# intervals because all displayed variables were active in the clustering.

if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is required. Install it with install.packages('tidyr').")
}

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop(
    "Package 'patchwork' is required. Install it with ",
    "install.packages('patchwork')."
  )
}

if (!requireNamespace("scales", quietly = TRUE)) {
  stop("Package 'scales' is required. Install it with install.packages('scales').")
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(patchwork)
library(scales)

figure_dir <- "figures"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# Use the same colour-blind-safe palette in both panels.
cluster_colours <- c(
  "Cluster 1" = "#0072B2",
  "Cluster 2" = "#D55E00"
)

# Standardise cluster order once so it remains consistent throughout the figure.
figure_data <- fishers_clustering_2clust %>%
  mutate(
    cluster_2 = factor(
      as.character(cluster_2),
      levels = c("Cluster 1", "Cluster 2")
    )
  )

cluster_sizes_2 <- figure_data %>%
  count(cluster_2, name = "sample_size", .drop = FALSE) %>%
  mutate(
    axis_label = paste0(cluster_2, "\n(n = ", sample_size, ")"),
    legend_label = paste0(cluster_2, " (n = ", sample_size, ")")
  )

cluster_axis_labels <- setNames(
  cluster_sizes_2$axis_label,
  as.character(cluster_sizes_2$cluster_2)
)

cluster_legend_labels <- setNames(
  cluster_sizes_2$legend_label,
  as.character(cluster_sizes_2$cluster_2)
)

# ------------------------------------------------------------------------------
# PANEL A: CONTINUOUS CHARACTERISTICS AS RAINCLOUD PLOTS
# ------------------------------------------------------------------------------
continuous_plot_data <- figure_data %>%
  select(
    cluster_2,
    mean_distance,
    cost_per_trip_usd,
    Fishing.experience,
    vessel_length_m
  ) %>%
  pivot_longer(
    cols = -cluster_2,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    characteristic = recode(
      variable,
      mean_distance = "Average fishing distance (km)",
      cost_per_trip_usd = "Cost per trip ($)",
      Fishing.experience = "Fishing experience (years)",
      vessel_length_m = "Vessel length (m)"
    ),
    characteristic = factor(
      characteristic,
      levels = c(
        "Average fishing distance (km)",
        "Cost per trip ($)",
        "Fishing experience (years)",
        "Vessel length (m)"
      )
    )
  ) %>%
  filter(!is.na(value), !is.na(cluster_2))

# The same mean +/- SD summary is used in all four facets.
continuous_summary <- continuous_plot_data %>%
  group_by(characteristic, cluster_2) %>%
  summarise(
    n = n(),
    cluster_mean = mean(value),
    cluster_sd = sd(value),
    summary_lower = cluster_mean - cluster_sd,
    summary_upper = cluster_mean + cluster_sd,
    .groups = "drop"
  )

plot_continuous_profiles <- ggplot(continuous_plot_data, aes(x = cluster_2, y = value)) +
  # 1. Half-violins aligned to the right
  ggdist::stat_halfeye(
    aes(fill = cluster_2),
    justification = -0.15,
    .width = 0,
    point_colour = NA,
    slab_alpha = 0.42,
    slab_colour = NA,
    scale = 0.68,
    normalize = "panels",  # <--- THIS IS THE FIX
    na.rm = TRUE
  ) +
  # 2. Raw data points jittered to the left
  geom_point(
    aes(colour = cluster_2),
    position = position_jitter(width = 0.08, seed = 2026),
    size = 0.75,
    alpha = 0.48,
    stroke = 0,
    na.rm = TRUE
  ) +
  # 3. Error bars nudged slightly to align cleanly
  geom_errorbar(
    data = continuous_summary,
    aes(
      y = cluster_mean,
      ymin = summary_lower,
      ymax = summary_upper,
      colour = cluster_2
    ),
    width = 0.075,
    linewidth = 0.55,
    position = position_nudge(x = 0.15)
  ) +
  # 4. Mean points placed over the error bars
  geom_point(
    data = continuous_summary,
    aes(
      y = cluster_mean,
      colour = cluster_2
    ),
    shape = 21,
    fill = "white",
    size = 2.4,
    stroke = 0.75,
    position = position_nudge(x = 0.15)
  ) +
  facet_wrap(
    ~ characteristic,
    scales = "free_y",
    ncol = 2
  ) +
  scale_x_discrete(
    labels = cluster_axis_labels,
    expand = expansion(add = c(0.4, 0.4))
  ) +
  scale_y_continuous(
    labels = scales::label_number(big.mark = ",", accuracy = 0.1),
    expand = expansion(mult = c(0.06, 0.10))
  ) +
  scale_fill_manual(
    values = cluster_colours,
    breaks = names(cluster_colours),
    drop = FALSE
  ) +
  scale_colour_manual(
    values = cluster_colours,
    breaks = names(cluster_colours),
    labels = cluster_legend_labels,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Observed value and cluster mean (± SD)",
    fill = "Cluster",
    colour = "Cluster",
    title = "Continuous operational characteristics"
  ) +
  guides(fill = "none", colour = "none") +
  theme_classic(base_size = 10.5) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9.5),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 8.5, lineheight = 0.9),
    axis.title.y = element_text(size = 9.5),
    plot.title = element_text(face = "bold", size = 11),
    panel.spacing = grid::unit(0.85, "lines"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# ------------------------------------------------------------------------------
# PANEL B: CATEGORICAL CHARACTERISTICS (GRID LINES REMOVED)
# ------------------------------------------------------------------------------
plot_categorical_profiles <- ggplot() +
  geom_segment(
    data = categorical_segments,
    aes(
      x = segment_start,
      xend = segment_end,
      y = characteristic,
      yend = characteristic
    ),
    colour = "grey72",
    linewidth = 0.75
  ) +
  geom_point(
    data = categorical_plot_data,
    aes(
      x = percentage,
      y = characteristic,
      colour = cluster_2
    ),
    size = 2.6
  ) +
  scale_colour_manual(
    values = cluster_colours,
    breaks = names(cluster_colours),
    labels = cluster_legend_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Fishers within cluster (%)",
    y = NULL,
    colour = "Cluster",
    title = "Categorical operational characteristics"
  ) +
  guides(
    colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 3.2, alpha = 1)
    )
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.y = element_text(size = 8.6),
    axis.title.x = element_text(size = 9.5),
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.major = element_blank(), # Removed major grid lines
    panel.grid.minor = element_blank(), # Removed minor grid lines
    legend.title = element_blank(),
    legend.text = element_text(size = 9.5),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )
# ------------------------------------------------------------------------------
# COMBINE PANELS
# ------------------------------------------------------------------------------
cluster_profiles_figure <- (
  (plot_continuous_profiles | plot_categorical_profiles) +
    patchwork::plot_layout(
      widths = c(1.05, 1.55),
      guides = "collect"
    ) +
    patchwork::plot_annotation(tag_levels = "A")
) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    plot.tag = element_text(face = "bold", size = 13)
  )

# ------------------------------------------------------------------------------
# PREVIEW THE FIGURES BEFORE SAVING
# ------------------------------------------------------------------------------
# 1. Update Panel A internal theme slightly for better facet spacing
plot_continuous_profiles <- plot_continuous_profiles +
  theme(
    panel.spacing.x = grid::unit(1.2, "lines"),  # Give space between facet columns
    axis.text.x = element_text(size = 8, lineheight = 0.95),
    strip.text = element_text(face = "bold", size = 9)
  )

# 2. Combine with updated patchwork layout widths
cluster_profiles_figure <- (
  (plot_continuous_profiles | plot_categorical_profiles) +
    patchwork::plot_layout(
      widths = c(1.3, 1), # Gives Panel A more width for its 2x2 grid
      guides = "collect"
    ) +
    patchwork::plot_annotation(tag_levels = "A")
) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    plot.tag = element_text(face = "bold", size = 14),
    plot.margin = margin(10, 10, 10, 10)
  )

# ------------------------------------------------------------------------------
# SAVE WITH UPDATED PUBLICATION DIMENSIONS
# ------------------------------------------------------------------------------
# Using a 14 x 7.5 inch canvas ensures zero label overlap across full double-column layouts

png_file <- file.path(output_dir, "Figure_5_cluster_operational_profiles.png")
tiff_file <- file.path(output_dir, "Figure_5_cluster_operational_profiles.tiff")

# Save PNG
ggplot2::ggsave(
  filename = png_file,
  plot = cluster_profiles_figure,
  device = "png",
  width = 14,       # Expanded from 7.5 to 14 inches
  height = 7.5,     # Expanded from 5.5 to 7.5 inches
  units = "in",
  dpi = 600,
  bg = "white"
)

# Save TIFF for journal submission
ggplot2::ggsave(
  filename = tiff_file,
  plot = cluster_profiles_figure,
  device = "tiff",
  width = 14,       # Expanded from 7.5 to 14 inches
  height = 7.5,     # Expanded from 5.5 to 7.5 inches
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

