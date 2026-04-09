# =====================================================================
# BLOCK C
# VALIDATION
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("mixture_pipeline_functions.R")

# INPUT FILES
evaporation_file <- "validation_evaporation.csv"
repellency_file  <- "validation_repellency.csv"

# OUTPUT FILES
out_results_csv     <- "validation_results.csv"
out_workspace_rdata <- "validation_analysis_outputs.RData"

# HELPERS
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

parse_logical_flag <- function(x) {
  if (is.logical(x)) return(x)
  x_chr <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x_chr %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

find_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

compute_alr_df <- function(df, x1_col, x2_col, x3_col, prefix) {
  alr_mat <- t(mapply(
    FUN = function(a, b, c) alr_from_simplex(a, b, c),
    df[[x1_col]],
    df[[x2_col]],
    df[[x3_col]]
  ))
  df[[paste0(prefix, "_z1")]] <- as.numeric(alr_mat[, "z1"])
  df[[paste0(prefix, "_z2")]] <- as.numeric(alr_mat[, "z2"])
  df
}

safe_t_ci_mean <- function(x, level = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  m <- mean(x)
  if (n == 1) return(c(lower = NA_real_, upper = NA_real_))
  se <- stats::sd(x) / sqrt(n)
  crit <- stats::qt((1 + level) / 2, df = n - 1)
  c(lower = m - crit * se, upper = m + crit * se)
}

make_plot_label <- function(validation_id, objective_class = NA, selection_rank = NA, selection_label = NA) {
  validation_id <- as.character(validation_id)
  objective_class <- as.character(objective_class)
  selection_rank <- as.character(selection_rank)
  selection_label <- as.character(selection_label)

  has_meta <- !(is.na(objective_class) | objective_class == "") |
              !(is.na(selection_label) | selection_label == "")

  out <- ifelse(
    has_meta,
    paste0(
      validation_id, "\n",
      ifelse(is.na(objective_class) | objective_class == "", "", objective_class),
      ifelse(is.na(selection_rank) | selection_rank == "", "", paste0(" #", selection_rank)),
      ifelse(is.na(selection_label) | selection_label == "", "", paste0(" ", selection_label))
    ),
    validation_id
  )
  trimws(out)
}

# READ + CLEAN
stop_if_missing(evaporation_file)
stop_if_missing(repellency_file)

evap_dat <- read.csv(evaporation_file, stringsAsFactors = FALSE)
rep_dat  <- read.csv(repellency_file, stringsAsFactors = FALSE)

required_evap <- c(
  "Validation_ID", "Time_min",
  "Pred_remaining_dose_mg", "Obs_Total_Mass_mg",
  "Pred_x_citronellal_t", "Pred_x_citronellol_t", "Pred_x_geraniol_t",
  "Obs_x_citronellal", "Obs_x_citronellol", "Obs_x_geraniol"
)
missing_evap <- setdiff(required_evap, names(evap_dat))
if (length(missing_evap) > 0) {
  stop("validation_evaporation.csv is missing required columns: ",
       paste(missing_evap, collapse = ", "))
}

pred_pt50_col <- find_first_existing_col(rep_dat, c("Pred_PT50_min", "Pred_PT50", "Predicted_PT50_min"))
obs_pt50_col  <- find_first_existing_col(rep_dat, c("Obs_PT50", "Obs_PT50_min", "Observed_PT50_min"))
if (is.na(pred_pt50_col) || is.na(obs_pt50_col)) {
  stop("validation_repellency.csv must contain predicted and observed PT50 columns.")
}

required_rep <- c("Validation_ID", "Replicate", pred_pt50_col, obs_pt50_col)
missing_rep <- setdiff(required_rep, names(rep_dat))
if (length(missing_rep) > 0) {
  stop("validation_repellency.csv is missing required columns: ",
       paste(missing_rep, collapse = ", "))
}

evap_dat <- evap_dat %>%
  mutate(
    Validation_ID = as.character(Validation_ID),
    Time_min = as.numeric(Time_min),
    Use_for_weight_validation = parse_logical_flag(Use_for_weight_validation),
    Use_for_composition_validation = parse_logical_flag(Use_for_composition_validation),
    Pred_remaining_dose_mg = as.numeric(Pred_remaining_dose_mg),
    Obs_Total_Mass_mg = as.numeric(Obs_Total_Mass_mg),
    Pred_x_citronellal_t = as.numeric(Pred_x_citronellal_t),
    Pred_x_citronellol_t = as.numeric(Pred_x_citronellol_t),
    Pred_x_geraniol_t = as.numeric(Pred_x_geraniol_t),
    Obs_x_citronellal = as.numeric(Obs_x_citronellal),
    Obs_x_citronellol = as.numeric(Obs_x_citronellol),
    Obs_x_geraniol = as.numeric(Obs_x_geraniol)
  )

rep_dat <- rep_dat %>%
  mutate(
    Validation_ID = as.character(Validation_ID),
    Replicate = as.character(Replicate),
    Pred_PT50_min = as.numeric(.data[[pred_pt50_col]]),
    Obs_PT50_min = as.numeric(.data[[obs_pt50_col]])
  )

meta_cols <- c("Validation_ID", "Objective_Class", "Selection_Rank", "Selection_Label", "Objective", "Grid_ID")
meta_lookup <- bind_rows(
  evap_dat %>% select(any_of(meta_cols)),
  rep_dat %>% select(any_of(meta_cols))
) %>%
  distinct()

# 1) WEIGHT VALIDATION
weight_dat <- evap_dat %>%
  filter(is.na(Use_for_weight_validation) | Use_for_weight_validation) %>%
  arrange(Validation_ID, Time_min)

weight_metrics_by_solution <- weight_dat %>%
  group_by(Validation_ID) %>%
  summarise(
    n_weight = sum(is.finite(Obs_Total_Mass_mg) & is.finite(Pred_remaining_dose_mg)),
    Weight_RMSE_mg = rmse(Obs_Total_Mass_mg, Pred_remaining_dose_mg),
    Weight_MAE_mg  = mae(Obs_Total_Mass_mg, Pred_remaining_dose_mg),
    .groups = "drop"
  )

# 2) COMPOSITION VALIDATION (ALR SPACE)
composition_dat <- evap_dat %>%
  filter(!is.na(Use_for_composition_validation) & Use_for_composition_validation) %>%
  arrange(Validation_ID, Time_min)

if (nrow(composition_dat) == 0) {
  stop("No rows were flagged for composition validation.")
}

if (all(c("Obs_z1", "Obs_z2") %in% names(composition_dat))) {
  composition_dat <- composition_dat %>%
    mutate(
      Obs_z1 = as.numeric(Obs_z1),
      Obs_z2 = as.numeric(Obs_z2)
    )
} else {
  composition_dat <- compute_alr_df(
    composition_dat,
    x1_col = "Obs_x_citronellal",
    x2_col = "Obs_x_citronellol",
    x3_col = "Obs_x_geraniol",
    prefix = "Obs"
  )
}

composition_dat <- compute_alr_df(
  composition_dat,
  x1_col = "Pred_x_citronellal_t",
  x2_col = "Pred_x_citronellol_t",
  x3_col = "Pred_x_geraniol_t",
  prefix = "Pred"
) %>%
  mutate(
    Composition_Error_ALR = sqrt((Obs_z1 - Pred_z1)^2 + (Obs_z2 - Pred_z2)^2)
  )

composition_metrics_by_solution <- composition_dat %>%
  group_by(Validation_ID) %>%
  summarise(
    n_composition = sum(is.finite(Composition_Error_ALR)),
    Composition_RMSE_ALR = sqrt(mean(Composition_Error_ALR^2, na.rm = TRUE)),
    Composition_MAE_ALR  = mean(Composition_Error_ALR, na.rm = TRUE),
    .groups = "drop"
  )

# 3) PT50 VALIDATION
pt50_by_solution <- rep_dat %>%
  group_by(Validation_ID) %>%
  summarise(
    n_pt50 = sum(is.finite(Obs_PT50_min) & is.finite(Pred_PT50_min)),
    Pred_PT50_min = dplyr::first(Pred_PT50_min[is.finite(Pred_PT50_min)]),
    Obs_PT50_mean_min = mean(Obs_PT50_min, na.rm = TRUE),
    Obs_PT50_sd_min = ifelse(sum(is.finite(Obs_PT50_min)) > 1, sd(Obs_PT50_min, na.rm = TRUE), NA_real_),
    Obs_PT50_se_min = ifelse(sum(is.finite(Obs_PT50_min)) > 1, sd(Obs_PT50_min, na.rm = TRUE) / sqrt(sum(is.finite(Obs_PT50_min))), NA_real_),
    PT50_RMSE_min = rmse(Obs_PT50_min, Pred_PT50_min),
    PT50_MAE_min = mae(Obs_PT50_min, Pred_PT50_min),
    PT50_Bias_min = mean(Obs_PT50_min - Pred_PT50_min, na.rm = TRUE),
    PT50_Mean_Error_min = mean(Obs_PT50_min, na.rm = TRUE) - dplyr::first(Pred_PT50_min[is.finite(Pred_PT50_min)]),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    Obs_PT50_CI_lower_min = safe_t_ci_mean(rep_dat$Obs_PT50_min[rep_dat$Validation_ID == Validation_ID])["lower"],
    Obs_PT50_CI_upper_min = safe_t_ci_mean(rep_dat$Obs_PT50_min[rep_dat$Validation_ID == Validation_ID])["upper"],
    Pred_in_Obs95CI_mean = ifelse(
      is.finite(Pred_PT50_min) & is.finite(Obs_PT50_CI_lower_min) & is.finite(Obs_PT50_CI_upper_min),
      Pred_PT50_min >= Obs_PT50_CI_lower_min & Pred_PT50_min <= Obs_PT50_CI_upper_min,
      NA
    )
  ) %>%
  ungroup()


# 4) SINGLE SIX-ROW SUMMARY TABLE
validation_results <- meta_lookup %>%
  distinct(Validation_ID, .keep_all = TRUE) %>%
  left_join(weight_metrics_by_solution, by = "Validation_ID") %>%
  left_join(composition_metrics_by_solution, by = "Validation_ID") %>%
  left_join(pt50_by_solution, by = "Validation_ID") %>%
  mutate(
    Plot_Label = make_plot_label(Validation_ID, Objective_Class, Selection_Rank, Selection_Label),
    PT50_Consistent = Pred_in_Obs95CI_mean,
    Model_Consistency_Note = dplyr::case_when(
      isTRUE(PT50_Consistent) ~ "Predicted PT50 is inside observed 95% CI of mean.",
      identical(PT50_Consistent, FALSE) ~ "Predicted PT50 is outside observed 95% CI of mean.",
      TRUE ~ "PT50 consistency could not be determined."
    )
  ) %>%
  select(
    Validation_ID,
    Objective_Class,
    Selection_Rank,
    Selection_Label,
    Objective,
    Grid_ID,
    Plot_Label,
    n_weight,
    Weight_RMSE_mg,
    Weight_MAE_mg,
    n_composition,
    Composition_RMSE_ALR,
    Composition_MAE_ALR,
    n_pt50,
    Pred_PT50_min,
    Obs_PT50_mean_min,
    Obs_PT50_sd_min,
    Obs_PT50_se_min,
    Obs_PT50_CI_lower_min,
    Obs_PT50_CI_upper_min,
    PT50_RMSE_min,
    PT50_MAE_min,
    PT50_Bias_min,
    PT50_Mean_Error_min,
    PT50_Consistent,
    Model_Consistency_Note
  ) %>%
  arrange(Objective_Class, suppressWarnings(as.numeric(Selection_Rank)), Validation_ID)

# 5) PLOTS
weight_plot_dat <- weight_dat %>%
  distinct(Validation_ID, Time_min, Pred_remaining_dose_mg, Obs_Total_Mass_mg,
           Objective_Class, Selection_Rank, Selection_Label, Objective) %>%
  mutate(Plot_Label = make_plot_label(Validation_ID, Objective_Class, Selection_Rank, Selection_Label))

p_weight <- ggplot(weight_plot_dat, aes(x = Time_min)) +
  geom_line(aes(y = Pred_remaining_dose_mg), linewidth = 0.9) +
  geom_point(aes(y = Obs_Total_Mass_mg), size = 2) +
  geom_line(aes(y = Obs_Total_Mass_mg), linewidth = 0.5, alpha = 0.7) +
  facet_wrap(~ Plot_Label, scales = "free_y") +
  theme_manuscript() +
  labs(
    x = "Time (min)",
    y = "Remaining mass (mg)",
    title = "Weight trajectory validation"
  )

composition_plot_dat <- composition_dat %>%
  distinct(Validation_ID, Time_min,
           Pred_x_citronellal_t, Pred_x_citronellol_t, Pred_x_geraniol_t,
           Obs_x_citronellal, Obs_x_citronellol, Obs_x_geraniol,
           Objective_Class, Selection_Rank, Selection_Label, Objective) %>%
  mutate(Plot_Label = make_plot_label(Validation_ID, Objective_Class, Selection_Rank, Selection_Label))

composition_long <- bind_rows(
  composition_plot_dat %>%
    transmute(Validation_ID, Plot_Label, Time_min, Component = "Citronellal", Type = "Predicted", Value = Pred_x_citronellal_t),
  composition_plot_dat %>%
    transmute(Validation_ID, Plot_Label, Time_min, Component = "Citronellal", Type = "Observed", Value = Obs_x_citronellal),
  composition_plot_dat %>%
    transmute(Validation_ID, Plot_Label, Time_min, Component = "Citronellol", Type = "Predicted", Value = Pred_x_citronellol_t),
  composition_plot_dat %>%
    transmute(Validation_ID, Plot_Label, Time_min, Component = "Citronellol", Type = "Observed", Value = Obs_x_citronellol),
  composition_plot_dat %>%
    transmute(Validation_ID, Plot_Label, Time_min, Component = "Geraniol", Type = "Predicted", Value = Pred_x_geraniol_t),
  composition_plot_dat %>%
    transmute(Validation_ID, Plot_Label, Time_min, Component = "Geraniol", Type = "Observed", Value = Obs_x_geraniol)
)

p_comp <- ggplot(composition_long, aes(x = Time_min, y = Value)) +
  geom_line(data = subset(composition_long, Type == "Predicted"), linewidth = 0.9) +
  geom_point(data = subset(composition_long, Type == "Observed"), size = 2) +
  geom_line(data = subset(composition_long, Type == "Observed"), linewidth = 0.5, alpha = 0.7) +
  facet_grid(Component ~ Plot_Label, scales = "free_y") +
  theme_manuscript() +
  labs(
    x = "Time (min)",
    y = "Composition fraction",
    title = "Composition trajectory validation"
  )

pt50_plot_dat <- validation_results
rep_plot_dat <- rep_dat %>%
  left_join(pt50_plot_dat %>% select(Validation_ID, Plot_Label), by = "Validation_ID")

p_pt50_solution <- ggplot(rep_plot_dat, aes(x = Plot_Label, y = Obs_PT50_min)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2, alpha = 0.75) +
  geom_errorbar(
    data = pt50_plot_dat,
    aes(x = Plot_Label, ymin = Obs_PT50_CI_lower_min, ymax = Obs_PT50_CI_upper_min),
    width = 0.12,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = pt50_plot_dat,
    aes(x = Plot_Label, y = Obs_PT50_mean_min),
    size = 3,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = pt50_plot_dat,
    aes(x = Plot_Label, y = Pred_PT50_min),
    shape = 4,
    size = 3.5,
    stroke = 1.1,
    inherit.aes = FALSE
  ) +
  theme_manuscript() +
  labs(
    x = NULL,
    y = "PT50 (min)",
    title = "Observed replicate PT50 vs predicted PT50"
  )

p_pt50_scatter <- ggplot(pt50_plot_dat, aes(x = Pred_PT50_min, y = Obs_PT50_mean_min)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = Obs_PT50_CI_lower_min, ymax = Obs_PT50_CI_upper_min), width = 0) +
  geom_point(size = 3) +
  theme_manuscript() +
  labs(
    x = "Predicted PT50 (min)",
    y = "Observed mean PT50 (min)",
    title = "PT50 validation"
  )

print(p_weight)
print(p_comp)
print(p_pt50_solution)
print(p_pt50_scatter)

# 6) SAVE
write.csv(validation_results, out_results_csv, row.names = FALSE)

save(
  evap_dat, rep_dat,
  weight_dat, composition_dat,
  weight_metrics_by_solution,
  composition_metrics_by_solution,
  pt50_by_solution,
  validation_results,
  p_weight, p_comp, p_pt50_solution, p_pt50_scatter,
  file = out_workspace_rdata
)

cat("Validation analysis finished.\n")
cat("Saved:\n")
cat(" - ", out_results_csv, "\n", sep = "")
cat(" - ", out_workspace_rdata, "\n", sep = "")
cat("Rows in validation_results.csv: ", nrow(validation_results), "\n", sep = "")
