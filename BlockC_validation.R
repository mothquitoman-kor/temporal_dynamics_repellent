
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

evaporation_file <- "validation_evaporation.csv"
repellency_file  <- "validation_repellency.csv"
blockB_boot_rdata <- "blockB_bootstrap_results.RData"

out_results_csv            <- "validation_results.csv"
out_workspace_rdata        <- "validation_analysis_outputs.RData"
out_evap_completed_csv     <- "validation_evaporation_completed.csv"
out_rep_completed_csv      <- "validation_repellency_completed.csv"

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

rmse <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((obs[ok] - pred[ok])^2))
}

mae <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  mean(abs(obs[ok] - pred[ok]))
}

theme_manuscript <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.35),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 1.0),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black")
    )
}

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
  stop(
    "validation_evaporation.csv is missing required columns: ",
    paste(missing_evap, collapse = ", ")
  )
}
if (!all(c("Validation_ID", "Replicate") %in% names(rep_dat))) {
  stop("validation_repellency.csv must contain Validation_ID and Replicate.")
}

pred_pt50_col <- find_first_existing_col(rep_dat, c("Pred_PT50_min", "Pred_PT50", "Predicted_PT50_min"))
pred_pt50_lower_col <- find_first_existing_col(rep_dat, c("Pred_PT50_lower_min", "Pred_PT50_lower", "PT50_lower"))
pred_pt50_upper_col <- find_first_existing_col(rep_dat, c("Pred_PT50_upper_min", "Pred_PT50_upper", "PT50_upper"))
obs_pt50_col <- find_first_existing_col(rep_dat, c("Obs_PT50", "Obs_PT50_min", "Observed_PT50_min"))

if (is.na(obs_pt50_col)) {
  stop(
    "validation_repellency.csv must contain observed PT50 column named one of: ",
    "Obs_PT50, Obs_PT50_min, Observed_PT50_min."
  )
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
    Obs_PT50_min = as.numeric(.data[[obs_pt50_col]])
  )

need_fill_pred <- is.na(pred_pt50_col) || is.na(pred_pt50_lower_col) || is.na(pred_pt50_upper_col)

if (need_fill_pred) {
  stop_if_missing(blockB_boot_rdata)
  load(blockB_boot_rdata)
  if (!exists("selected_solution_fixed_boot_result")) {
    stop("selected_solution_fixed_boot_result was not found in ", blockB_boot_rdata)
  }
  obj <- selected_solution_fixed_boot_result
  if (is.null(obj$selected_solutions) || nrow(obj$selected_solutions) == 0) {
    stop("selected_solution_fixed_boot_result$selected_solutions is empty.")
  }
  
  sel <- obj$selected_solutions %>%
    mutate(
      Objective = if ("Objective" %in% names(.)) as.character(Objective) else NA_character_,
      Objective_Class = if ("Objective_Class" %in% names(.)) as.character(Objective_Class) else NA_character_,
      Selection_Rank = if ("Selection_Rank" %in% names(.)) suppressWarnings(as.integer(Selection_Rank)) else NA_integer_,
      Validation_ID_lookup = if ("Validation_ID" %in% names(.)) as.character(Validation_ID) else NA_character_
    )
  
  if (all(c("Objective_Class", "Selection_Rank") %in% names(rep_dat)) &&
      all(c("Objective_Class", "Selection_Rank") %in% names(sel))) {
    
    pt_lookup <- sel %>%
      select(Objective_Class, Selection_Rank, PT50, PT50_lower, PT50_upper) %>%
      distinct()
    
    rep_dat <- rep_dat %>%
      mutate(Selection_Rank = suppressWarnings(as.integer(Selection_Rank))) %>%
      left_join(pt_lookup, by = c("Objective_Class", "Selection_Rank"))
    
  } else if (all(c("Objective") %in% names(rep_dat)) &&
             all(c("Objective") %in% names(sel))) {
    
    pt_lookup <- sel %>%
      select(Objective, PT50, PT50_lower, PT50_upper) %>%
      distinct()
    
    rep_dat <- rep_dat %>%
      left_join(pt_lookup, by = "Objective")
    
  } else if (all(c("Validation_ID_lookup") %in% names(sel))) {
    
    pt_lookup <- sel %>%
      filter(!is.na(Validation_ID_lookup) & Validation_ID_lookup != "") %>%
      select(Validation_ID_lookup, PT50, PT50_lower, PT50_upper) %>%
      distinct()
    
    rep_dat <- rep_dat %>%
      left_join(pt_lookup, by = c("Validation_ID" = "Validation_ID_lookup"))
    
  } else {
    stop(
      "Predicted PT50 lower/upper CI columns are missing, and no safe join key was found.\n",
      "Provide Objective_Class + Selection_Rank in validation_repellency.csv, ",
      "or include predicted PT50 lower/upper columns directly."
    )
  }
  
  rep_dat <- rep_dat %>%
    mutate(
      Pred_PT50_min = if (!is.na(pred_pt50_col)) as.numeric(.data[[pred_pt50_col]]) else as.numeric(PT50),
      Pred_PT50_lower_min = if (!is.na(pred_pt50_lower_col)) as.numeric(.data[[pred_pt50_lower_col]]) else as.numeric(PT50_lower),
      Pred_PT50_upper_min = if (!is.na(pred_pt50_upper_col)) as.numeric(.data[[pred_pt50_upper_col]]) else as.numeric(PT50_upper)
    )
} else {
  rep_dat <- rep_dat %>%
    mutate(
      Pred_PT50_min = as.numeric(.data[[pred_pt50_col]]),
      Pred_PT50_lower_min = as.numeric(.data[[pred_pt50_lower_col]]),
      Pred_PT50_upper_min = as.numeric(.data[[pred_pt50_upper_col]])
    )
}

missing_pred_rows <- rep_dat %>%
  filter(!(is.finite(Pred_PT50_min) & is.finite(Pred_PT50_lower_min) & is.finite(Pred_PT50_upper_min)))

if (nrow(missing_pred_rows) > 0) {
  stop(
    "Could not complete predicted PT50 columns for these rows:\n",
    paste(unique(missing_pred_rows$Validation_ID), collapse = ", ")
  )
}

meta_cols <- c("Validation_ID", "Objective_Class", "Selection_Rank", "Selection_Label", "Objective", "Grid_ID")
meta_lookup <- bind_rows(
  evap_dat %>% select(any_of(meta_cols)),
  rep_dat %>% select(any_of(meta_cols))
) %>%
  distinct()

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

composition_dat <- evap_dat %>%
  filter(is.na(Use_for_composition_validation) | Use_for_composition_validation) %>%
  arrange(Validation_ID, Time_min) %>%
  compute_alr_df("Obs_x_citronellal", "Obs_x_citronellol", "Obs_x_geraniol", "Obs") %>%
  compute_alr_df("Pred_x_citronellal_t", "Pred_x_citronellol_t", "Pred_x_geraniol_t", "Pred") %>%
  mutate(
    Comp_Error_ALR = sqrt((Obs_z1 - Pred_z1)^2 + (Obs_z2 - Pred_z2)^2)
  )

composition_metrics_by_solution <- composition_dat %>%
  group_by(Validation_ID) %>%
  summarise(
    n_composition = sum(is.finite(Comp_Error_ALR)),
    Composition_RMSE_ALR = sqrt(mean(Comp_Error_ALR[is.finite(Comp_Error_ALR)]^2)),
    Composition_MAE_ALR  = mean(Comp_Error_ALR[is.finite(Comp_Error_ALR)]),
    .groups = "drop"
  )

pt50_rep_level <- rep_dat %>%
  mutate(
    Obs_in_Pred95CI = ifelse(
      is.finite(Obs_PT50_min) &
        is.finite(Pred_PT50_lower_min) &
        is.finite(Pred_PT50_upper_min),
      Obs_PT50_min >= Pred_PT50_lower_min & Obs_PT50_min <= Pred_PT50_upper_min,
      NA
    )
  )

pt50_by_solution <- pt50_rep_level %>%
  group_by(Validation_ID) %>%
  summarise(
    n_pt50 = sum(is.finite(Obs_PT50_min) & is.finite(Pred_PT50_min)),
    n_pt50_checked = sum(!is.na(Obs_in_Pred95CI)),
    n_pt50_within_pred95CI = sum(Obs_in_Pred95CI, na.rm = TRUE),
    Prop_pt50_within_pred95CI = ifelse(
      n_pt50_checked > 0,
      n_pt50_within_pred95CI / n_pt50_checked,
      NA_real_
    ),
    Pred_PT50_min = dplyr::first(Pred_PT50_min[is.finite(Pred_PT50_min)]),
    Pred_PT50_lower_min = dplyr::first(Pred_PT50_lower_min[is.finite(Pred_PT50_lower_min)]),
    Pred_PT50_upper_min = dplyr::first(Pred_PT50_upper_min[is.finite(Pred_PT50_upper_min)]),
    Obs_PT50_mean_min = mean(Obs_PT50_min, na.rm = TRUE),
    Obs_PT50_sd_min = ifelse(
      sum(is.finite(Obs_PT50_min)) > 1,
      sd(Obs_PT50_min, na.rm = TRUE),
      NA_real_
    ),
    Obs_PT50_se_min = ifelse(
      sum(is.finite(Obs_PT50_min)) > 1,
      sd(Obs_PT50_min, na.rm = TRUE) / sqrt(sum(is.finite(Obs_PT50_min))),
      NA_real_
    ),
    PT50_RMSE_min = rmse(Obs_PT50_min, Pred_PT50_min),
    PT50_MAE_min  = mae(Obs_PT50_min, Pred_PT50_min),
    PT50_Bias_min = mean(Obs_PT50_min - Pred_PT50_min, na.rm = TRUE),
    PT50_All_Replicates_Within_Pred95CI = ifelse(
      n_pt50_checked > 0,
      all(Obs_in_Pred95CI[!is.na(Obs_in_Pred95CI)]),
      NA
    ),
    .groups = "drop"
  )

validation_results <- meta_lookup %>%
  distinct(Validation_ID, .keep_all = TRUE) %>%
  left_join(weight_metrics_by_solution, by = "Validation_ID") %>%
  left_join(composition_metrics_by_solution, by = "Validation_ID") %>%
  left_join(pt50_by_solution, by = "Validation_ID") %>%
  mutate(
    Plot_Label = make_plot_label(Validation_ID, Objective_Class, Selection_Rank, Selection_Label),
    PT50_Consistent = PT50_All_Replicates_Within_Pred95CI,
    Model_Consistency_Note = dplyr::case_when(
      isTRUE(PT50_Consistent) ~ "All observed PT50 replicates fell within the predicted 95% CI.",
      identical(PT50_Consistent, FALSE) ~ "At least one observed PT50 replicate fell outside the predicted 95% CI.",
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
    n_pt50_checked,
    n_pt50_within_pred95CI,
    Prop_pt50_within_pred95CI,
    Pred_PT50_min,
    Pred_PT50_lower_min,
    Pred_PT50_upper_min,
    Obs_PT50_mean_min,
    Obs_PT50_sd_min,
    Obs_PT50_se_min,
    PT50_RMSE_min,
    PT50_MAE_min,
    PT50_Bias_min,
    PT50_Consistent,
    Model_Consistency_Note
  ) %>%
  arrange(Validation_ID)

weight_plot_dat <- weight_dat %>%
  distinct(
    Validation_ID, Time_min, Pred_remaining_dose_mg, Obs_Total_Mass_mg,
    Objective_Class, Selection_Rank, Selection_Label, Objective
  ) %>%
  mutate(Plot_Label = make_plot_label(Validation_ID, Objective_Class, Selection_Rank, Selection_Label))

p_weight <- ggplot(weight_plot_dat, aes(x = Time_min, y = Pred_remaining_dose_mg)) +
  geom_line(linewidth = 1) +
  geom_point(aes(y = Obs_Total_Mass_mg), size = 2) +
  facet_wrap(~ Plot_Label, scales = "free_y") +
  theme_manuscript() +
  labs(
    x = "Time (min)",
    y = "Remaining dose (mg)",
    title = "Validation of remaining dose trajectories"
  )

composition_plot_dat <- composition_dat %>%
  select(
    Validation_ID, Time_min,
    Pred_x_citronellal_t, Pred_x_citronellol_t, Pred_x_geraniol_t,
    Obs_x_citronellal, Obs_x_citronellol, Obs_x_geraniol,
    Objective_Class, Selection_Rank, Selection_Label
  ) %>%
  pivot_longer(
    cols = c(
      Pred_x_citronellal_t, Pred_x_citronellol_t, Pred_x_geraniol_t,
      Obs_x_citronellal, Obs_x_citronellol, Obs_x_geraniol
    ),
    names_to = c("Type", "Component"),
    names_pattern = "(Pred|Obs)_(.*)"
  ) %>%
  mutate(
    Plot_Label = make_plot_label(Validation_ID, Objective_Class, Selection_Rank, Selection_Label),
    Component = recode(
      Component,
      citronellal_t = "Citronellal",
      citronellol_t = "α-Citronellol",
      geraniol_t = "Geraniol",
      citronellal = "Citronellal",
      citronellol = "α-Citronellol",
      geraniol = "Geraniol"
    )
  )

p_comp <- ggplot(composition_plot_dat, aes(x = Time_min, y = value, colour = Component)) +
  geom_line(data = subset(composition_plot_dat, Type == "Pred"), linewidth = 1) +
  geom_point(data = subset(composition_plot_dat, Type == "Obs"), size = 2) +
  facet_wrap(~ Plot_Label, scales = "free_y") +
  theme_manuscript() +
  labs(
    x = "Time (min)",
    y = "Composition fraction",
    title = "Validation of composition trajectories"
  )

pt50_plot_dat <- validation_results
rep_plot_dat <- pt50_rep_level %>%
  left_join(pt50_plot_dat %>% select(Validation_ID, Plot_Label), by = "Validation_ID")

p_pt50_solution <- ggplot(rep_plot_dat, aes(x = Plot_Label, y = Obs_PT50_min)) +
  geom_point(
    aes(shape = Obs_in_Pred95CI),
    position = position_jitter(width = 0.08, height = 0),
    size = 2,
    alpha = 0.8
  ) +
  geom_errorbar(
    data = pt50_plot_dat,
    aes(x = Plot_Label, ymin = Pred_PT50_lower_min, ymax = Pred_PT50_upper_min),
    width = 0.12,
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
  geom_point(
    data = pt50_plot_dat,
    aes(x = Plot_Label, y = Obs_PT50_mean_min),
    size = 3,
    inherit.aes = FALSE
  ) +
  theme_manuscript() +
  labs(
    x = NULL,
    y = "PT50 (min)",
    title = "Observed replicate PT50 vs predicted PT50 interval",
    shape = "Inside\npredicted 95% CI"
  )

p_pt50_scatter <- ggplot(pt50_plot_dat, aes(x = Pred_PT50_min, y = Obs_PT50_mean_min)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_segment(
    aes(x = Pred_PT50_lower_min, xend = Pred_PT50_upper_min,
        y = Obs_PT50_mean_min, yend = Obs_PT50_mean_min),
    linewidth = 0.7
  ) +
  geom_point(size = 3) +
  theme_manuscript() +
  labs(
    x = "Predicted PT50 (min)",
    y = "Observed mean PT50 (min)",
    title = "PT50 validation"
  )

write.csv(evap_dat, out_evap_completed_csv, row.names = FALSE)
write.csv(rep_dat, out_rep_completed_csv, row.names = FALSE)
write.csv(validation_results, out_results_csv, row.names = FALSE)

save(
  evap_dat, rep_dat,
  weight_dat, composition_dat,
  pt50_rep_level,
  weight_metrics_by_solution,
  composition_metrics_by_solution,
  pt50_by_solution,
  validation_results,
  p_weight, p_comp, p_pt50_solution, p_pt50_scatter,
  file = out_workspace_rdata
)

cat("Validation analysis finished.\n")
cat("Saved:\n")
cat(" - ", out_evap_completed_csv, "\n", sep = "")
cat(" - ", out_rep_completed_csv, "\n", sep = "")
cat(" - ", out_results_csv, "\n", sep = "")
cat(" - ", out_workspace_rdata, "\n", sep = "")
cat("Rows in validation_results.csv: ", nrow(validation_results), "\n", sep = "")
print(p_weight)
print(p_comp)
print(p_pt50_solution)
print(p_pt50_scatter)
