# =====================================================================
# BLOCK A
# POINT ESTIMATION + MODEL SELECTION + DETERMINISTIC OPTIMIZATION
# =====================================================================

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtern)

source("mixture_pipeline_functions.R")


# USER SETTINGS

weight_csv <- "1_evaporation_weight_template.csv"
gcms_csv <- "2_evaporation_composition_template.csv"
ed50_csv <- "3_inherent_repellency_ed50_template.csv"

blockA_model_rdata <- "blockA_point_model_results.RData"
blockA_opt_rdata <- "blockA_point_optimization_results.RData"

my_lowers <- c(0.1494, 0.0816, 0.1407)
my_uppers <- c(0.6270, 0.5275, 0.5377)

my_initial_dose_mg <- 50
my_time_seq <- seq(0, 720, by = 1)

opt_grid_res <- 80
opt_time_seq_search <- seq(0, 720, by = 1)

evap_trace <- 0
verbose_A <- TRUE

selection_threshold_frac <- 0.90


# LOCAL HELPERS
add_alr_coordinates <- function(df) {
  alr_mat <- t(vapply(
    seq_len(nrow(df)),
    FUN = function(i) {
      alr_from_simplex(
        df$x_citronellal[i],
        df$x_citronellol[i],
        df$x_geraniol[i]
      )
    },
    FUN.VALUE = c(z1 = 0, z2 = 0)
  ))

  dplyr::bind_cols(
    df,
    data.frame(
      z1 = alr_mat[, "z1"],
      z2 = alr_mat[, "z2"],
      stringsAsFactors = FALSE
    )
  )
}

min_alr_distance_to_selected <- function(candidate_df, selected_df) {
  cand_mat <- as.matrix(candidate_df[, c("z1", "z2")])
  sel_mat <- as.matrix(selected_df[, c("z1", "z2")])

  vapply(
    seq_len(nrow(cand_mat)),
    FUN = function(i) {
      dvec <- sqrt((cand_mat[i, 1] - sel_mat[, 1])^2 + (cand_mat[i, 2] - sel_mat[, 2])^2)
      min(dvec)
    },
    FUN.VALUE = numeric(1)
  )
}

select_three_diverse_solutions <- function(surface_df,
                                           value_col,
                                           objective_class,
                                           threshold_frac = 0.90,
                                           label_prefix = objective_class,
                                           verbose = TRUE) {
  cand_all <- surface_df %>%
    dplyr::filter(is.finite(.data[[value_col]]))

  if (nrow(cand_all) < 3) {
    stop("Fewer than 3 finite candidates found for ", objective_class, ".")
  }

  value_max <- max(cand_all[[value_col]], na.rm = TRUE)
  value_cut <- threshold_frac * value_max

  cand_keep <- cand_all %>%
    dplyr::filter(.data[[value_col]] >= value_cut)

  if (nrow(cand_keep) < 3) {
    warning(
      objective_class,
      ": fewer than 3 candidates remained after the ",
      threshold_frac * 100,
      "% filter. Falling back to all finite candidates."
    )
    cand_keep <- cand_all
  }

  cand_keep <- add_alr_coordinates(cand_keep)

  best_row <- cand_keep %>%
    dplyr::slice_max(order_by = .data[[value_col]], n = 1, with_ties = FALSE)

  remaining_1 <- cand_keep %>%
    dplyr::filter(Grid_ID != best_row$Grid_ID[1])

  if (nrow(remaining_1) < 2) {
    stop("Not enough remaining candidates for ", objective_class, " after selecting the best row.")
  }

  remaining_1$min_dist <- min_alr_distance_to_selected(remaining_1, best_row)
  second_row <- remaining_1 %>%
    dplyr::slice_max(order_by = min_dist, n = 1, with_ties = FALSE)

  selected_12 <- dplyr::bind_rows(best_row, second_row)

  remaining_2 <- cand_keep %>%
    dplyr::filter(!(Grid_ID %in% selected_12$Grid_ID))

  if (nrow(remaining_2) < 1) {
    stop("Not enough remaining candidates for ", objective_class, " to select the third row.")
  }

  remaining_2$min_dist <- min_alr_distance_to_selected(remaining_2, selected_12)
  third_row <- remaining_2 %>%
    dplyr::slice_max(order_by = min_dist, n = 1, with_ties = FALSE)

  selected_out <- dplyr::bind_rows(best_row, second_row, third_row) %>%
    dplyr::mutate(
      Objective_Class = objective_class,
      Selection_Rank = c(1L, 2L, 3L),
      Selection_Label = c("Best", "Alt1", "Alt2"),
      Objective = paste0(label_prefix, "_", Selection_Label)
    ) %>%
    dplyr::select(
      Objective,
      Objective_Class,
      Selection_Rank,
      Selection_Label,
      Grid_ID,
      x_citronellal,
      x_citronellol,
      x_geraniol,
      Initial_Dose_mg,
      dplyr::all_of(value_col),
      z1,
      z2
    )

  if (verbose) {
    cat("\n--- ", objective_class, " candidate filter ---\n", sep = "")
    cat("Maximum ", value_col, ": ", signif(value_max, 6), "\n", sep = "")
    cat("Threshold (>= ", threshold_frac * 100, "% of max): ", signif(value_cut, 6), "\n", sep = "")
    cat("Candidates kept: ", nrow(cand_keep), " of ", nrow(cand_all), "\n", sep = "")
    cat("Selected Grid_IDs: ", paste(selected_out$Grid_ID, collapse = ", "), "\n", sep = "")
  }

  selected_out
}

build_six_optima <- function(lower_bounds,
                             upper_bounds,
                             evap_model,
                             ed50_fit,
                             initial_dose_mg,
                             grid_res = 80,
                             time_seq_search = seq(0, 720, by = 10),
                             time_seq_final = seq(0, 720, by = 1),
                             threshold_frac = 0.90,
                             verbose = TRUE) {
  feasible_grid <- generate_feasible_mixture_grid(
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    res = grid_res
  )

  if (verbose) {
    cat("\nNumber of feasible candidate blends:", nrow(feasible_grid), "\n")
    cat("Initial dose fixed at:", initial_dose_mg, "mg\n")
  }

  surface_df <- evaluate_mixture_grid(
    feasible_grid = feasible_grid,
    evap_model = evap_model,
    ed50_fit = ed50_fit,
    initial_dose_mg = initial_dose_mg,
    time_seq = time_seq_search,
    verbose = verbose
  )

  pt50_selected <- select_three_diverse_solutions(
    surface_df = surface_df,
    value_col = "PT50",
    objective_class = "PT50",
    threshold_frac = threshold_frac,
    label_prefix = "PT50",
    verbose = verbose
  )

  ru_selected <- select_three_diverse_solutions(
    surface_df = surface_df,
    value_col = "RU_t0",
    objective_class = "RU0",
    threshold_frac = threshold_frac,
    label_prefix = "RU0",
    verbose = verbose
  )

  selected_df <- dplyr::bind_rows(pt50_selected, ru_selected) %>%
    dplyr::select(
      Objective,
      Objective_Class,
      Selection_Rank,
      Selection_Label,
      Grid_ID,
      x_citronellal,
      x_citronellol,
      x_geraniol,
      Initial_Dose_mg
    )

  curve_list <- vector("list", nrow(selected_df))
  summary_list <- vector("list", nrow(selected_df))

  for (i in seq_len(nrow(selected_df))) {
    init_comp <- c(
      selected_df$x_citronellal[i],
      selected_df$x_citronellol[i],
      selected_df$x_geraniol[i]
    )
    names(init_comp) <- c("x_citronellal", "x_citronellol", "x_geraniol")

    eval_i <- evaluate_one_solution_from_composition(
      init_comp = init_comp,
      evap_model = evap_model,
      ed50_fit = ed50_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq_final
    )

    curve_list[[i]] <- eval_i$curve %>%
      dplyr::mutate(
        Objective = selected_df$Objective[i],
        Objective_Class = selected_df$Objective_Class[i],
        Selection_Rank = selected_df$Selection_Rank[i],
        Selection_Label = selected_df$Selection_Label[i]
      )

    summary_list[[i]] <- data.frame(
      Objective = selected_df$Objective[i],
      Objective_Class = selected_df$Objective_Class[i],
      Selection_Rank = selected_df$Selection_Rank[i],
      Selection_Label = selected_df$Selection_Label[i],
      Grid_ID = selected_df$Grid_ID[i],
      x_citronellal = selected_df$x_citronellal[i],
      x_citronellol = selected_df$x_citronellol[i],
      x_geraniol = selected_df$x_geraniol[i],
      Initial_Dose_mg = initial_dose_mg,
      ED50_t0 = eval_i$ED50_t0,
      log_ED50_t0 = eval_i$log_ED50_t0,
      RU_t0 = eval_i$RU_t0,
      PT50 = eval_i$PT50,
      log_PT50 = eval_i$log_PT50,
      stringsAsFactors = FALSE
    )
  }

  point_summary <- dplyr::bind_rows(summary_list)
  point_curves <- dplyr::bind_rows(curve_list)

  if (verbose) {
    cat("\n--- Six selected solutions (point estimates) ---\n")
    print_full_df(point_summary)
  }

  list(
    feasible_grid = feasible_grid,
    surface_df = surface_df,
    selected_solutions = point_summary,
    point_curves = point_curves,
    selection_threshold_frac = threshold_frac
  )
}


# 1. FIT ED50 MODEL
my_ed50_fit <- fit_ed50_surface(
  ed50_csv = ed50_csv,
  verbose = verbose_A
)

plot_ed50_surface(
  model = my_ed50_fit$final_model,
  lower_bounds = my_lowers,
  upper_bounds = my_uppers
)

plot_log_ed50_surface(
  model = my_ed50_fit$final_model,
  lower_bounds = my_lowers,
  upper_bounds = my_uppers
)

# 2. SELECT EVAPORATION MODEL
model_selection_result <- select_evaporation_model(
  weight_csv = weight_csv,
  gcms_csv = gcms_csv,
  trace = evap_trace,
  verbose = verbose_A
)

final_evaporation_model <- model_selection_result$final_model

# 3. POINT-ESTIMATE OPTIMIZATION
three_opt_result <- build_six_optima(
  lower_bounds = my_lowers,
  upper_bounds = my_uppers,
  evap_model = final_evaporation_model,
  ed50_fit = my_ed50_fit,
  initial_dose_mg = my_initial_dose_mg,
  grid_res = opt_grid_res,
  time_seq_search = opt_time_seq_search,
  time_seq_final = my_time_seq,
  threshold_frac = selection_threshold_frac,
  verbose = TRUE
)

cat("\n--- Block A selected solutions ---\n")
print_full_df(three_opt_result$selected_solutions)

selected_for_ru_surface <- three_opt_result$selected_solutions %>%
  dplyr::filter(Objective_Class == "RU0") %>%
  dplyr::mutate(Objective = "Max RU_t0")

selected_for_pt50_surface <- three_opt_result$selected_solutions %>%
  dplyr::filter(Objective_Class == "PT50") %>%
  dplyr::mutate(Objective = "Max PT50")

plot_surface_ru_t0_with_markers(
  surface_df = three_opt_result$surface_df,
  selected_solutions = selected_for_ru_surface
)

plot_surface_pt50_with_markers(
  surface_df = three_opt_result$surface_df,
  selected_solutions = selected_for_pt50_surface
)

plot_selected_weight_profiles_point(three_opt_result$point_curves)
plot_selected_composition_profiles_point(three_opt_result$point_curves)
plot_selected_ed50_profiles_point(three_opt_result$point_curves)
plot_selected_ru_profiles_point(three_opt_result$point_curves)
plot_selected_endpoint_points(three_opt_result$selected_solutions)


# 4. SAVE POINT-ESTIMATE OBJECTS
save(
  weight_csv,
  gcms_csv,
  ed50_csv,
  my_lowers,
  my_uppers,
  my_initial_dose_mg,
  my_time_seq,
  opt_grid_res,
  opt_time_seq_search,
  selection_threshold_frac,
  my_ed50_fit,
  model_selection_result,
  final_evaporation_model,
  file = blockA_model_rdata
)

save(
  opt_grid_res,
  opt_time_seq_search,
  selection_threshold_frac,
  three_opt_result,
  file = blockA_opt_rdata
)
