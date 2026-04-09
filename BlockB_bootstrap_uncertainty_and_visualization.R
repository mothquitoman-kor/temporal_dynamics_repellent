# =====================================================================
# BLOCK B
# BOOTSTRAP FOR 95% CI + VISUALIZATIONS
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

blockA_model_rdata <- "blockA_point_model_results.RData"
blockA_opt_rdata <- "blockA_point_optimization_results.RData"
blockB_boot_rdata <- "blockB_bootstrap_results.RData"

my_boot_B <- 200
bootstrap_seed <- 990202
bootstrap_trace <- 1
verbose_B <- TRUE

# Use the same grid as Block A for final analysis.
# For quick testing only, you may change this, e.g. seq(0, 720, by = 10).
my_time_seq_boot <- NULL

# Keep NULL to include all six selected solutions.
# You may filter using either unique Objective labels (e.g. "PT50_Best")
# or Objective_Class values ("PT50", "RU0").
objectives_to_keep <- NULL

# 1. LOAD POINT-ESTIMATE OBJECTS
load(blockA_model_rdata)
load(blockA_opt_rdata)

if (!exists("my_ed50_fit")) stop("my_ed50_fit missing from blockA model RData")
if (!exists("final_evaporation_model")) stop("final_evaporation_model missing from blockA model RData")
if (!exists("three_opt_result")) stop("three_opt_result missing from blockA optimization RData")

if (is.null(my_time_seq_boot)) {
  my_time_seq_boot <- my_time_seq
}

selected_solutions_use <- three_opt_result$selected_solutions
point_curves_use <- three_opt_result$point_curves

if (!is.null(objectives_to_keep)) {
  keep_solution <- selected_solutions_use$Objective %in% objectives_to_keep
  if ("Objective_Class" %in% names(selected_solutions_use)) {
    keep_solution <- keep_solution | selected_solutions_use$Objective_Class %in% objectives_to_keep
  }

  selected_solutions_use <- selected_solutions_use[keep_solution, , drop = FALSE]

  keep_curve <- point_curves_use$Objective %in% objectives_to_keep
  if ("Objective_Class" %in% names(point_curves_use)) {
    keep_curve <- keep_curve | point_curves_use$Objective_Class %in% objectives_to_keep
  }

  point_curves_use <- point_curves_use[keep_curve, , drop = FALSE]
}

if (nrow(selected_solutions_use) == 0) {
  stop("No selected solutions remained after objectives_to_keep filtering.")
}

# Recompute point curves on the bootstrap grid if needed.
if (!identical(as.numeric(my_time_seq_boot), as.numeric(my_time_seq))) {
  curve_list <- vector("list", nrow(selected_solutions_use))
  summary_list <- vector("list", nrow(selected_solutions_use))

  for (i in seq_len(nrow(selected_solutions_use))) {
    init_comp <- c(
      selected_solutions_use$x_citronellal[i],
      selected_solutions_use$x_citronellol[i],
      selected_solutions_use$x_geraniol[i]
    )
    names(init_comp) <- c("x_citronellal", "x_citronellol", "x_geraniol")

    eval_i <- evaluate_one_solution_from_composition(
      init_comp = init_comp,
      evap_model = final_evaporation_model,
      ed50_fit = my_ed50_fit,
      initial_dose_mg = my_initial_dose_mg,
      time_seq = my_time_seq_boot
    )

    curve_list[[i]] <- eval_i$curve %>%
      dplyr::mutate(
        Objective = selected_solutions_use$Objective[i],
        Objective_Class = if ("Objective_Class" %in% names(selected_solutions_use)) selected_solutions_use$Objective_Class[i] else NA_character_,
        Selection_Rank = if ("Selection_Rank" %in% names(selected_solutions_use)) selected_solutions_use$Selection_Rank[i] else NA_integer_,
        Selection_Label = if ("Selection_Label" %in% names(selected_solutions_use)) selected_solutions_use$Selection_Label[i] else NA_character_
      )

    summary_list[[i]] <- data.frame(
      Objective = selected_solutions_use$Objective[i],
      Objective_Class = if ("Objective_Class" %in% names(selected_solutions_use)) selected_solutions_use$Objective_Class[i] else NA_character_,
      Selection_Rank = if ("Selection_Rank" %in% names(selected_solutions_use)) selected_solutions_use$Selection_Rank[i] else NA_integer_,
      Selection_Label = if ("Selection_Label" %in% names(selected_solutions_use)) selected_solutions_use$Selection_Label[i] else NA_character_,
      Grid_ID = selected_solutions_use$Grid_ID[i],
      x_citronellal = selected_solutions_use$x_citronellal[i],
      x_citronellol = selected_solutions_use$x_citronellol[i],
      x_geraniol = selected_solutions_use$x_geraniol[i],
      Initial_Dose_mg = my_initial_dose_mg,
      ED50_t0 = eval_i$ED50_t0,
      log_ED50_t0 = eval_i$log_ED50_t0,
      RU_t0 = eval_i$RU_t0,
      PT50 = eval_i$PT50,
      log_PT50 = eval_i$log_PT50,
      stringsAsFactors = FALSE
    )
  }

  point_curves_use <- dplyr::bind_rows(curve_list)
  selected_solutions_use <- dplyr::bind_rows(summary_list)
}

# 2. SHARED BOOTSTRAP CI FOR TESTED BLENDS AND FIXED OPTIMIZED SOLUTIONS
shared_boot_result <- joint_bootstrap_shared_fixed_model(
  selected_solutions = selected_solutions_use,
  point_curves = point_curves_use,
  ed50_fit = my_ed50_fit,
  evap_fit = final_evaporation_model,
  weight_csv = weight_csv,
  gcms_csv = gcms_csv,
  initial_dose_mg = my_initial_dose_mg,
  time_seq = my_time_seq_boot,
  B = my_boot_B,
  seed = bootstrap_seed,
  trace = bootstrap_trace,
  verbose = verbose_B
)

all_blend_fixed_boot_result <- shared_boot_result$all_blend_fixed_boot_result
selected_solution_fixed_boot_result <- shared_boot_result$selected_solution_fixed_boot_result

n_success_shared <- length(shared_boot_result$success_ids)
cat("\nSuccessful shared bootstrap replicates:", n_success_shared, "of", my_boot_B, "\n")
cat("Successful all-blend bootstrap replicates:",
    dplyr::n_distinct(all_blend_fixed_boot_result$boot_all$bootstrap_id),
    "of", my_boot_B, "\n")
cat("Successful selected-solution bootstrap replicates:",
    dplyr::n_distinct(selected_solution_fixed_boot_result$boot_all$bootstrap_id),
    "of", my_boot_B, "\n")

cat("\n--- Selected optimized solutions with PT50 CI ---\n")
print_full_df(
  selected_solution_fixed_boot_result$selected_solutions %>%
    dplyr::select(
      dplyr::any_of(c("Objective_Class", "Selection_Rank", "Selection_Label")),
      Objective,
      Grid_ID,
      x_citronellal,
      x_citronellol,
      x_geraniol,
      ED50_t0,
      RU_t0,
      PT50,
      PT50_lower,
      PT50_upper
    )
)

# 4. TESTED-BLEND VISUALIZATIONS
plot_all_blend_dose_ci(all_blend_fixed_boot_result)
plot_all_blend_composition_ci(all_blend_fixed_boot_result)
plot_all_blend_ed50_ci(all_blend_fixed_boot_result)
plot_all_blend_ru_ci(all_blend_fixed_boot_result)
plot_all_blend_pt50_ci(all_blend_fixed_boot_result)

# 5. OPTIMIZED-SOLUTION VISUALIZATIONS
if ("Objective_Class" %in% names(selected_solutions_use)) {
  selected_for_ru_surface <- selected_solutions_use %>%
    dplyr::filter(Objective_Class == "RU0") %>%
    dplyr::mutate(Objective = "Max RU_t0")

  selected_for_pt50_surface <- selected_solutions_use %>%
    dplyr::filter(Objective_Class == "PT50") %>%
    dplyr::mutate(Objective = "Max PT50")
} else {
  selected_for_ru_surface <- selected_solutions_use %>%
    dplyr::filter(grepl("^RU0", Objective)) %>%
    dplyr::mutate(Objective = "Max RU_t0")

  selected_for_pt50_surface <- selected_solutions_use %>%
    dplyr::filter(grepl("^PT50", Objective)) %>%
    dplyr::mutate(Objective = "Max PT50")
}

if (nrow(selected_for_ru_surface) > 0) {
  plot_surface_ru_t0_with_markers(
    surface_df = three_opt_result$surface_df,
    selected_solutions = selected_for_ru_surface
  )
}

if (nrow(selected_for_pt50_surface) > 0) {
  plot_surface_pt50_with_markers(
    surface_df = three_opt_result$surface_df,
    selected_solutions = selected_for_pt50_surface
  )
}

plot_selected_weight_profiles_ci(selected_solution_fixed_boot_result)
plot_selected_composition_profiles_ci(selected_solution_fixed_boot_result)
plot_selected_ed50_profiles_ci(selected_solution_fixed_boot_result)
plot_selected_ru_profiles_ci(selected_solution_fixed_boot_result)
plot_selected_pt50_ci(selected_solution_fixed_boot_result)

# 6. SAVE BOOTSTRAP OBJECTS
save(
  all_blend_fixed_boot_result,
  selected_solution_fixed_boot_result,
  file = blockB_boot_rdata
)
