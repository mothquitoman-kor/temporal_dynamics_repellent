# =====================================================================
# SHARED FUNCTIONS FOR POINT ESTIMATION, OPTIMIZATION, BOOTSTRAP, AND
# VALIDATION SUPPORT
# =====================================================================

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, prob = prob, na.rm = TRUE))
}

all_subsets <- function(x) {
  out <- vector("list", 0)
  for (m in 0:length(x)) {
    if (m == 0) {
      out[[length(out) + 1]] <- character(0)
    } else {
      out <- c(out, combn(x, m, simplify = FALSE))
    }
  }
  out
}

clip_simplex <- function(x, eps = 1e-8) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- eps
  x <- pmax(x, eps)
  x / sum(x)
}

alr_from_simplex <- function(x1, x2, x3, eps = 1e-8) {
  v <- clip_simplex(c(x1, x2, x3), eps = eps)
  c(
    z1 = log(v[1] / v[3]),
    z2 = log(v[2] / v[3])
  )
}

inv_alr_3c <- function(z1, z2) {
  e1 <- exp(z1)
  e2 <- exp(z2)
  denom <- 1 + e1 + e2
  c(
    x_citronellal = e1 / denom,
    x_citronellol = e2 / denom,
    x_geraniol    = 1 / denom
  )
}

format_interaction_set <- function(x) {
  if (length(x) == 0) "None" else paste(x, collapse = ", ")
}

print_full_df <- function(x, digits = 6) {
  old_width <- getOption("width")
  on.exit(options(width = old_width), add = TRUE)
  options(width = 10000)
  print(as.data.frame(x), row.names = FALSE, digits = digits)
}

theme_manuscript <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.3, colour = "grey85"),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot"
    )
}

# DATA HELPERS

prepare_composition_rows <- function(gcms_data) {
  out <- gcms_data %>%
    dplyr::filter(!is.na(citronellal_frac), !is.na(citronellol_frac)) %>%
    dplyr::mutate(
      geraniol_frac = dplyr::if_else(
        !is.na(geraniol_frac),
        geraniol_frac,
        1 - citronellal_frac - citronellol_frac
      )
    )

  if (nrow(out) == 0) {
    out$z1_obs <- numeric(0)
    out$z2_obs <- numeric(0)
    return(out)
  }

  alr_mat <- t(vapply(
    seq_len(nrow(out)),
    FUN = function(i) {
      alr_from_simplex(
        out$citronellal_frac[i],
        out$citronellol_frac[i],
        out$geraniol_frac[i]
      )
    },
    FUN.VALUE = c(z1 = 0, z2 = 0)
  ))
  colnames(alr_mat) <- c("z1", "z2")

  out$z1_obs <- alr_mat[, "z1"]
  out$z2_obs <- alr_mat[, "z2"]
  out
}

extract_blend_initial_conditions <- function(weight_data, gcms_data, blend_id) {
  w0 <- weight_data %>%
    dplyr::filter(Blend_ID == blend_id, Time_min == 0, !is.na(Total_Mass_mg))

  g0 <- gcms_data %>%
    dplyr::filter(
      Blend_ID == blend_id,
      Time_min == 0,
      !is.na(citronellal_frac),
      !is.na(citronellol_frac)
    ) %>%
    dplyr::mutate(
      geraniol_frac = dplyr::if_else(
        !is.na(geraniol_frac),
        geraniol_frac,
        1 - citronellal_frac - citronellol_frac
      )
    )

  if (nrow(w0) == 0) stop(paste("No time-0 weight data for blend:", blend_id))
  if (nrow(g0) == 0) stop(paste("No time-0 composition data for blend:", blend_id))

  initial_mass <- mean(w0$Total_Mass_mg, na.rm = TRUE)
  init_comp_raw <- c(
    mean(g0$citronellal_frac, na.rm = TRUE),
    mean(g0$citronellol_frac, na.rm = TRUE),
    mean(g0$geraniol_frac, na.rm = TRUE)
  )
  init_comp <- clip_simplex(init_comp_raw)

  list(
    initial_mass = initial_mass,
    init_comp = c(
      citronellal = init_comp[1],
      citronellol = init_comp[2],
      geraniol    = init_comp[3]
    )
  )
}

get_initial_composition_for_blend <- function(blend_id, gcms_data) {
  g0 <- gcms_data %>%
    dplyr::filter(
      Blend_ID == blend_id,
      Time_min == 0,
      !is.na(citronellal_frac),
      !is.na(citronellol_frac)
    ) %>%
    dplyr::mutate(
      geraniol_frac = dplyr::if_else(
        !is.na(geraniol_frac),
        geraniol_frac,
        1 - citronellal_frac - citronellol_frac
      )
    )

  if (nrow(g0) == 0) stop("No time-0 composition row found for this blend.")

  init_comp_raw <- c(
    mean(g0$citronellal_frac, na.rm = TRUE),
    mean(g0$citronellol_frac, na.rm = TRUE),
    mean(g0$geraniol_frac, na.rm = TRUE)
  )
  init_comp <- clip_simplex(init_comp_raw)

  c(
    citronellal = init_comp[1],
    citronellol = init_comp[2],
    geraniol    = init_comp[3]
  )
}

# ED50 MODEL

ed50_interaction_terms <- c(
  "x_citronellal:x_citronellol",
  "x_citronellal:x_geraniol",
  "x_citronellol:x_geraniol"
)

build_ed50_formula <- function(retained_interactions) {
  main_terms <- c("x_citronellal", "x_citronellol", "x_geraniol")
  rhs_terms <- c(main_terms, retained_interactions)
  stats::as.formula(
    paste("log(ED50_mg) ~ -1 +", paste(rhs_terms, collapse = " + "))
  )
}

fit_one_ed50_model <- function(ed50_data, retained_interactions) {
  stats::lm(
    formula = build_ed50_formula(retained_interactions),
    data = ed50_data
  )
}

fit_ed50_surface_from_data <- function(ed50_data, verbose = TRUE) {
  required_cols <- c(
    "Blend_ID", "x_citronellal", "x_citronellol", "x_geraniol", "ED50_mg"
  )
  if (!all(required_cols %in% names(ed50_data))) {
    stop("ED50 data is missing required columns.")
  }

  ed50_data <- ed50_data %>%
    dplyr::filter(!is.na(ED50_mg), ED50_mg > 0)

  if (nrow(ed50_data) == 0) stop("No valid ED50 data found.")

  candidate_sets <- all_subsets(ed50_interaction_terms)
  candidate_fits <- vector("list", length(candidate_sets))
  candidate_table <- vector("list", length(candidate_sets))

  for (i in seq_along(candidate_sets)) {
    retained_i <- candidate_sets[[i]]
    fit_i <- fit_one_ed50_model(ed50_data, retained_i)

    logLik_i <- as.numeric(stats::logLik(fit_i))
    n_par_i <- attr(stats::logLik(fit_i), "df")
    aic_i <- stats::AIC(fit_i)
    bic_i <- stats::BIC(fit_i)

    candidate_fits[[i]] <- fit_i
    candidate_table[[i]] <- data.frame(
      model_id = paste0("ED50_model_", i),
      retained_interactions = format_interaction_set(retained_i),
      n_par = n_par_i,
      logLik = logLik_i,
      AIC = aic_i,
      BIC = bic_i,
      stringsAsFactors = FALSE
    )
  }

  candidate_table_df <- dplyr::bind_rows(candidate_table) %>%
    dplyr::arrange(AIC, BIC)

  best_id <- candidate_table_df$model_id[1]
  best_idx <- match(best_id, sapply(candidate_table, function(x) x$model_id[1]))
  best_fit <- candidate_fits[[best_idx]]
  best_retained <- candidate_sets[[best_idx]]
  sigma_log_ed50 <- summary(best_fit)$sigma

  final_summary <- data.frame(
    selected_model_id = best_id,
    retained_interactions = format_interaction_set(best_retained),
    n_par = candidate_table_df$n_par[1],
    logLik = candidate_table_df$logLik[1],
    AIC = candidate_table_df$AIC[1],
    BIC = candidate_table_df$BIC[1],
    sigma_log_ed50 = sigma_log_ed50,
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("\n--- ED50 candidate model table (ranked by AIC) ---\n")
    print(candidate_table_df)
    cat("\n--- Selected ED50 model ---\n")
    print(stats::formula(best_fit))
    print(summary(best_fit))
    print(final_summary)
  }

  list(
    final_model = best_fit,
    final_retained_interactions = best_retained,
    sigma_log_ed50 = sigma_log_ed50,
    candidate_table = candidate_table_df,
    final_summary = final_summary,
    data = ed50_data
  )
}

fit_ed50_surface <- function(ed50_csv, verbose = TRUE) {
  cat("\n--- Fitting quadratic Scheffé ED50 model on log scale ---\n")
  ed50_data <- read.csv(ed50_csv)
  fit_ed50_surface_from_data(ed50_data = ed50_data, verbose = verbose)
}

plot_ed50_surface <- function(model, lower_bounds, upper_bounds) {
  res <- 200

  plot_grid <- expand.grid(
    x_citronellal = seq(lower_bounds[1], upper_bounds[1], length.out = res),
    x_citronellol = seq(lower_bounds[2], upper_bounds[2], length.out = res)
  ) %>%
    dplyr::mutate(x_geraniol = 1 - x_citronellal - x_citronellol) %>%
    dplyr::filter(
      x_geraniol >= lower_bounds[3] & x_geraniol <= upper_bounds[3],
      x_citronellal >= lower_bounds[1] & x_citronellal <= upper_bounds[1],
      x_citronellol >= lower_bounds[2] & x_citronellol <= upper_bounds[2]
    )

  plot_grid$Predicted_ED50 <- exp(stats::predict(model, newdata = plot_grid))

  p <- ggtern::ggtern(
    data = plot_grid,
    ggplot2::aes(x = x_citronellal, y = x_citronellol, z = x_geraniol)
  ) +
    ggplot2::theme_bw() +
    ggtern::theme_noarrows() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 16, color = "black"),
      tern.axis.text = ggplot2::element_text(size = 12, color = "black"),
      legend.position = "right",
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid = ggplot2::element_line(linewidth = 0.6, color = "gray80"),
      tern.axis.line = ggplot2::element_line(linewidth = 1.0, color = "black")
    ) +
    ggplot2::geom_point(ggplot2::aes(color = Predicted_ED50), size = 1.1) +
    ggplot2::scale_color_viridis_c(
      name = "Predicted\nED50 (mg)",
      option = "plasma",
      direction = -1
    ) +
    ggplot2::labs(x = "Citronellal", y = "Citronellol", z = "Geraniol")

  suppressWarnings(print(p))
}

plot_log_ed50_surface <- function(model, lower_bounds, upper_bounds) {
  res <- 200

  plot_grid <- expand.grid(
    x_citronellal = seq(lower_bounds[1], upper_bounds[1], length.out = res),
    x_citronellol = seq(lower_bounds[2], upper_bounds[2], length.out = res)
  ) %>%
    dplyr::mutate(x_geraniol = 1 - x_citronellal - x_citronellol) %>%
    dplyr::filter(
      x_geraniol >= lower_bounds[3] & x_geraniol <= upper_bounds[3],
      x_citronellal >= lower_bounds[1] & x_citronellal <= upper_bounds[1],
      x_citronellol >= lower_bounds[2] & x_citronellol <= upper_bounds[2]
    )

  plot_grid$Predicted_log_ED50 <- stats::predict(model, newdata = plot_grid)

  p <- ggtern::ggtern(
    data = plot_grid,
    ggplot2::aes(x = x_citronellal, y = x_citronellol, z = x_geraniol)
  ) +
    ggplot2::theme_bw() +
    ggtern::theme_noarrows() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 16, color = "black"),
      tern.axis.text = ggplot2::element_text(size = 12, color = "black"),
      legend.position = "right",
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid = ggplot2::element_line(linewidth = 0.6, color = "gray80"),
      tern.axis.line = ggplot2::element_line(linewidth = 1.0, color = "black")
    ) +
    ggplot2::geom_point(ggplot2::aes(color = Predicted_log_ED50), size = 1.1) +
    ggplot2::scale_color_viridis_c(
      name = "Predicted\nlog(ED50)",
      option = "plasma",
      direction = -1
    ) +
    ggplot2::labs(x = "Citronellal", y = "Citronellol", z = "Geraniol")

  suppressWarnings(print(p))
}

# EVAPORATION MODEL

all_param_names <- c(
  "k_citronellal", "k_citronellol", "k_geraniol",
  "A_cit_citol", "A_cit_ger", "A_citol_ger"
)

interaction_names <- c("A_cit_citol", "A_cit_ger", "A_citol_ger")

evaporation_ode_3c <- function(time, state, parameters) {
  m1 <- max(state["m_citronellal"], 0)
  m2 <- max(state["m_citronellol"], 0)
  m3 <- max(state["m_geraniol"], 0)

  total_mass <- m1 + m2 + m3

  if (!is.finite(total_mass) || total_mass <= 1e-12) {
    return(list(c(0, 0, 0)))
  }

  x1 <- m1 / total_mass
  x2 <- m2 / total_mass
  x3 <- m3 / total_mass

  k1 <- parameters["k_citronellal"]
  k2 <- parameters["k_citronellol"]
  k3 <- parameters["k_geraniol"]

  A12 <- parameters["A_cit_citol"]
  A13 <- parameters["A_cit_ger"]
  A23 <- parameters["A_citol_ger"]

  ln_gamma1 <- (x2^2) * A12 + (x3^2) * A13 + x2 * x3 * (A12 + A13 - A23)
  ln_gamma2 <- (x1^2) * A12 + (x3^2) * A23 + x1 * x3 * (A12 + A23 - A13)
  ln_gamma3 <- (x1^2) * A13 + (x2^2) * A23 + x1 * x2 * (A13 + A23 - A12)

  gamma1 <- exp(ln_gamma1)
  gamma2 <- exp(ln_gamma2)
  gamma3 <- exp(ln_gamma3)

  dm1_dt <- -k1 * m1 * gamma1
  dm2_dt <- -k2 * m2 * gamma2
  dm3_dt <- -k3 * m3 * gamma3

  if (!all(is.finite(c(dm1_dt, dm2_dt, dm3_dt)))) {
    return(list(c(0, 0, 0)))
  }

  list(c(dm1_dt, dm2_dt, dm3_dt))
}

build_full_params <- function(free_params, fixed_params = NULL) {
  full <- rep(NA_real_, length(all_param_names))
  names(full) <- all_param_names

  if (!is.null(fixed_params) && length(fixed_params) > 0) {
    full[names(fixed_params)] <- fixed_params
  }

  full[names(free_params)] <- free_params
  full
}

simulate_one_blend <- function(initial_state, sim_times, params_full) {
  out <- tryCatch(
    as.data.frame(
      deSolve::ode(
        y = initial_state,
        times = sim_times,
        func = evaporation_ode_3c,
        parms = params_full,
        method = "lsoda",
        maxsteps = 50000,
        rtol = 1e-8,
        atol = 1e-10
      )
    ),
    error = function(e) NULL
  )

  if (is.null(out)) return(NULL)

  req <- c("time", "m_citronellal", "m_citronellol", "m_geraniol")
  if (!all(req %in% names(out))) return(NULL)
  if (any(!is.finite(as.matrix(out[, req])))) return(NULL)

  out %>%
    dplyr::mutate(
      m_citronellal = pmax(m_citronellal, 0),
      m_citronellol = pmax(m_citronellol, 0),
      m_geraniol    = pmax(m_geraniol, 0),
      Total_Mass_Sim = m_citronellal + m_citronellol + m_geraniol,
      x_cit_sim   = dplyr::if_else(Total_Mass_Sim > 1e-12, m_citronellal / Total_Mass_Sim, NA_real_),
      x_citol_sim = dplyr::if_else(Total_Mass_Sim > 1e-12, m_citronellol / Total_Mass_Sim, NA_real_),
      x_ger_sim   = dplyr::if_else(Total_Mass_Sim > 1e-12, m_geraniol / Total_Mass_Sim, NA_real_)
    )
}

negloglik_evaporation <- function(theta_free,
                                  weight_data,
                                  gcms_data,
                                  fixed_params = NULL,
                                  free_param_names = NULL) {
  names(theta_free) <- free_param_names

  log_sigma_mass <- theta_free["log_sigma_mass"]
  log_sigma_z1   <- theta_free["log_sigma_z1"]
  log_sigma_z2   <- theta_free["log_sigma_z2"]

  sigma_mass <- exp(log_sigma_mass)
  sigma_z1   <- exp(log_sigma_z1)
  sigma_z2   <- exp(log_sigma_z2)

  if (any(!is.finite(c(sigma_mass, sigma_z1, sigma_z2))) ||
      any(c(sigma_mass, sigma_z1, sigma_z2) <= 0)) {
    return(1e12)
  }

  structural_free <- theta_free[
    setdiff(names(theta_free), c("log_sigma_mass", "log_sigma_z1", "log_sigma_z2"))
  ]
  params_full <- build_full_params(structural_free, fixed_params)

  if (any(is.na(params_full))) return(1e12)
  if (any(params_full[c("k_citronellal", "k_citronellol", "k_geraniol")] <= 0)) return(1e12)

  total_nll <- 0
  blend_ids <- union(unique(weight_data$Blend_ID), unique(gcms_data$Blend_ID))

  for (blend in blend_ids) {
    init_info <- extract_blend_initial_conditions(weight_data, gcms_data, blend)

    initial_state <- c(
      init_info$initial_mass * init_info$init_comp["citronellal"],
      init_info$initial_mass * init_info$init_comp["citronellol"],
      init_info$initial_mass * init_info$init_comp["geraniol"]
    )
    names(initial_state) <- c("m_citronellal", "m_citronellol", "m_geraniol")

    w_blend <- weight_data %>%
      dplyr::filter(Blend_ID == blend, !is.na(Total_Mass_mg))

    g_blend <- gcms_data %>%
      dplyr::filter(Blend_ID == blend, !is.na(z1_obs), !is.na(z2_obs))

    sim_times <- sort(unique(c(w_blend$Time_min, g_blend$Time_min, 0)))
    sim_path <- simulate_one_blend(initial_state, sim_times, params_full)

    if (is.null(sim_path)) return(1e12)

    if (nrow(w_blend) > 0) {
      idx_w <- match(w_blend$Time_min, sim_path$time)
      if (any(is.na(idx_w))) return(1e12)

      mu_mass <- sim_path$Total_Mass_Sim[idx_w]
      if (any(!is.finite(mu_mass))) return(1e12)

      total_nll <- total_nll - sum(stats::dnorm(
        x = w_blend$Total_Mass_mg,
        mean = mu_mass,
        sd = sigma_mass,
        log = TRUE
      ))
    }

    if (nrow(g_blend) > 0) {
      idx_g <- match(g_blend$Time_min, sim_path$time)
      if (any(is.na(idx_g))) return(1e12)

      z_pred <- t(vapply(idx_g, function(ii) {
        alr_from_simplex(
          sim_path$x_cit_sim[ii],
          sim_path$x_citol_sim[ii],
          sim_path$x_ger_sim[ii]
        )
      }, FUN.VALUE = c(z1 = 0, z2 = 0)))

      if (any(!is.finite(z_pred))) return(1e12)

      total_nll <- total_nll - sum(stats::dnorm(
        x = g_blend$z1_obs,
        mean = z_pred[, "z1"],
        sd = sigma_z1,
        log = TRUE
      ))

      total_nll <- total_nll - sum(stats::dnorm(
        x = g_blend$z2_obs,
        mean = z_pred[, "z2"],
        sd = sigma_z2,
        log = TRUE
      ))
    }
  }

  if (!is.finite(total_nll)) return(1e12)
  total_nll
}

make_theta_start <- function(weight_data, gcms_data, retained_interactions) {
  w_non0 <- weight_data %>% dplyr::filter(Time_min > 0, !is.na(Total_Mass_mg))
  avg_mass <- stats::sd(w_non0$Total_Mass_mg, na.rm = TRUE)
  if (!is.finite(avg_mass) || avg_mass <= 0) avg_mass <- 1

  gcms_p <- prepare_composition_rows(gcms_data)
  avg_z1 <- stats::sd(gcms_p$z1_obs, na.rm = TRUE)
  avg_z2 <- stats::sd(gcms_p$z2_obs, na.rm = TRUE)
  if (!is.finite(avg_z1) || avg_z1 <= 0) avg_z1 <- 0.2
  if (!is.finite(avg_z2) || avg_z2 <= 0) avg_z2 <- 0.2

  theta <- c(
    k_citronellal = 0.005,
    k_citronellol = 0.002,
    k_geraniol = 0.0025
  )

  if (length(retained_interactions) > 0) {
    theta <- c(theta, stats::setNames(rep(0, length(retained_interactions)), retained_interactions))
  }

  theta <- c(
    theta,
    log_sigma_mass = log(avg_mass),
    log_sigma_z1 = log(avg_z1),
    log_sigma_z2 = log(avg_z2)
  )
  theta
}

make_theta_bounds <- function(retained_interactions) {
  lower <- c(
    k_citronellal = 1e-8,
    k_citronellol = 1e-8,
    k_geraniol = 1e-8
  )
  upper <- c(
    k_citronellal = 0.2,
    k_citronellol = 0.2,
    k_geraniol = 0.2
  )

  if (length(retained_interactions) > 0) {
    lower <- c(lower, stats::setNames(rep(-10, length(retained_interactions)), retained_interactions))
    upper <- c(upper, stats::setNames(rep(10, length(retained_interactions)), retained_interactions))
  }

  lower <- c(lower, log_sigma_mass = log(1e-6), log_sigma_z1 = log(1e-6), log_sigma_z2 = log(1e-6))
  upper <- c(upper, log_sigma_mass = log(1e3),  log_sigma_z1 = log(1e2),  log_sigma_z2 = log(1e2))

  list(lower = lower, upper = upper)
}


fit_evaporation_model_from_data <- function(weight_data,
                                            gcms_data,
                                            fixed_params = NULL,
                                            free_structural = NULL,
                                            init_structural = NULL,
                                            trace = 0) {
  if (is.null(fixed_params)) fixed_params <- numeric(0)

  if (is.null(free_structural)) {
    free_structural <- setdiff(all_param_names, names(fixed_params))
  }

  default_inits <- c(
    k_citronellal = 0.005,
    k_citronellol = 0.002,
    k_geraniol    = 0.001,
    A_cit_citol   = 0,
    A_cit_ger     = 0,
    A_citol_ger   = 0
  )

  if (!is.null(init_structural) && length(init_structural) > 0) {
    default_inits[names(init_structural)] <- init_structural
  }

  free_init <- default_inits[free_structural]

  theta_init <- c(
    free_init,
    log_sigma_mass = log(0.5),
    log_sigma_z1   = log(0.1),
    log_sigma_z2   = log(0.1)
  )

  lower_b <- rep(-Inf, length(theta_init))
  upper_b <- rep( Inf, length(theta_init))
  names(lower_b) <- names(theta_init)
  names(upper_b) <- names(theta_init)

  for (nm in intersect(names(theta_init), c("k_citronellal", "k_citronellol", "k_geraniol"))) {
    lower_b[nm] <- 1e-8
    upper_b[nm] <- 10
  }

  for (nm in intersect(names(theta_init), interaction_names)) {
    lower_b[nm] <- -10
    upper_b[nm] <- 10
  }

  lower_b["log_sigma_mass"] <- log(1e-6)
  upper_b["log_sigma_mass"] <- log(1e3)
  lower_b["log_sigma_z1"] <- log(1e-6)
  upper_b["log_sigma_z1"] <- log(1e3)
  lower_b["log_sigma_z2"] <- log(1e-6)
  upper_b["log_sigma_z2"] <- log(1e3)

  opt_result <- stats::optim(
    par = theta_init,
    fn = negloglik_evaporation,
    weight_data = weight_data,
    gcms_data = gcms_data,
    fixed_params = fixed_params,
    free_param_names = names(theta_init),
    method = "L-BFGS-B",
    lower = lower_b,
    upper = upper_b,
    control = list(trace = trace, maxit = 1000),
    hessian = TRUE
  )

  theta_hat <- opt_result$par
  names(theta_hat) <- names(theta_init)

  structural_hat <- theta_hat[
    setdiff(names(theta_hat), c("log_sigma_mass", "log_sigma_z1", "log_sigma_z2"))
  ]

  sigma_hat <- c(
    sigma_mass = unname(exp(theta_hat["log_sigma_mass"])),
    sigma_z1   = unname(exp(theta_hat["log_sigma_z1"])),
    sigma_z2   = unname(exp(theta_hat["log_sigma_z2"]))
  )

  params_full_hat <- build_full_params(structural_hat, fixed_params)
  logLik_val <- -opt_result$value
  n_par <- length(theta_hat)

  n_mass <- sum(!is.na(weight_data$Total_Mass_mg))
  g_use <- gcms_data %>%
    dplyr::filter(!is.na(z1_obs), !is.na(z2_obs))
  n_comp_rows <- nrow(g_use)
  n_obs <- n_mass + 2 * n_comp_rows

  out <- list(
    estimates_structural = params_full_hat,
    estimates_error = sigma_hat,
    estimates_free = theta_hat,
    fixed_params = fixed_params,
    free_structural = free_structural,
    objective_value = opt_result$value,
    logLik = logLik_val,
    n_par = n_par,
    AIC = -2 * logLik_val + 2 * n_par,
    BIC = -2 * logLik_val + log(n_obs) * n_par,
    convergence = opt_result$convergence,
    hessian = opt_result$hessian,
    n_obs = n_obs
  )

  retained_interactions <- intersect(interaction_names, free_structural)
  out$retained_interactions <- retained_interactions
  out$converged <- isTRUE(opt_result$convergence == 0)
  out
}

fit_one_evaporation_model <- function(weight_data,
                                      gcms_data,
                                      retained_interactions,
                                      trace = 0,
                                      verbose = FALSE) {
  fixed_i <- stats::setNames(
    rep(0, length(setdiff(interaction_names, retained_interactions))),
    setdiff(interaction_names, retained_interactions)
  )

  fit_i <- fit_evaporation_model_from_data(
    weight_data = weight_data,
    gcms_data = prepare_composition_rows(gcms_data),
    fixed_params = fixed_i,
    free_structural = c(
      "k_citronellal", "k_citronellol", "k_geraniol",
      retained_interactions
    ),
    trace = trace
  )

  fit_i$retained_interactions <- retained_interactions
  fit_i$converged <- isTRUE(fit_i$convergence == 0)
  fit_i
}

fit_evaporation_model <- function(weight_csv,
                                  gcms_csv,
                                  fixed_params = NULL,
                                  free_structural = NULL,
                                  init_structural = NULL,
                                  trace = 0) {
  cat("\n--- Fitting mass-proportional evaporation model ---\n")

  w_data <- read.csv(weight_csv)
  g_data <- prepare_composition_rows(read.csv(gcms_csv))

  fit_evaporation_model_from_data(
    weight_data = w_data,
    gcms_data = g_data,
    fixed_params = fixed_params,
    free_structural = free_structural,
    init_structural = init_structural,
    trace = trace
  )
}

select_evaporation_model_from_data <- function(weight_data,
                                               gcms_data,
                                               trace = 0,
                                               verbose = TRUE) {
  candidate_sets <- all_subsets(interaction_names)
  fit_store <- vector("list", length(candidate_sets))
  table_store <- vector("list", length(candidate_sets))

  for (i in seq_along(candidate_sets)) {
    retained_i <- candidate_sets[[i]]
    fixed_i <- stats::setNames(
      rep(0, length(setdiff(interaction_names, retained_i))),
      setdiff(interaction_names, retained_i)
    )

    fit_i <- fit_evaporation_model_from_data(
      weight_data = weight_data,
      gcms_data = gcms_data,
      fixed_params = fixed_i,
      free_structural = c(
        "k_citronellal", "k_citronellol", "k_geraniol",
        retained_i
      ),
      trace = trace
    )

    fit_i$retained_interactions <- retained_i
    fit_i$converged <- isTRUE(fit_i$convergence == 0)

    fit_store[[i]] <- fit_i
    table_store[[i]] <- data.frame(
      model_id = paste0("EVAP_model_", i),
      retained_interactions = format_interaction_set(retained_i),
      dropped_interactions = format_interaction_set(setdiff(interaction_names, retained_i)),
      n_par = fit_i$n_par,
      logLik = fit_i$logLik,
      AIC = fit_i$AIC,
      BIC = fit_i$BIC,
      convergence = fit_i$convergence,
      stringsAsFactors = FALSE
    )
  }

  candidate_table_df <- dplyr::bind_rows(table_store) %>%
    dplyr::filter(convergence == 0) %>%
    dplyr::arrange(AIC, BIC)

  if (nrow(candidate_table_df) == 0) {
    stop("No converged evaporation candidate models were found.")
  }

  best_id <- candidate_table_df$model_id[1]
  best_idx <- match(best_id, sapply(table_store, function(x) x$model_id[1]))
  best_fit <- fit_store[[best_idx]]
  best_retained <- candidate_sets[[best_idx]]

  final_summary <- data.frame(
    selected_model_id = best_id,
    retained_interactions = format_interaction_set(best_retained),
    dropped_interactions = format_interaction_set(setdiff(interaction_names, best_retained)),
    n_par = candidate_table_df$n_par[1],
    logLik = candidate_table_df$logLik[1],
    AIC = candidate_table_df$AIC[1],
    BIC = candidate_table_df$BIC[1],
    convergence = candidate_table_df$convergence[1],
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("\n--- Evaporation candidate model table (converged models only; ranked by AIC) ---\n")
    print(candidate_table_df)
    cat("\n--- Selected evaporation model ---\n")
    print(best_fit$estimates_structural)
    print(best_fit$estimates_error)
    print(final_summary)
  }

  list(
    final_model = best_fit,
    candidate_table = candidate_table_df,
    final_summary = final_summary,
    weight_data = weight_data,
    gcms_data = gcms_data
  )
}

select_evaporation_model <- function(weight_csv,
                                     gcms_csv,
                                     trace = 0,
                                     verbose = TRUE) {
  w_data <- read.csv(weight_csv)
  g_data <- prepare_composition_rows(read.csv(gcms_csv))

  select_evaporation_model_from_data(
    weight_data = w_data,
    gcms_data = g_data,
    trace = trace,
    verbose = verbose
  )
}

simulate_evaporation_curve_from_composition <- function(init_comp,
                                                        evap_model,
                                                        initial_dose_mg,
                                                        time_seq = seq(0, 720, by = 1)) {
  init_comp <- clip_simplex(init_comp)
  names(init_comp) <- c("x_citronellal", "x_citronellol", "x_geraniol")

  initial_state <- c(
    unname(initial_dose_mg * init_comp["x_citronellal"]),
    unname(initial_dose_mg * init_comp["x_citronellol"]),
    unname(initial_dose_mg * init_comp["x_geraniol"])
  )
  names(initial_state) <- c("m_citronellal", "m_citronellol", "m_geraniol")

  sim_path <- simulate_one_blend(
    initial_state = initial_state,
    sim_times = sort(unique(c(0, time_seq))),
    params_full = evap_model$estimates_structural
  )

  if (is.null(sim_path)) stop("Evaporation simulation failed.")

  sim_path %>%
    dplyr::transmute(
      time = time,
      D = Total_Mass_Sim,
      x_citronellal = x_cit_sim,
      x_citronellol = x_citol_sim,
      x_geraniol = x_ger_sim
    )
}

predict_ed50_along_curve <- function(evap_curve, ed50_fit) {
  pred_log <- stats::predict(ed50_fit$final_model, newdata = evap_curve)

  evap_curve %>%
    dplyr::mutate(
      log_ED50 = as.numeric(pred_log),
      ED50 = exp(log_ED50)
    )
}

compute_ru_curve_from_composition <- function(init_comp,
                                              evap_model,
                                              ed50_fit,
                                              initial_dose_mg,
                                              time_seq = seq(0, 720, by = 1)) {
  evap_curve <- simulate_evaporation_curve_from_composition(
    init_comp = init_comp,
    evap_model = evap_model,
    initial_dose_mg = initial_dose_mg,
    time_seq = time_seq
  )

  predict_ed50_along_curve(evap_curve, ed50_fit) %>%
    dplyr::mutate(RU = D / ED50)
}

compute_ru_curve <- function(blend_id,
                             gcms_data,
                             evap_model,
                             ed50_fit,
                             initial_dose_mg,
                             time_seq = seq(0, 720, by = 1)) {
  init_comp <- get_initial_composition_for_blend(blend_id, gcms_data)

  compute_ru_curve_from_composition(
    init_comp = c(
      x_citronellal = unname(init_comp["citronellal"]),
      x_citronellol = unname(init_comp["citronellol"]),
      x_geraniol    = unname(init_comp["geraniol"])
    ),
    evap_model = evap_model,
    ed50_fit = ed50_fit,
    initial_dose_mg = initial_dose_mg,
    time_seq = time_seq
  )
}

compute_pt50 <- function(ru_curve, threshold = 1) {
  above <- ru_curve$RU >= threshold

  if (all(above, na.rm = TRUE)) {
    return(list(
      PT50 = max(ru_curve$time, na.rm = TRUE),
      crossed = FALSE
    ))
  }

  idx <- which(ru_curve$RU < threshold)[1]

  if (idx <= 1) {
    return(list(
      PT50 = ru_curve$time[idx],
      crossed = TRUE
    ))
  }

  t1 <- ru_curve$time[idx - 1]
  t2 <- ru_curve$time[idx]
  r1 <- ru_curve$RU[idx - 1]
  r2 <- ru_curve$RU[idx]

  if (!is.finite(r1) || !is.finite(r2) || r1 == r2) {
    pt <- t2
  } else {
    pt <- t1 + (threshold - r1) * (t2 - t1) / (r2 - r1)
  }

  list(
    PT50 = pt,
    crossed = TRUE
  )
}

# BOOTSTRAP

simulate_bootstrap_ed50_data <- function(ed50_fit) {
  ed50_data <- ed50_fit$data
  mu_log <- stats::predict(ed50_fit$final_model, newdata = ed50_data)
  sigma_log <- ed50_fit$sigma_log_ed50

  sim_log_ed50 <- stats::rnorm(nrow(ed50_data), mean = mu_log, sd = sigma_log)

  ed50_boot <- ed50_data
  ed50_boot$ED50_mg <- exp(sim_log_ed50)
  ed50_boot
}

refit_bootstrap_ed50_fixed_model <- function(ed50_boot_data, ed50_fit_template) {
  fit_i <- fit_one_ed50_model(
    ed50_data = ed50_boot_data,
    retained_interactions = ed50_fit_template$final_retained_interactions
  )

  list(
    final_model = fit_i,
    final_retained_interactions = ed50_fit_template$final_retained_interactions,
    sigma_log_ed50 = summary(fit_i)$sigma,
    data = ed50_boot_data
  )
}

extract_sigma <- function(x, key) {
  if (is.null(x)) return(NA_real_)
  nms <- names(x)
  if (is.null(nms)) return(NA_real_)
  if (key %in% nms) return(unname(x[[key]]))
  hit <- grep(paste0("^", key, "(\\.|$)"), nms, value = TRUE)
  if (length(hit) >= 1) return(unname(x[[hit[1]]]))
  NA_real_
}

simulate_bootstrap_evaporation_data <- function(evap_fit,
                                                weight_data,
                                                gcms_data) {
  sigma_mass <- extract_sigma(evap_fit$estimates_error, "sigma_mass")
  sigma_z1   <- extract_sigma(evap_fit$estimates_error, "sigma_z1")
  sigma_z2   <- extract_sigma(evap_fit$estimates_error, "sigma_z2")

  if (any(!is.finite(c(sigma_mass, sigma_z1, sigma_z2))) ||
      any(c(sigma_mass, sigma_z1, sigma_z2) <= 0)) {
    stop("Bootstrap evaporation sigmas are missing or non-finite.")
  }

  gcms_data <- if (all(c("z1_obs", "z2_obs") %in% names(gcms_data))) gcms_data else prepare_composition_rows(gcms_data)

  w_boot <- weight_data
  g_boot <- gcms_data

  blend_ids <- union(unique(weight_data$Blend_ID), unique(gcms_data$Blend_ID))

  for (blend in blend_ids) {
    init_info <- extract_blend_initial_conditions(weight_data, gcms_data, blend)

    initial_state <- c(
      init_info$initial_mass * unname(init_info$init_comp["citronellal"]),
      init_info$initial_mass * unname(init_info$init_comp["citronellol"]),
      init_info$initial_mass * unname(init_info$init_comp["geraniol"])
    )
    names(initial_state) <- c("m_citronellal", "m_citronellol", "m_geraniol")

    w_blend_idx <- which(w_boot$Blend_ID == blend & !is.na(w_boot$Total_Mass_mg))
    g_blend_idx <- which(g_boot$Blend_ID == blend & !is.na(g_boot$z1_obs) & !is.na(g_boot$z2_obs))

    sim_times <- sort(unique(c(
      weight_data$Time_min[weight_data$Blend_ID == blend & !is.na(weight_data$Total_Mass_mg)],
      gcms_data$Time_min[gcms_data$Blend_ID == blend & !is.na(gcms_data$z1_obs) & !is.na(gcms_data$z2_obs)],
      0
    )))

    sim_path <- simulate_one_blend(
      initial_state = initial_state,
      sim_times = sim_times,
      params_full = evap_fit$estimates_structural
    )

    if (is.null(sim_path)) stop("Bootstrap evaporation simulation failed.")

    if (length(w_blend_idx) > 0) {
      idx_w <- match(w_boot$Time_min[w_blend_idx], sim_path$time)
      if (any(is.na(idx_w))) stop("Bootstrap evaporation simulation times did not match weight data.")
      mu_mass <- sim_path$Total_Mass_Sim[idx_w]
      if (any(!is.finite(mu_mass))) stop("Bootstrap evaporation produced non-finite predicted mass.")
      y_mass <- stats::rnorm(length(mu_mass), mean = mu_mass, sd = sigma_mass)
      w_boot$Total_Mass_mg[w_blend_idx] <- pmax(y_mass, 1e-8)
    }

    if (length(g_blend_idx) > 0) {
      idx_g <- match(g_boot$Time_min[g_blend_idx], sim_path$time)
      if (any(is.na(idx_g))) stop("Bootstrap evaporation simulation times did not match composition data.")

      z_pred <- t(vapply(
        idx_g,
        FUN = function(ii) {
          alr_from_simplex(
            sim_path$x_cit_sim[ii],
            sim_path$x_citol_sim[ii],
            sim_path$x_ger_sim[ii]
          )
        },
        FUN.VALUE = c(z1 = 0, z2 = 0)
      ))
      colnames(z_pred) <- c("z1", "z2")

      if (any(!is.finite(z_pred))) stop("Bootstrap evaporation produced non-finite predicted ALR means.")

      z1_star <- stats::rnorm(nrow(z_pred), mean = z_pred[, "z1"], sd = sigma_z1)
      z2_star <- stats::rnorm(nrow(z_pred), mean = z_pred[, "z2"], sd = sigma_z2)

      comp_star <- t(vapply(
        seq_along(z1_star),
        FUN = function(i) inv_alr_3c(z1_star[i], z2_star[i]),
        FUN.VALUE = c(x_citronellal = 0, x_citronellol = 0, x_geraniol = 0)
      ))
      colnames(comp_star) <- c("x_citronellal", "x_citronellol", "x_geraniol")

      g_boot$citronellal_frac[g_blend_idx] <- comp_star[, "x_citronellal"]
      g_boot$citronellol_frac[g_blend_idx] <- comp_star[, "x_citronellol"]
      if ("geraniol_frac" %in% names(g_boot)) g_boot$geraniol_frac[g_blend_idx] <- comp_star[, "x_geraniol"]
      g_boot$z1_obs[g_blend_idx] <- z1_star
      g_boot$z2_obs[g_blend_idx] <- z2_star
    }
  }

  list(weight_boot = w_boot, gcms_boot = g_boot)
}


refit_bootstrap_evaporation_fixed_model <- function(weight_boot,
                                                    gcms_boot,
                                                    evap_fit_template,
                                                    trace = 0) {
  fit_boot <- fit_evaporation_model_from_data(
    weight_data = weight_boot,
    gcms_data = prepare_composition_rows(gcms_boot),
    fixed_params = evap_fit_template$fixed_params,
    free_structural = evap_fit_template$free_structural,
    init_structural = evap_fit_template$estimates_structural,
    trace = trace
  )

  fit_boot$retained_interactions <- intersect(interaction_names, evap_fit_template$free_structural)
  fit_boot$converged <- isTRUE(fit_boot$convergence == 0)
  fit_boot
}

bootstrap_fixed_model_all_blends_once <- function(ed50_fit,
                                                  evap_fit,
                                                  weight_data,
                                                  gcms_data,
                                                  blend_ids,
                                                  initial_dose_mg,
                                                  time_seq = seq(0, 720, by = 1),
                                                  trace = 0) {
  ed50_boot_data <- simulate_bootstrap_ed50_data(ed50_fit)
  ed50_boot_fit <- refit_bootstrap_ed50_fixed_model(
    ed50_boot_data = ed50_boot_data,
    ed50_fit_template = ed50_fit
  )

  evap_boot_data <- simulate_bootstrap_evaporation_data(
    evap_fit = evap_fit,
    weight_data = weight_data,
    gcms_data = gcms_data
  )

  evap_boot_fit <- refit_bootstrap_evaporation_fixed_model(
    weight_boot = evap_boot_data$weight_boot,
    gcms_boot = evap_boot_data$gcms_boot,
    evap_fit_template = evap_fit,
    trace = trace
  )

  evap_ok <- !is.null(evap_boot_fit) &&
    !is.null(evap_boot_fit$estimates_structural) &&
    all(is.finite(evap_boot_fit$estimates_structural)) &&
    !is.null(evap_boot_fit$estimates_error) &&
    all(is.finite(evap_boot_fit$estimates_error)) &&
    all(evap_boot_fit$estimates_error > 0)

  if (!evap_ok) stop("Bootstrap evaporation refit failed.")

  traj_list <- vector("list", length(blend_ids))

  for (i in seq_along(blend_ids)) {
    blend <- blend_ids[i]

    ru_i <- compute_ru_curve(
      blend_id = blend,
      gcms_data = evap_boot_data$gcms_boot,
      evap_model = evap_boot_fit,
      ed50_fit = ed50_boot_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq
    )

    pt_i <- compute_pt50(ru_i)

    traj_list[[i]] <- ru_i %>%
      dplyr::mutate(
        Blend_ID = blend,
        PT50 = pt_i$PT50,
        crossed = pt_i$crossed
      )
  }

  dplyr::bind_rows(traj_list)
}

joint_bootstrap_all_blends_fixed_model <- function(ed50_fit,
                                                   evap_fit,
                                                   weight_csv = "1_evaporation_weight_template.csv",
                                                   gcms_csv = "2_evaporation_composition_template.csv",
                                                   initial_dose_mg,
                                                   time_seq = seq(0, 720, by = 1),
                                                   B = 200,
                                                   seed = 123,
                                                   trace = 0,
                                                   verbose = TRUE) {
  set.seed(seed)

  weight_data <- read.csv(weight_csv)
  gcms_data <- prepare_composition_rows(read.csv(gcms_csv))
  blend_ids <- sort(unique(weight_data$Blend_ID))

  point_list <- vector("list", length(blend_ids))
  pt_list <- vector("list", length(blend_ids))

  for (i in seq_along(blend_ids)) {
    blend <- blend_ids[i]

    ru_i <- compute_ru_curve(
      blend_id = blend,
      gcms_data = gcms_data,
      evap_model = evap_fit,
      ed50_fit = ed50_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq
    )
    pt_i <- compute_pt50(ru_i)

    point_list[[i]] <- ru_i %>%
      dplyr::mutate(
        Blend_ID = blend,
        PT50 = pt_i$PT50,
        crossed = pt_i$crossed
      )

    pt_list[[i]] <- data.frame(
      Blend_ID = blend,
      PT50_point = pt_i$PT50,
      crossed_point = pt_i$crossed
    )
  }

  point_all <- dplyr::bind_rows(point_list)
  point_pt50 <- dplyr::bind_rows(pt_list)

  boot_all <- vector("list", B)

  for (b in seq_len(B)) {
    if (verbose && (b == 1 || b %% max(1, floor(B / 10)) == 0)) {
      cat("All-blend bootstrap replicate", b, "of", B, "\n")
    }

    res_b <- tryCatch(
      bootstrap_fixed_model_all_blends_once(
        ed50_fit = ed50_fit,
        evap_fit = evap_fit,
        weight_data = weight_data,
        gcms_data = gcms_data,
        blend_ids = blend_ids,
        initial_dose_mg = initial_dose_mg,
        time_seq = time_seq,
        trace = trace
      ),
      error = function(e) { message("Bootstrap replicate ", b, " failed: ", conditionMessage(e)); NULL }
    )

    if (!is.null(res_b)) {
      boot_all[[b]] <- res_b %>%
        dplyr::mutate(bootstrap_id = b)
    }
  }

  boot_all_df <- dplyr::bind_rows(boot_all)
  if (nrow(boot_all_df) == 0) {
    stop("No successful all-blend bootstrap replicates were obtained.")
  }

  traj_ci <- boot_all_df %>%
    dplyr::group_by(Blend_ID, time) %>%
    dplyr::summarise(
      D_lower    = safe_quantile(D, 0.025),
      D_upper    = safe_quantile(D, 0.975),
      x1_lower   = safe_quantile(x_citronellal, 0.025),
      x1_upper   = safe_quantile(x_citronellal, 0.975),
      x2_lower   = safe_quantile(x_citronellol, 0.025),
      x2_upper   = safe_quantile(x_citronellol, 0.975),
      x3_lower   = safe_quantile(x_geraniol, 0.025),
      x3_upper   = safe_quantile(x_geraniol, 0.975),
      log_ED50_lower = safe_quantile(log_ED50, 0.025),
      log_ED50_upper = safe_quantile(log_ED50, 0.975),
      ED50_lower = safe_quantile(ED50, 0.025),
      ED50_upper = safe_quantile(ED50, 0.975),
      RU_lower   = safe_quantile(RU, 0.025),
      RU_upper   = safe_quantile(RU, 0.975),
      PT50_lower = safe_quantile(PT50, 0.025),
      PT50_upper = safe_quantile(PT50, 0.975),
      .groups = "drop"
    )

  pt50_summary <- boot_all_df %>%
    dplyr::group_by(Blend_ID) %>%
    dplyr::summarise(
      PT50_lower = safe_quantile(PT50, 0.025),
      PT50_upper = safe_quantile(PT50, 0.975),
      PT50_median = stats::median(PT50, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(point_pt50, by = "Blend_ID")

  list(
    point_all = point_all,
    traj_ci = traj_ci,
    pt50_summary = pt50_summary,
    boot_all = boot_all_df
  )
}

bootstrap_fixed_model_selected_solutions_once <- function(selected_solutions,
                                                          ed50_fit,
                                                          evap_fit,
                                                          weight_data,
                                                          gcms_data,
                                                          initial_dose_mg,
                                                          time_seq = seq(0, 720, by = 1),
                                                          trace = 0) {
  ed50_boot_data <- simulate_bootstrap_ed50_data(ed50_fit)

  ed50_boot_fit <- refit_bootstrap_ed50_fixed_model(
    ed50_boot_data = ed50_boot_data,
    ed50_fit_template = ed50_fit
  )

  evap_boot_data <- simulate_bootstrap_evaporation_data(
    evap_fit = evap_fit,
    weight_data = weight_data,
    gcms_data = prepare_composition_rows(gcms_data)
  )

  evap_boot_fit <- refit_bootstrap_evaporation_fixed_model(
    weight_boot = evap_boot_data$weight_boot,
    gcms_boot = evap_boot_data$gcms_boot,
    evap_fit_template = evap_fit,
    trace = trace
  )

  if (!isTRUE(evap_boot_fit$converged)) stop("Bootstrap evaporation refit failed.")

  traj_list <- vector("list", nrow(selected_solutions))

  for (i in seq_len(nrow(selected_solutions))) {
    init_comp <- c(
      selected_solutions$x_citronellal[i],
      selected_solutions$x_citronellol[i],
      selected_solutions$x_geraniol[i]
    )
    names(init_comp) <- c("x_citronellal", "x_citronellol", "x_geraniol")

    ru_i <- compute_ru_curve_from_composition(
      init_comp = init_comp,
      evap_model = evap_boot_fit,
      ed50_fit = ed50_boot_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq
    )

    pt_i <- compute_pt50(ru_i)

    traj_list[[i]] <- ru_i %>%
      dplyr::mutate(
        Objective = selected_solutions$Objective[i],
        Grid_ID = selected_solutions$Grid_ID[i],
        x0_citronellal = selected_solutions$x_citronellal[i],
        x0_citronellol = selected_solutions$x_citronellol[i],
        x0_geraniol = selected_solutions$x_geraniol[i],
        PT50 = pt_i$PT50,
        crossed = pt_i$crossed
      )
  }

  dplyr::bind_rows(traj_list)
}

joint_bootstrap_selected_solutions_fixed_model <- function(selected_solutions,
                                                           point_curves,
                                                           ed50_fit,
                                                           evap_fit,
                                                           weight_csv,
                                                           gcms_csv,
                                                           initial_dose_mg,
                                                           time_seq = seq(0, 720, by = 1),
                                                           B = 200,
                                                           seed = 123,
                                                           trace = 0,
                                                           verbose = TRUE) {
  set.seed(seed)

  weight_data <- read.csv(weight_csv)
  gcms_data <- prepare_composition_rows(read.csv(gcms_csv))

  point_pt50 <- selected_solutions %>%
    dplyr::select(
      Objective,
      Grid_ID,
      x_citronellal,
      x_citronellol,
      x_geraniol,
      Initial_Dose_mg,
      ED50_t0,
      log_ED50_t0,
      RU_t0,
      PT50,
      log_PT50
    ) %>%
    dplyr::rename(
      PT50_point = PT50,
      log_PT50_point = log_PT50
    )

  boot_all <- vector("list", B)

  for (b in seq_len(B)) {
    if (verbose && (b == 1 || b %% max(1, floor(B / 10)) == 0)) {
      cat("Selected-solution bootstrap replicate", b, "of", B, "\n")
    }

    res_b <- tryCatch(
      bootstrap_fixed_model_selected_solutions_once(
        selected_solutions = selected_solutions,
        ed50_fit = ed50_fit,
        evap_fit = evap_fit,
        weight_data = weight_data,
        gcms_data = gcms_data,
        initial_dose_mg = initial_dose_mg,
        time_seq = time_seq,
        trace = trace
      ),
      error = function(e) { message("Bootstrap replicate ", b, " failed: ", conditionMessage(e)); NULL }
    )

    if (!is.null(res_b)) {
      boot_all[[b]] <- res_b %>%
        dplyr::mutate(bootstrap_id = b)
    }
  }

  boot_all_df <- dplyr::bind_rows(boot_all)
  if (nrow(boot_all_df) == 0) {
    stop("No successful selected-solution bootstrap replicates were obtained.")
  }

  traj_ci <- boot_all_df %>%
    dplyr::group_by(Objective, time) %>%
    dplyr::summarise(
      D_lower    = safe_quantile(D, 0.025),
      D_upper    = safe_quantile(D, 0.975),
      x1_lower   = safe_quantile(x_citronellal, 0.025),
      x1_upper   = safe_quantile(x_citronellal, 0.975),
      x2_lower   = safe_quantile(x_citronellol, 0.025),
      x2_upper   = safe_quantile(x_citronellol, 0.975),
      x3_lower   = safe_quantile(x_geraniol, 0.025),
      x3_upper   = safe_quantile(x_geraniol, 0.975),
      log_ED50_lower = safe_quantile(log_ED50, 0.025),
      log_ED50_upper = safe_quantile(log_ED50, 0.975),
      ED50_lower = safe_quantile(ED50, 0.025),
      ED50_upper = safe_quantile(ED50, 0.975),
      RU_lower   = safe_quantile(RU, 0.025),
      RU_upper   = safe_quantile(RU, 0.975),
      PT50_lower = safe_quantile(PT50, 0.025),
      PT50_upper = safe_quantile(PT50, 0.975),
      .groups = "drop"
    )

  pt50_summary <- boot_all_df %>%
    dplyr::group_by(Objective) %>%
    dplyr::summarise(
      PT50_lower = safe_quantile(PT50, 0.025),
      PT50_upper = safe_quantile(PT50, 0.975),
      PT50_median = stats::median(PT50, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(point_pt50, by = "Objective")

  selected_solutions_with_ci <- selected_solutions %>%
    dplyr::left_join(
      pt50_summary %>%
        dplyr::select(Objective, PT50_lower, PT50_upper, PT50_median),
      by = "Objective"
    )

  list(
    selected_solutions = selected_solutions_with_ci,
    point_curves = point_curves,
    traj_ci = traj_ci,
    pt50_summary = pt50_summary,
    boot_all = boot_all_df
  )
}

# OPTIMIZATION

generate_feasible_mixture_grid <- function(lower_bounds,
                                           upper_bounds,
                                           res = 100) {
  grid_df <- expand.grid(
    x_citronellal = seq(lower_bounds[1], upper_bounds[1], length.out = res),
    x_citronellol = seq(lower_bounds[2], upper_bounds[2], length.out = res)
  ) %>%
    dplyr::mutate(
      x_geraniol = 1 - x_citronellal - x_citronellol
    ) %>%
    dplyr::filter(
      x_citronellal >= lower_bounds[1],
      x_citronellal <= upper_bounds[1],
      x_citronellol >= lower_bounds[2],
      x_citronellol <= upper_bounds[2],
      x_geraniol >= lower_bounds[3],
      x_geraniol <= upper_bounds[3]
    ) %>%
    dplyr::mutate(
      Grid_ID = paste0("Grid_", dplyr::row_number())
    ) %>%
    dplyr::select(
      Grid_ID,
      x_citronellal,
      x_citronellol,
      x_geraniol
    )

  if (nrow(grid_df) == 0) {
    stop("No feasible mixtures were generated. Check lower/upper bounds.")
  }

  grid_df
}

evaluate_one_solution_from_composition <- function(init_comp,
                                                   evap_model,
                                                   ed50_fit,
                                                   initial_dose_mg,
                                                   time_seq = seq(0, 720, by = 1)) {
  init_comp <- as.numeric(init_comp)
  names(init_comp) <- c("x_citronellal", "x_citronellol", "x_geraniol")
  init_comp <- clip_simplex(init_comp)
  names(init_comp) <- c("x_citronellal", "x_citronellol", "x_geraniol")

  ru_curve <- tryCatch(
    compute_ru_curve_from_composition(
      init_comp = init_comp,
      evap_model = evap_model,
      ed50_fit = ed50_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq
    ),
    error = function(e) NULL
  )

  if (is.null(ru_curve) || nrow(ru_curve) == 0) {
    return(list(
      ED50_t0 = NA_real_,
      log_ED50_t0 = NA_real_,
      RU_t0 = NA_real_,
      PT50 = NA_real_,
      log_PT50 = NA_real_,
      crossed = NA,
      curve = NULL
    ))
  }

  pt <- compute_pt50(ru_curve)
  pt50_val <- pt$PT50

  list(
    ED50_t0 = ru_curve$ED50[1],
    log_ED50_t0 = ru_curve$log_ED50[1],
    RU_t0 = ru_curve$RU[1],
    PT50 = pt50_val,
    log_PT50 = ifelse(is.finite(pt50_val) && pt50_val > 0, log(pt50_val), NA_real_),
    crossed = pt$crossed,
    curve = ru_curve
  )
}

evaluate_mixture_grid <- function(feasible_grid,
                                  evap_model,
                                  ed50_fit,
                                  initial_dose_mg,
                                  time_seq = seq(0, 720, by = 1),
                                  verbose = TRUE) {
  n_grid <- nrow(feasible_grid)
  result_list <- vector("list", n_grid)

  for (i in seq_len(n_grid)) {
    if (verbose && (i == 1 || i %% max(1, floor(n_grid / 10)) == 0)) {
      cat("Evaluating candidate", i, "of", n_grid, "\n")
    }

    init_comp <- c(
      feasible_grid$x_citronellal[i],
      feasible_grid$x_citronellol[i],
      feasible_grid$x_geraniol[i]
    )
    names(init_comp) <- c("x_citronellal", "x_citronellol", "x_geraniol")

    eval_i <- evaluate_one_solution_from_composition(
      init_comp = init_comp,
      evap_model = evap_model,
      ed50_fit = ed50_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq
    )

    result_list[[i]] <- data.frame(
      Grid_ID = feasible_grid$Grid_ID[i],
      x_citronellal = feasible_grid$x_citronellal[i],
      x_citronellol = feasible_grid$x_citronellol[i],
      x_geraniol = feasible_grid$x_geraniol[i],
      Initial_Dose_mg = initial_dose_mg,
      ED50_t0 = eval_i$ED50_t0,
      log_ED50_t0 = eval_i$log_ED50_t0,
      RU_t0 = eval_i$RU_t0,
      PT50 = eval_i$PT50,
      log_PT50 = eval_i$log_PT50,
      crossed = eval_i$crossed,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(result_list)
}

find_three_optima <- function(lower_bounds,
                              upper_bounds,
                              evap_model,
                              ed50_fit,
                              initial_dose_mg,
                              grid_res = 80,
                              time_seq_search = seq(0, 720, by = 10),
                              time_seq_final = seq(0, 720, by = 1),
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

  ed50_row <- surface_df %>%
    dplyr::filter(is.finite(ED50_t0)) %>%
    dplyr::slice_min(order_by = ED50_t0, n = 1, with_ties = FALSE)

  ru_row <- surface_df %>%
    dplyr::filter(is.finite(RU_t0)) %>%
    dplyr::slice_max(order_by = RU_t0, n = 1, with_ties = FALSE)

  pt50_row <- surface_df %>%
    dplyr::filter(is.finite(PT50)) %>%
    dplyr::slice_max(order_by = PT50, n = 1, with_ties = FALSE)

  if (nrow(ed50_row) == 0) stop("ED50(t=0) optimum could not be identified.")
  if (nrow(ru_row) == 0) stop("RU_t0 optimum could not be identified.")
  if (nrow(pt50_row) == 0) stop("PT50 optimum could not be identified.")

  selected_df <- dplyr::bind_rows(
    ed50_row %>% dplyr::mutate(Objective = "Min ED50 at t=0"),
    ru_row %>% dplyr::mutate(Objective = "Max RU_t0"),
    pt50_row %>% dplyr::mutate(Objective = "Max PT50")
  ) %>%
    dplyr::select(
      Objective,
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
      dplyr::mutate(Objective = selected_df$Objective[i])

    summary_list[[i]] <- data.frame(
      Objective = selected_df$Objective[i],
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
    cat("\n--- Three selected solutions (point estimates) ---\n")
    print_full_df(point_summary)
  }

  list(
    feasible_grid = feasible_grid,
    surface_df = surface_df,
    selected_solutions = point_summary,
    point_curves = point_curves
  )
}

# PLOTTING

build_selected_marker_df <- function(selected_solutions) {
  selected_solutions %>%
    dplyr::select(
      Objective,
      x_citronellal,
      x_citronellol,
      x_geraniol
    )
}

plot_surface_generic <- function(surface_df,
                                 selected_solutions,
                                 fill_var,
                                 fill_lab,
                                 objective_name,
                                 title_text) {
  marker_df <- build_selected_marker_df(selected_solutions) %>%
    dplyr::filter(Objective == objective_name)

  p <- ggtern::ggtern(
    data = surface_df,
    ggplot2::aes(x = x_citronellal, y = x_citronellol, z = x_geraniol)
  ) +
    ggplot2::theme_bw() +
    ggtern::theme_noarrows() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 16, color = "black"),
      tern.axis.text = ggplot2::element_text(size = 12, color = "black"),
      legend.position = "right",
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid = ggplot2::element_line(linewidth = 0.6, color = "gray80"),
      tern.axis.line = ggplot2::element_line(linewidth = 1.0, color = "black")
    ) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[fill_var]]), size = 1.1) +
    ggplot2::scale_color_viridis_c(
      name = fill_lab,
      option = "plasma",
      direction = -1
    ) +
    ggplot2::geom_point(
      data = marker_df,
      size = 4,
      shape = 21,
      color = "black",
      fill = "white",
      stroke = 1.2
    ) +
    ggplot2::labs(
      x = "Citronellal",
      y = "Citronellol",
      z = "Geraniol",
      title = title_text
    )

  print(p)
  invisible(p)
}

plot_surface_ed50_t0_with_markers <- function(surface_df, selected_solutions) {
  plot_surface_generic(
    surface_df = surface_df,
    selected_solutions = selected_solutions,
    fill_var = "ED50_t0",
    fill_lab = "Predicted\nED50 at t=0",
    objective_name = "Min ED50 at t=0",
    title_text = "ED50(t=0) surface with selected solution"
  )
}

plot_surface_ru_t0_with_markers <- function(surface_df, selected_solutions) {
  plot_surface_generic(
    surface_df = surface_df,
    selected_solutions = selected_solutions,
    fill_var = "RU_t0",
    fill_lab = "Predicted\nRU at t=0",
    objective_name = "Max RU_t0",
    title_text = "RU(t=0) surface with selected solution"
  )
}

plot_surface_pt50_with_markers <- function(surface_df, selected_solutions) {
  plot_surface_generic(
    surface_df = surface_df,
    selected_solutions = selected_solutions,
    fill_var = "PT50",
    fill_lab = "Predicted\nPT50",
    objective_name = "Max PT50",
    title_text = "PT50 surface with selected solution"
  )
}

plot_selected_weight_profiles_point <- function(point_curves) {
  p <- ggplot2::ggplot(point_curves, ggplot2::aes(x = time / 60, y = D, color = Objective)) +
    ggplot2::geom_line(linewidth = 1.1) +
    theme_manuscript() +
    ggplot2::labs(x = "Time (hours)", y = "Remaining total weight (mg)", title = "Evaporation profile: total weight")
  print(p)
}

plot_selected_composition_profiles_point <- function(point_curves) {
  plot_df <- point_curves %>%
    dplyr::select(Objective, time, x_citronellal, x_citronellol, x_geraniol) %>%
    tidyr::pivot_longer(
      cols = c(x_citronellal, x_citronellol, x_geraniol),
      names_to = "Component",
      values_to = "Value"
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = Value, color = Component)) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::facet_wrap(~ Objective, scales = "fixed") +
    theme_manuscript() +
    ggplot2::labs(x = "Time (hours)", y = "Composition / mass fraction", title = "Evaporation profile: composition")
  print(p)
}

plot_selected_ed50_profiles_point <- function(point_curves) {
  p <- ggplot2::ggplot(point_curves, ggplot2::aes(x = time / 60, y = ED50, color = Objective)) +
    ggplot2::geom_line(linewidth = 1.1) +
    theme_manuscript() +
    ggplot2::labs(x = "Time (hours)", y = "ED50(t)", title = "ED50 trajectories for the selected solutions")
  print(p)
}

plot_selected_ru_profiles_point <- function(point_curves) {
  p <- ggplot2::ggplot(point_curves, ggplot2::aes(x = time / 60, y = RU, color = Objective)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_line(linewidth = 1.1) +
    theme_manuscript() +
    ggplot2::labs(x = "Time (hours)", y = "RU(t)", title = "RU trajectories for the selected solutions")
  print(p)
}

plot_selected_endpoint_points <- function(selected_solutions) {
  plot_df <- selected_solutions %>%
    dplyr::select(Objective, RU_t0, PT50) %>%
    tidyr::pivot_longer(
      cols = c(RU_t0, PT50),
      names_to = "Endpoint",
      values_to = "Value"
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Objective, y = Value)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~ Endpoint, scales = "free_y") +
    theme_manuscript() +
    ggplot2::labs(x = NULL, y = "Value", title = "Selected-solution endpoints") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  print(p)
}

bootstrap_fixed_model_shared_once <- function(selected_solutions,
                                              ed50_fit,
                                              evap_fit,
                                              weight_data,
                                              gcms_data,
                                              blend_ids,
                                              initial_dose_mg,
                                              time_seq = seq(0, 720, by = 1),
                                              trace = 0) {
  ed50_boot_data <- simulate_bootstrap_ed50_data(ed50_fit)
  ed50_boot_fit <- refit_bootstrap_ed50_fixed_model(
    ed50_boot_data = ed50_boot_data,
    ed50_fit_template = ed50_fit
  )

  evap_boot_data <- simulate_bootstrap_evaporation_data(
    evap_fit = evap_fit,
    weight_data = weight_data,
    gcms_data = gcms_data
  )

  evap_boot_fit <- refit_bootstrap_evaporation_fixed_model(
    weight_boot = evap_boot_data$weight_boot,
    gcms_boot = evap_boot_data$gcms_boot,
    evap_fit_template = evap_fit,
    trace = trace
  )

  evap_ok <- !is.null(evap_boot_fit) &&
    !is.null(evap_boot_fit$estimates_structural) &&
    all(is.finite(evap_boot_fit$estimates_structural)) &&
    !is.null(evap_boot_fit$estimates_error) &&
    all(is.finite(evap_boot_fit$estimates_error)) &&
    all(evap_boot_fit$estimates_error > 0)

  if (!evap_ok) stop("Bootstrap evaporation refit failed.")

  tested_list <- vector("list", length(blend_ids))
  for (i in seq_along(blend_ids)) {
    blend <- blend_ids[i]

    ru_i <- compute_ru_curve(
      blend_id = blend,
      gcms_data = evap_boot_data$gcms_boot,
      evap_model = evap_boot_fit,
      ed50_fit = ed50_boot_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq
    )

    pt_i <- compute_pt50(ru_i)

    tested_list[[i]] <- ru_i %>%
      dplyr::mutate(
        Blend_ID = blend,
        PT50 = pt_i$PT50,
        crossed = pt_i$crossed
      )
  }

  selected_list <- vector("list", nrow(selected_solutions))
  for (i in seq_len(nrow(selected_solutions))) {
    init_comp <- c(
      selected_solutions$x_citronellal[i],
      selected_solutions$x_citronellol[i],
      selected_solutions$x_geraniol[i]
    )
    names(init_comp) <- c("x_citronellal", "x_citronellol", "x_geraniol")

    ru_i <- compute_ru_curve_from_composition(
      init_comp = init_comp,
      evap_model = evap_boot_fit,
      ed50_fit = ed50_boot_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq
    )

    pt_i <- compute_pt50(ru_i)

    selected_list[[i]] <- ru_i %>%
      dplyr::mutate(
        Objective = selected_solutions$Objective[i],
        Grid_ID = selected_solutions$Grid_ID[i],
        x0_citronellal = selected_solutions$x_citronellal[i],
        x0_citronellol = selected_solutions$x_citronellol[i],
        x0_geraniol = selected_solutions$x_geraniol[i],
        PT50 = pt_i$PT50,
        crossed = pt_i$crossed
      )
  }

  list(
    tested_blends = dplyr::bind_rows(tested_list),
    selected_solutions = dplyr::bind_rows(selected_list)
  )
}

joint_bootstrap_shared_fixed_model <- function(selected_solutions,
                                               point_curves,
                                               ed50_fit,
                                               evap_fit,
                                               weight_csv = "1_evaporation_weight_template.csv",
                                               gcms_csv = "2_evaporation_composition_template.csv",
                                               initial_dose_mg,
                                               time_seq = seq(0, 720, by = 1),
                                               B = 200,
                                               seed = 123,
                                               trace = 0,
                                               verbose = TRUE) {
  set.seed(seed)

  weight_data <- read.csv(weight_csv)
  gcms_data <- prepare_composition_rows(read.csv(gcms_csv))
  blend_ids <- sort(unique(weight_data$Blend_ID))

  point_list <- vector("list", length(blend_ids))
  pt_list <- vector("list", length(blend_ids))

  for (i in seq_along(blend_ids)) {
    blend <- blend_ids[i]

    ru_i <- compute_ru_curve(
      blend_id = blend,
      gcms_data = gcms_data,
      evap_model = evap_fit,
      ed50_fit = ed50_fit,
      initial_dose_mg = initial_dose_mg,
      time_seq = time_seq
    )
    pt_i <- compute_pt50(ru_i)

    point_list[[i]] <- ru_i %>%
      dplyr::mutate(
        Blend_ID = blend,
        PT50 = pt_i$PT50,
        crossed = pt_i$crossed
      )

    pt_list[[i]] <- data.frame(
      Blend_ID = blend,
      PT50_point = pt_i$PT50,
      crossed_point = pt_i$crossed
    )
  }

  point_all <- dplyr::bind_rows(point_list)
  point_pt50 <- dplyr::bind_rows(pt_list)

  selected_point_pt50 <- selected_solutions %>%
    dplyr::select(
      Objective,
      Grid_ID,
      x_citronellal,
      x_citronellol,
      x_geraniol,
      Initial_Dose_mg,
      ED50_t0,
      log_ED50_t0,
      RU_t0,
      PT50,
      log_PT50
    ) %>%
    dplyr::rename(
      PT50_point = PT50,
      log_PT50_point = log_PT50
    )

  tested_boot <- vector("list", B)
  selected_boot <- vector("list", B)
  success_ids <- integer(0)

  for (b in seq_len(B)) {
    if (verbose && (b == 1 || b %% max(1, floor(B / 10)) == 0)) {
      cat("Shared bootstrap replicate", b, "of", B, "\n")
    }

    res_b <- tryCatch(
      bootstrap_fixed_model_shared_once(
        selected_solutions = selected_solutions,
        ed50_fit = ed50_fit,
        evap_fit = evap_fit,
        weight_data = weight_data,
        gcms_data = gcms_data,
        blend_ids = blend_ids,
        initial_dose_mg = initial_dose_mg,
        time_seq = time_seq,
        trace = trace
      ),
      error = function(e) { message("Bootstrap replicate ", b, " failed: ", conditionMessage(e)); NULL }
    )

    if (!is.null(res_b)) {
      success_ids <- c(success_ids, b)
      tested_boot[[b]] <- res_b$tested_blends %>%
        dplyr::mutate(bootstrap_id = b)
      selected_boot[[b]] <- res_b$selected_solutions %>%
        dplyr::mutate(bootstrap_id = b)
    }
  }

  tested_boot_df <- dplyr::bind_rows(tested_boot)
  selected_boot_df <- dplyr::bind_rows(selected_boot)

  if (nrow(tested_boot_df) == 0 || nrow(selected_boot_df) == 0) {
    stop("No successful shared bootstrap replicates were obtained.")
  }

  tested_traj_ci <- tested_boot_df %>%
    dplyr::group_by(Blend_ID, time) %>%
    dplyr::summarise(
      D_lower    = safe_quantile(D, 0.025),
      D_upper    = safe_quantile(D, 0.975),
      x1_lower   = safe_quantile(x_citronellal, 0.025),
      x1_upper   = safe_quantile(x_citronellal, 0.975),
      x2_lower   = safe_quantile(x_citronellol, 0.025),
      x2_upper   = safe_quantile(x_citronellol, 0.975),
      x3_lower   = safe_quantile(x_geraniol, 0.025),
      x3_upper   = safe_quantile(x_geraniol, 0.975),
      log_ED50_lower = safe_quantile(log_ED50, 0.025),
      log_ED50_upper = safe_quantile(log_ED50, 0.975),
      ED50_lower = safe_quantile(ED50, 0.025),
      ED50_upper = safe_quantile(ED50, 0.975),
      RU_lower   = safe_quantile(RU, 0.025),
      RU_upper   = safe_quantile(RU, 0.975),
      PT50_lower = safe_quantile(PT50, 0.025),
      PT50_upper = safe_quantile(PT50, 0.975),
      .groups = "drop"
    )

  tested_pt50_summary <- tested_boot_df %>%
    dplyr::group_by(Blend_ID) %>%
    dplyr::summarise(
      PT50_lower = safe_quantile(PT50, 0.025),
      PT50_upper = safe_quantile(PT50, 0.975),
      PT50_median = stats::median(PT50, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(point_pt50, by = "Blend_ID")

  selected_traj_ci <- selected_boot_df %>%
    dplyr::group_by(Objective, time) %>%
    dplyr::summarise(
      D_lower    = safe_quantile(D, 0.025),
      D_upper    = safe_quantile(D, 0.975),
      x1_lower   = safe_quantile(x_citronellal, 0.025),
      x1_upper   = safe_quantile(x_citronellal, 0.975),
      x2_lower   = safe_quantile(x_citronellol, 0.025),
      x2_upper   = safe_quantile(x_citronellol, 0.975),
      x3_lower   = safe_quantile(x_geraniol, 0.025),
      x3_upper   = safe_quantile(x_geraniol, 0.975),
      log_ED50_lower = safe_quantile(log_ED50, 0.025),
      log_ED50_upper = safe_quantile(log_ED50, 0.975),
      ED50_lower = safe_quantile(ED50, 0.025),
      ED50_upper = safe_quantile(ED50, 0.975),
      RU_lower   = safe_quantile(RU, 0.025),
      RU_upper   = safe_quantile(RU, 0.975),
      PT50_lower = safe_quantile(PT50, 0.025),
      PT50_upper = safe_quantile(PT50, 0.975),
      .groups = "drop"
    )

  selected_pt50_summary <- selected_boot_df %>%
    dplyr::group_by(Objective, Grid_ID, x0_citronellal, x0_citronellol, x0_geraniol) %>%
    dplyr::summarise(
      PT50_lower = safe_quantile(PT50, 0.025),
      PT50_upper = safe_quantile(PT50, 0.975),
      PT50_median = stats::median(PT50, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      selected_point_pt50,
      by = c(
        "Objective",
        "Grid_ID",
        "x0_citronellal" = "x_citronellal",
        "x0_citronellol" = "x_citronellol",
        "x0_geraniol" = "x_geraniol"
      )
    )

  selected_solutions_out <- selected_solutions %>%
    dplyr::left_join(
      selected_pt50_summary %>%
        dplyr::select(Objective, Grid_ID, PT50_lower, PT50_upper, PT50_median),
      by = c("Objective", "Grid_ID")
    )

  list(
    all_blend_fixed_boot_result = list(
      point_all = point_all,
      traj_ci = tested_traj_ci,
      pt50_summary = tested_pt50_summary,
      boot_all = tested_boot_df
    ),
    selected_solution_fixed_boot_result = list(
      selected_solutions = selected_solutions_out,
      point_curves = point_curves,
      traj_ci = selected_traj_ci,
      pt50_summary = selected_pt50_summary,
      boot_all = selected_boot_df
    ),
    success_ids = success_ids
  )
}


plot_all_blend_dose_ci <- function(all_blend_boot_result) {
  plot_df <- all_blend_boot_result$point_all %>%
    dplyr::select(Blend_ID, time, D) %>%
    dplyr::left_join(
      all_blend_boot_result$traj_ci %>%
        dplyr::select(Blend_ID, time, D_lower, D_upper),
      by = c("Blend_ID", "time")
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = D)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = D_lower, ymax = D_upper), alpha = 0.18) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~ Blend_ID, scales = "free_y") +
    theme_manuscript() +
    ggplot2::labs(x = "Time (h)", y = "Remaining dose (mg)", title = "Dose trajectories of experimental blends")
  print(p)
}

plot_all_blend_composition_ci <- function(all_blend_boot_result) {
  point_long <- all_blend_boot_result$point_all %>%
    dplyr::select(Blend_ID, time, x_citronellal, x_citronellol, x_geraniol) %>%
    tidyr::pivot_longer(
      cols = c(x_citronellal, x_citronellol, x_geraniol),
      names_to = "Component",
      values_to = "Value"
    )

  ci_long <- all_blend_boot_result$traj_ci %>%
    dplyr::select(Blend_ID, time, x1_lower, x1_upper, x2_lower, x2_upper, x3_lower, x3_upper) %>%
    tidyr::pivot_longer(
      cols = c(x1_lower, x1_upper, x2_lower, x2_upper, x3_lower, x3_upper),
      names_to = c("Component", ".value"),
      names_pattern = "x([123])_(lower|upper)"
    ) %>%
    dplyr::mutate(
      Component = dplyr::recode(Component, "1" = "x_citronellal", "2" = "x_citronellol", "3" = "x_geraniol")
    )

  plot_df <- point_long %>%
    dplyr::left_join(ci_long, by = c("Blend_ID", "time", "Component")) %>%
    dplyr::mutate(
      Component = factor(Component,
        levels = c("x_citronellal", "x_citronellol", "x_geraniol"),
        labels = c("Citronellal", "Citronellol", "Geraniol")
      )
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = Value, color = Component, fill = Component)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~ Blend_ID, scales = "fixed") +
    theme_manuscript() +
    ggplot2::labs(x = "Time (h)", y = "Mass fraction", title = "Composition trajectories of experimental blends")
  print(p)
}

plot_all_blend_ed50_ci <- function(all_blend_boot_result) {
  plot_df <- all_blend_boot_result$point_all %>%
    dplyr::select(Blend_ID, time, ED50) %>%
    dplyr::left_join(
      all_blend_boot_result$traj_ci %>%
        dplyr::select(Blend_ID, time, ED50_lower, ED50_upper),
      by = c("Blend_ID", "time")
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = ED50)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ED50_lower, ymax = ED50_upper), alpha = 0.18) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~ Blend_ID, scales = "free_y") +
    theme_manuscript() +
    ggplot2::labs(x = "Time (h)", y = expression(ED[50](t)~"(mg)"), title = expression(ED[50](t)~" trajectories of experimental blends"))
  print(p)
}

plot_all_blend_ru_ci <- function(all_blend_boot_result) {
  plot_df <- all_blend_boot_result$point_all %>%
    dplyr::select(Blend_ID, time, RU) %>%
    dplyr::left_join(
      all_blend_boot_result$traj_ci %>%
        dplyr::select(Blend_ID, time, RU_lower, RU_upper),
      by = c("Blend_ID", "time")
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = RU)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = RU_lower, ymax = RU_upper), alpha = 0.18) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(~ Blend_ID, scales = "free_y") +
    theme_manuscript() +
    ggplot2::labs(x = "Time (h)", y = "RU(t)", title = "Repellent-unit trajectories of experimental blends")
  print(p)
}

plot_all_blend_pt50_ci <- function(all_blend_boot_result) {
  plot_df <- all_blend_boot_result$pt50_summary %>%
    dplyr::arrange(dplyr::desc(PT50_point)) %>%
    dplyr::mutate(Blend_ID = factor(Blend_ID, levels = Blend_ID))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Blend_ID, y = PT50_point)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = PT50_lower, ymax = PT50_upper), width = 0.2, linewidth = 0.6) +
    ggplot2::geom_point(size = 2.5) +
    theme_manuscript() +
    ggplot2::labs(x = "Blend ID", y = "PT50 (min)", title = "PT50 of experimental blends") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  print(p)
}

plot_selected_weight_profiles_ci <- function(selected_solution_boot_result) {
  plot_df <- selected_solution_boot_result$point_curves %>%
    dplyr::select(Objective, time, D) %>%
    dplyr::left_join(
      selected_solution_boot_result$traj_ci %>%
        dplyr::select(Objective, time, D_lower, D_upper),
      by = c("Objective", "time")
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = D, color = Objective, fill = Objective)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = D_lower, ymax = D_upper), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 1.1) +
    theme_manuscript() +
    ggplot2::labs(x = "Time (hours)", y = "Remaining total weight (mg)", title = "Predicted dose trajectories with bootstrap 95% CI")
  print(p)
}

plot_selected_composition_profiles_ci <- function(selected_solution_boot_result) {
  point_long <- selected_solution_boot_result$point_curves %>%
    dplyr::select(Objective, time, x_citronellal, x_citronellol, x_geraniol) %>%
    tidyr::pivot_longer(
      cols = c(x_citronellal, x_citronellol, x_geraniol),
      names_to = "Component",
      values_to = "Value"
    )

  ci_long <- selected_solution_boot_result$traj_ci %>%
    dplyr::select(Objective, time, x1_lower, x1_upper, x2_lower, x2_upper, x3_lower, x3_upper) %>%
    tidyr::pivot_longer(
      cols = c(x1_lower, x1_upper, x2_lower, x2_upper, x3_lower, x3_upper),
      names_to = c("Component", ".value"),
      names_pattern = "x([123])_(lower|upper)"
    ) %>%
    dplyr::mutate(
      Component = dplyr::recode(Component, "1" = "x_citronellal", "2" = "x_citronellol", "3" = "x_geraniol")
    )

  plot_df <- point_long %>%
    dplyr::left_join(ci_long, by = c("Objective", "time", "Component")) %>%
    dplyr::mutate(
      Component = factor(Component,
        levels = c("x_citronellal", "x_citronellol", "x_geraniol"),
        labels = c("Citronellal", "Citronellol", "Geraniol")
      )
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = Value, color = Component, fill = Component)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::facet_wrap(~ Objective, scales = "fixed") +
    theme_manuscript() +
    ggplot2::labs(x = "Time (hours)", y = "Composition / mass fraction", title = "Predicted composition trajectories with bootstrap 95% CI")
  print(p)
}

plot_selected_ed50_profiles_ci <- function(selected_solution_boot_result) {
  plot_df <- selected_solution_boot_result$point_curves %>%
    dplyr::select(Objective, time, ED50) %>%
    dplyr::left_join(
      selected_solution_boot_result$traj_ci %>%
        dplyr::select(Objective, time, ED50_lower, ED50_upper),
      by = c("Objective", "time")
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = ED50, color = Objective, fill = Objective)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ED50_lower, ymax = ED50_upper), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 1.1) +
    theme_manuscript() +
    ggplot2::labs(x = "Time (hours)", y = "ED50(t)", title = "Predicted ED50 trajectories with bootstrap 95% CI")
  print(p)
}

plot_selected_ru_profiles_ci <- function(selected_solution_boot_result) {
  plot_df <- selected_solution_boot_result$point_curves %>%
    dplyr::select(Objective, time, RU) %>%
    dplyr::left_join(
      selected_solution_boot_result$traj_ci %>%
        dplyr::select(Objective, time, RU_lower, RU_upper),
      by = c("Objective", "time")
    )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time / 60, y = RU, color = Objective, fill = Objective)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = RU_lower, ymax = RU_upper), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 1.1) +
    theme_manuscript() +
    ggplot2::labs(x = "Time (hours)", y = "RU(t)", title = "Predicted RU trajectories with bootstrap 95% CI")
  print(p)
}

plot_selected_pt50_ci <- function(selected_solution_boot_result) {
  plot_df <- selected_solution_boot_result$pt50_summary %>%
    dplyr::mutate(Objective = factor(Objective, levels = Objective))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Objective, y = PT50_point)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = PT50_lower, ymax = PT50_upper), width = 0.15, linewidth = 0.7) +
    ggplot2::geom_point(size = 3) +
    theme_manuscript() +
    ggplot2::labs(x = NULL, y = "PT50 (min)", title = "Predicted PT50 with bootstrap 95% CI") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  print(p)
}

# =====================================================================
# VALIDATION HELPERS
# =====================================================================

resolve_validation_dose <- function(selected_df,
                                    dose_mode = c("initial_dose", "predicted_ed50_t0", "manual"),
                                    ed50_multiplier = 1,
                                    manual_validation_dose_mg = NULL) {
  dose_mode <- match.arg(dose_mode)

  if (dose_mode == "initial_dose") {
    out <- selected_df$Initial_Dose_mg
  } else if (dose_mode == "predicted_ed50_t0") {
    out <- selected_df$ED50_t0 * ed50_multiplier
  } else {
    if (is.null(manual_validation_dose_mg)) {
      stop("manual_validation_dose_mg must be supplied when dose_mode = 'manual'.")
    }

    if (length(manual_validation_dose_mg) == 1) {
      out <- rep(manual_validation_dose_mg, nrow(selected_df))
    } else if (length(manual_validation_dose_mg) == nrow(selected_df)) {
      out <- manual_validation_dose_mg
    } else {
      stop("manual_validation_dose_mg must have length 1 or nrow(selected_df).")
    }
  }

  if (any(!is.finite(out)) || any(out <= 0)) {
    stop("Validation doses must be positive finite numbers.")
  }

  out
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

bias_mean <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  mean(obs[ok] - pred[ok])
}

wilson_ci <- function(x, n, conf.level = 0.95) {
  if (!is.finite(x) || !is.finite(n) || n <= 0 || x < 0 || x > n) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  z <- qnorm(1 - (1 - conf.level) / 2)
  phat <- x / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- (z / denom) * sqrt((phat * (1 - phat) / n) + (z^2 / (4 * n^2)))

  c(lower = max(0, center - half), upper = min(1, center + half))
}

compute_pt50_from_repellency_curve <- function(time_min,
                                               repellency_pct,
                                               threshold_pct = 50,
                                               use_right_censoring_rule = TRUE) {
  keep <- is.finite(time_min) & is.finite(repellency_pct)
  time_min <- time_min[keep]
  repellency_pct <- repellency_pct[keep]

  if (length(time_min) == 0) {
    return(list(PT50 = NA_real_, status = "no_data"))
  }

  ord <- order(time_min)
  time_min <- time_min[ord]
  repellency_pct <- repellency_pct[ord]

  if (all(repellency_pct >= threshold_pct)) {
    if (use_right_censoring_rule) {
      return(list(PT50 = max(time_min), status = "right_censored_above_threshold"))
    } else {
      return(list(PT50 = NA_real_, status = "never_crossed_above_threshold"))
    }
  }

  first_below <- which(repellency_pct < threshold_pct)[1]

  if (first_below == 1) {
    return(list(PT50 = time_min[1], status = "already_below_threshold_at_first_time"))
  }

  t1 <- time_min[first_below - 1]
  t2 <- time_min[first_below]
  p1 <- repellency_pct[first_below - 1]
  p2 <- repellency_pct[first_below]

  if (!is.finite(p1) || !is.finite(p2) || p1 == p2) {
    return(list(PT50 = t2, status = "crossed_no_interpolation"))
  }

  pt50 <- t1 + (threshold_pct - p1) * (t2 - t1) / (p2 - p1)
  list(PT50 = pt50, status = "interpolated_crossing")
}

bootstrap_one_solution_pt50 <- function(df_solution,
                                        threshold_pct = 50,
                                        B = 2000,
                                        seed = 123,
                                        use_right_censoring_rule = TRUE) {
  set.seed(seed)

  df_solution <- df_solution %>%
    dplyr::filter(
      is.finite(Time_min),
      is.finite(N_total),
      is.finite(N_repelled),
      N_total > 0,
      N_repelled >= 0,
      N_repelled <= N_total
    )

  if (nrow(df_solution) == 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  by_time <- split(df_solution, df_solution$Time_min)
  boot_pt50 <- numeric(B)

  for (b in seq_len(B)) {
    boot_time_summary <- lapply(by_time, function(df_t) {
      idx <- sample(seq_len(nrow(df_t)), size = nrow(df_t), replace = TRUE)
      boot_df <- df_t[idx, , drop = FALSE]

      total_n <- sum(boot_df$N_total, na.rm = TRUE)
      total_r <- sum(boot_df$N_repelled, na.rm = TRUE)

      data.frame(
        Time_min = unique(df_t$Time_min)[1],
        Repellency_pct = ifelse(total_n > 0, 100 * total_r / total_n, NA_real_)
      )
    }) %>%
      dplyr::bind_rows()

    pt <- compute_pt50_from_repellency_curve(
      time_min = boot_time_summary$Time_min,
      repellency_pct = boot_time_summary$Repellency_pct,
      threshold_pct = threshold_pct,
      use_right_censoring_rule = use_right_censoring_rule
    )

    boot_pt50[b] <- pt$PT50
  }

  boot_pt50 <- boot_pt50[is.finite(boot_pt50)]
  if (length(boot_pt50) == 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  stats::quantile(boot_pt50, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
