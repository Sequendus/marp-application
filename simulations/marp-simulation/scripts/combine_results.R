# scripts/combine_results.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- args[1]

summary_dir <- file.path(project_dir, "results", "summary")
figure_dir <- file.path(project_dir, "results", "figures")

dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)


# Combine point-estimate results
point_files <- list.files(
  file.path(project_dir, "results", "raw", "point"),
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(point_files) > 0) {
  
  results_point_raw <- map_dfr(point_files, read_csv, show_col_types = FALSE)
  
  results_point_clean <- results_point_raw %>%
    filter(!(scenario == "S3_Mixture_misspecified" & method == "Generating"))
  
  point_summary <- results_point_clean %>%
    filter(!is.na(estimate), !is.na(true)) %>%
    mutate(
      error = estimate - true,
      squared_error = error^2
    ) %>%
    group_by(scenario, n, method, quantity, time) %>%
    summarise(
      n_success = n(),
      bias = mean(error, na.rm = TRUE),
      mse = mean(squared_error, na.rm = TRUE),
      rmse = sqrt(mse),
      se = sd(estimate, na.rm = TRUE),
      .groups = "drop"
    )
  
  robust_summary <- results_point_clean %>%
    filter(!is.na(estimate), !is.na(true)) %>%
    mutate(
      error = estimate - true,
      abs_error = abs(error)
    ) %>%
    group_by(scenario, n, method, quantity, time) %>%
    summarise(
      median_estimate = median(estimate, na.rm = TRUE),
      median_error = median(error, na.rm = TRUE),
      mae = mean(abs_error, na.rm = TRUE),
      median_abs_error = median(abs_error, na.rm = TRUE),
      max_abs_error = max(abs_error, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(
    results_point_raw,
    file.path(summary_dir, "results_point_raw.csv")
  )
  
  write_csv(
    point_summary,
    file.path(summary_dir, "point_summary.csv")
  )
  
  write_csv(
    robust_summary,
    file.path(summary_dir, "robust_summary.csv")
  )
  

  # Figures for point estimates
  method_palette <- c(
    "Generating" = "#0072B2",
    "Best_AIC" = "#E69F00",
    "AIC_MA" = "#009E73"
  )
  
  rmse_mean_plot <- point_summary %>%
    filter(quantity == "mean") %>%
    arrange(desc(rmse)) %>%
    mutate(n = factor(n)) %>%
    ggplot(aes(x = n, y = rmse, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = method_palette) +
    scale_y_log10() +
    facet_wrap(~ scenario, scales = "free_y") +
    labs(
      title = "RMSE for mean inter-event time estimates",
      x = "Sample size",
      y = "RMSE, log scale",
      fill = "Method"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(figure_dir, "rmse_mean_interevent_time.pdf"),
    plot = rmse_mean_plot,
    width = 9,
    height = 5,
    units = "in"
  )
  
  median_error_summary <- results_point_clean %>%
    filter(quantity == "mean", !is.na(estimate), !is.na(true)) %>%
    mutate(error = estimate - true) %>%
    group_by(scenario, n, method) %>%
    summarise(
      median_error = median(error, na.rm = TRUE),
      .groups = "drop"
    )
  
  median_error_plot <- median_error_summary %>%
    mutate(n = factor(n)) %>%
    ggplot(aes(x = n, y = median_error, colour = method, group = method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = method_palette) +
    facet_wrap(~ scenario, scales = "free_y") +
    labs(
      title = "Median error for mean inter-event time estimates",
      x = "Sample size",
      y = "Median error",
      colour = "Method"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(figure_dir, "median_error_mean_interevent_time.pdf"),
    plot = median_error_plot,
    width = 9,
    height = 5,
    units = "in"
  )
  
  error_data <- results_point_clean %>%
    filter(!is.na(estimate), !is.na(true)) %>%
    mutate(
      error = estimate - true,
      n = factor(n)
    )
  
  mean_error_boxplot_zoom <- error_data %>%
    filter(quantity == "mean") %>%
    ggplot(aes(x = n, y = error, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.65,
      outlier.alpha = 0.35
    ) +
    scale_fill_manual(values = method_palette) +
    facet_wrap(~ scenario, scales = "free_y") +
    coord_cartesian(ylim = c(-150, 500)) +
    labs(
      title = "Distribution of estimation errors for mean inter-event time",
      subtitle = "Y-axis zoomed to show central error distribution",
      x = "Sample size",
      y = "Estimate - true value",
      fill = "Method"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(figure_dir, "mean_error_boxplot_zoom.pdf"),
    plot = mean_error_boxplot_zoom,
    width = 9,
    height = 5,
    units = "in"
  )
}


# Combine CI results
ci_files <- list.files(
  file.path(project_dir, "results", "raw", "ci"),
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(ci_files) > 0) {
  
  results_ci_raw <- map_dfr(ci_files, read_csv, show_col_types = FALSE)
  
  results_ci_clean <- results_ci_raw %>%
    filter(!(scenario == "S3_Mixture_misspecified" & method == "Generating"))
  
  ci_summary <- results_ci_clean %>%
    filter(!is.na(method), !is.na(quantity)) %>%
    mutate(
      error = estimate - true,
      squared_error = error^2,
      covered = lower <= true & upper >= true,
      ci_width = upper - lower,
      ci_success = !is.na(lower) & !is.na(upper)
    ) %>%
    group_by(scenario, n, method, quantity, time) %>%
    summarise(
      n_success = sum(!is.na(estimate)),
      n_ci_success = sum(ci_success, na.rm = TRUE),
      bias = mean(error, na.rm = TRUE),
      mse = mean(squared_error, na.rm = TRUE),
      rmse = sqrt(mse),
      se = sd(estimate, na.rm = TRUE),
      coverage = mean(covered[ci_success], na.rm = TRUE),
      avg_ci_width = mean(ci_width[ci_success], na.rm = TRUE),
      .groups = "drop"
    )
  
  ci_failure_summary <- results_ci_raw %>%
    group_by(scenario, n, rep) %>%
    summarise(
      ci_failed = any(ci_failed, na.rm = TRUE),
      error_message = paste(unique(na.omit(error_message)), collapse = " | "),
      .groups = "drop"
    ) %>%
    group_by(scenario, n) %>%
    summarise(
      total_reps = n(),
      ci_failures = sum(ci_failed, na.rm = TRUE),
      ci_failure_rate = mean(ci_failed, na.rm = TRUE),
      example_errors = paste(unique(na.omit(error_message)), collapse = " | "),
      .groups = "drop"
    )
  
  write_csv(
    results_ci_raw,
    file.path(summary_dir, "results_ci_raw.csv")
  )
  
  write_csv(
    ci_summary,
    file.path(summary_dir, "ci_summary.csv")
  )
  
  write_csv(
    ci_failure_summary,
    file.path(summary_dir, "ci_failure_summary.csv")
  )
  
  # CI coverage figure
  method_palette <- c(
    "Generating" = "#0072B2",
    "Best_AIC" = "#E69F00",
    "AIC_MA" = "#009E73"
  )
  
  method_shapes <- c(
    "AIC_MA" = 16,
    "Best_AIC" = 17,
    "Generating" = 15
  )
  
  pd <- position_dodge(width = 0.18)
  
  ci_coverage_plot <- ci_summary %>%
    mutate(
      n = factor(n),
      quantity_label = case_when(
        quantity == "mean" ~ "Mean",
        quantity == "probability" ~ "P(T \u2264 y)",
        quantity == "hazard" ~ paste0("Hazard at t = ", time),
        TRUE ~ quantity
      )
    ) %>%
    ggplot(aes(
      x = n,
      y = coverage,
      colour = method,
      shape = method,
      group = method
    )) +
    geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.7, position = pd) +
    geom_point(aes(size = n_ci_success), alpha = 0.9, position = pd) +
    scale_colour_manual(values = method_palette) +
    scale_shape_manual(values = method_shapes) +
    scale_size_continuous(
      breaks = c(18, 19, 20, 30, 40, 50),
      range = c(1.5, 3.2)
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0.6, 1.0),
      breaks = seq(0.6, 1.0, by = 0.1)
    ) +
    facet_grid(quantity_label ~ scenario) +
    labs(
      title = "Empirical CI coverage among successful CI computations",
      subtitle = "Dashed line indicates nominal 95% coverage; point size indicates number of successful CI computations",
      x = "Sample size",
      y = "Coverage",
      colour = "Method",
      shape = "Method",
      size = "Successful CIs"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(figure_dir, "ci_coverage_plot.pdf"),
    plot = ci_coverage_plot,
    width = 9,
    height = 6,
    units = "in"
  )
}