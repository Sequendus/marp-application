# scripts/simulation_functions.R

suppressPackageStartupMessages({
  library(marp)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
})

# -----------------------------
# Simulation settings
# -----------------------------

# Main point-estimate simulation
n_vec <- c(20, 50, 100)
t_vec <- c(150, 200)

# CI simulation
n_vec_ci <- c(50, 100)
t_vec_ci <- c(150, 200)

# Bootstrap settings
B <- 99
BB <- 99
alpha <- 0.05

# MARP settings
m <- 10
y <- 304

# MARP model index:
# 1 = Poisson
# 2 = Gamma
# 3 = Log-logistic
# 4 = Weibull
# 5 = Log-normal
# 6 = BPT


# Simulation scenarios
sim_gamma <- list(
  name = "S1_Gamma_well_specified",
  which_model = 2,
  
  rfun = function(n) {
    rgamma(n, shape = 3, rate = 0.01)
  },
  
  true_fun = function(t_vec, y) {
    shape <- 3
    rate <- 0.01
    
    mu_true <- shape / rate
    pr_true <- pgamma(y, shape = shape, rate = rate)
    
    haz_true <- dgamma(t_vec, shape = shape, rate = rate) /
      (1 - pgamma(t_vec, shape = shape, rate = rate))
    
    list(mu = mu_true, pr = pr_true, haz = haz_true)
  }
)

sim_weibull_challenging <- list(
  name = "S2_Weibull_challenging",
  which_model = 4,
  
  rfun = function(n) {
    rweibull(n, shape = 0.7, scale = 300)
  },
  
  true_fun = function(t_vec, y) {
    shape <- 0.7
    scale <- 300
    
    mu_true <- scale * gamma(1 + 1 / shape)
    pr_true <- pweibull(y, shape = shape, scale = scale)
    
    haz_true <- dweibull(t_vec, shape = shape, scale = scale) /
      (1 - pweibull(t_vec, shape = shape, scale = scale))
    
    list(mu = mu_true, pr = pr_true, haz = haz_true)
  }
)

sim_mixture <- list(
  name = "S3_Mixture_misspecified",
  
  # There is no true MARP generating model for the mixture.
  which_model = 2,
  
  rfun = function(n) {
    z <- rbinom(n, size = 1, prob = 0.7)
    
    x_gamma <- rgamma(n, shape = 3, rate = 0.01)
    x_weib <- rweibull(n, shape = 0.7, scale = 500)
    
    ifelse(z == 1, x_gamma, x_weib)
  },
  
  true_fun = function(t_vec, y) {
    w <- 0.7
    
    mu_gamma <- 3 / 0.01
    mu_weib <- 500 * gamma(1 + 1 / 0.7)
    mu_true <- w * mu_gamma + (1 - w) * mu_weib
    
    pr_true <- w * pgamma(y, shape = 3, rate = 0.01) +
      (1 - w) * pweibull(y, shape = 0.7, scale = 500)
    
    f_mix <- w * dgamma(t_vec, shape = 3, rate = 0.01) +
      (1 - w) * dweibull(t_vec, shape = 0.7, scale = 500)
    
    F_mix <- w * pgamma(t_vec, shape = 3, rate = 0.01) +
      (1 - w) * pweibull(t_vec, shape = 0.7, scale = 500)
    
    haz_true <- f_mix / (1 - F_mix)
    
    list(mu = mu_true, pr = pr_true, haz = haz_true)
  }
)

scenarios <- list(
  sim_gamma,
  sim_weibull_challenging,
  sim_mixture
)


# Helper functions: point estimates
extract_point_results <- function(
    fit,
    scenario_name,
    n,
    rep_id,
    t_vec,
    true_vals
) {
  
  mean_rows <- tibble(
    scenario = scenario_name,
    n = n,
    rep = rep_id,
    method = c("Generating", "Best_AIC", "AIC_MA"),
    quantity = "mean",
    time = NA_real_,
    estimate = c(
      as.numeric(fit$mu_gen),
      as.numeric(fit$mu_best),
      as.numeric(fit$mu_aic)
    ),
    true = true_vals$mu,
    lower = NA_real_,
    upper = NA_real_,
    ci_failed = NA,
    error_message = NA_character_
  )
  
  pr_rows <- tibble(
    scenario = scenario_name,
    n = n,
    rep = rep_id,
    method = c("Generating", "Best_AIC", "AIC_MA"),
    quantity = "probability",
    time = NA_real_,
    estimate = c(
      as.numeric(fit$pr_gen),
      as.numeric(fit$pr_best),
      as.numeric(fit$pr_aic)
    ),
    true = true_vals$pr,
    lower = NA_real_,
    upper = NA_real_,
    ci_failed = NA,
    error_message = NA_character_
  )
  
  # MARP returns hazard estimates on the log-hazard scale.
  # Transform to hazard scale with exp().
  haz_rows <- bind_rows(
    tibble(
      scenario = scenario_name,
      n = n,
      rep = rep_id,
      method = "Generating",
      quantity = "hazard",
      time = t_vec,
      estimate = exp(as.numeric(fit$haz_gen)),
      true = true_vals$haz,
      lower = NA_real_,
      upper = NA_real_,
      ci_failed = NA,
      error_message = NA_character_
    ),
    tibble(
      scenario = scenario_name,
      n = n,
      rep = rep_id,
      method = "Best_AIC",
      quantity = "hazard",
      time = t_vec,
      estimate = exp(as.numeric(fit$haz_best)),
      true = true_vals$haz,
      lower = NA_real_,
      upper = NA_real_,
      ci_failed = NA,
      error_message = NA_character_
    ),
    tibble(
      scenario = scenario_name,
      n = n,
      rep = rep_id,
      method = "AIC_MA",
      quantity = "hazard",
      time = t_vec,
      estimate = exp(as.numeric(fit$haz_aic)),
      true = true_vals$haz,
      lower = NA_real_,
      upper = NA_real_,
      ci_failed = NA,
      error_message = NA_character_
    )
  )
  
  bind_rows(mean_rows, pr_rows, haz_rows)
}

run_one_sim_point <- function(
    scenario,
    n,
    rep_id,
    t_vec,
    y,
    m
) {
  
  dat <- scenario$rfun(n)
  true_vals <- scenario$true_fun(t_vec = t_vec, y = y)
  
  fit <- tryCatch(
    marp::marp(
      data = dat,
      t = t_vec,
      m = m,
      y = y,
      which.model = scenario$which_model
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(
      tibble(
        scenario = scenario$name,
        n = n,
        rep = rep_id,
        method = NA_character_,
        quantity = NA_character_,
        time = NA_real_,
        estimate = NA_real_,
        true = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        ci_failed = NA,
        error_message = fit$message
      )
    )
  }
  
  extract_point_results(
    fit = fit,
    scenario_name = scenario$name,
    n = n,
    rep_id = rep_id,
    t_vec = t_vec,
    true_vals = true_vals
  )
}

# Helper functions: confidence intervals
safe_get <- function(x, name) {
  if (!is.null(x) && name %in% names(x)) {
    x[[name]]
  } else {
    NA
  }
}

make_rows_from_ci <- function(
    ci,
    scenario_name,
    n,
    rep_id,
    t_vec,
    true_vals
) {
  
  out <- ci$out
  student <- ci$student_CI
  
  mean_rows <- tibble(
    scenario = scenario_name,
    n = n,
    rep = rep_id,
    method = c("Generating", "Best_AIC", "AIC_MA"),
    quantity = "mean",
    time = NA_real_,
    estimate = c(
      as.numeric(safe_get(out, "mu_gen"))[1],
      as.numeric(safe_get(out, "mu_best"))[1],
      as.numeric(safe_get(out, "mu_aic"))[1]
    ),
    true = true_vals$mu,
    lower = c(
      as.numeric(safe_get(student, "mu_lower_gen"))[1],
      as.numeric(safe_get(student, "mu_lower_best"))[1],
      as.numeric(safe_get(student, "mu_lower_ma"))[1]
    ),
    upper = c(
      as.numeric(safe_get(student, "mu_upper_gen"))[1],
      as.numeric(safe_get(student, "mu_upper_best"))[1],
      as.numeric(safe_get(student, "mu_upper_ma"))[1]
    )
  )
  
  pr_rows <- tibble(
    scenario = scenario_name,
    n = n,
    rep = rep_id,
    method = c("Generating", "Best_AIC", "AIC_MA"),
    quantity = "probability",
    time = NA_real_,
    estimate = c(
      as.numeric(safe_get(out, "pr_gen"))[1],
      as.numeric(safe_get(out, "pr_best"))[1],
      as.numeric(safe_get(out, "pr_aic"))[1]
    ),
    true = true_vals$pr,
    lower = c(
      as.numeric(safe_get(student, "pr_lower_gen"))[1],
      as.numeric(safe_get(student, "pr_lower_best"))[1],
      as.numeric(safe_get(student, "pr_lower_ma"))[1]
    ),
    upper = c(
      as.numeric(safe_get(student, "pr_upper_gen"))[1],
      as.numeric(safe_get(student, "pr_upper_best"))[1],
      as.numeric(safe_get(student, "pr_upper_ma"))[1]
    )
  )
  
  # MARP returns hazard estimates on the log-hazard scale.
  haz_rows <- bind_rows(
    tibble(
      scenario = scenario_name,
      n = n,
      rep = rep_id,
      method = "Generating",
      quantity = "hazard",
      time = t_vec,
      estimate = exp(as.numeric(safe_get(out, "haz_gen")))[seq_along(t_vec)],
      true = true_vals$haz,
      lower = exp(as.numeric(safe_get(student, "haz_lower_gen")))[seq_along(t_vec)],
      upper = exp(as.numeric(safe_get(student, "haz_upper_gen")))[seq_along(t_vec)]
    ),
    tibble(
      scenario = scenario_name,
      n = n,
      rep = rep_id,
      method = "Best_AIC",
      quantity = "hazard",
      time = t_vec,
      estimate = exp(as.numeric(safe_get(out, "haz_best")))[seq_along(t_vec)],
      true = true_vals$haz,
      lower = exp(as.numeric(safe_get(student, "haz_lower_best")))[seq_along(t_vec)],
      upper = exp(as.numeric(safe_get(student, "haz_upper_best")))[seq_along(t_vec)]
    ),
    tibble(
      scenario = scenario_name,
      n = n,
      rep = rep_id,
      method = "AIC_MA",
      quantity = "hazard",
      time = t_vec,
      estimate = exp(as.numeric(safe_get(out, "haz_aic")))[seq_along(t_vec)],
      true = true_vals$haz,
      lower = exp(as.numeric(safe_get(student, "haz_lower_ma")))[seq_along(t_vec)],
      upper = exp(as.numeric(safe_get(student, "haz_upper_ma")))[seq_along(t_vec)]
    )
  )
  
  bind_rows(mean_rows, pr_rows, haz_rows) %>%
    mutate(
      ci_failed = FALSE,
      error_message = NA_character_
    )
}

run_one_sim_ci <- function(
    scenario,
    n,
    rep_id,
    t_vec,
    y,
    m,
    B,
    BB,
    alpha
) {
  
  dat <- scenario$rfun(n)
  true_vals <- scenario$true_fun(t_vec = t_vec, y = y)
  
  ci <- tryCatch(
    marp::marp_confint(
      data = dat,
      m = m,
      t = t_vec,
      B = B,
      BB = BB,
      alpha = alpha,
      y = y,
      which.model = scenario$which_model
    ),
    error = function(e) e
  )
  
  if (inherits(ci, "error")) {
    return(
      tibble(
        scenario = scenario$name,
        n = n,
        rep = rep_id,
        method = NA_character_,
        quantity = NA_character_,
        time = NA_real_,
        estimate = NA_real_,
        true = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        ci_failed = TRUE,
        error_message = ci$message
      )
    )
  }
  
  make_rows_from_ci(
    ci = ci,
    scenario_name = scenario$name,
    n = n,
    rep_id = rep_id,
    t_vec = t_vec,
    true_vals = true_vals
  )
}