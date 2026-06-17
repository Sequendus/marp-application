# scripts/run_sim_round.R

args <- commandArgs(trailingOnly = TRUE)

rep_id <- as.integer(args[1])
mode <- args[2]
project_dir <- args[3]

source(file.path(project_dir, "scripts", "simulation_functions.R"))

scenario_ids <- seq_along(scenarios)

base_seed <- 12345

run_point_round <- function(rep_id, project_dir) {
  
  mode_offset <- 0
  
  out <- purrr::map_dfr(
    scenario_ids,
    function(scenario_id) {
      
      scenario <- scenarios[[scenario_id]]
      
      purrr::map_dfr(
        n_vec,
        function(n) {
          
          set.seed(base_seed + mode_offset + rep_id * 1000 + scenario_id * 100 + n)
          
          message("Point simulation: scenario = ", scenario_id,
                  ", n = ", n,
                  ", rep = ", rep_id)
          
          run_one_sim_point(
            scenario = scenario,
            n = n,
            rep_id = rep_id,
            t_vec = t_vec,
            y = y,
            m = m
          )
        }
      )
    }
    
  )
  
  out_file <- file.path(
    project_dir,
    "results", "raw", "point",
    paste0("point_rep", rep_id, ".csv")
  )
  
  readr::write_csv(out, out_file)
}

run_ci_round <- function(rep_id, project_dir) {
  
  mode_offset <- 500000
  
  out <- purrr::map_dfr(
    scenario_ids,
    function(scenario_id) {
      scenario <- scenarios[[scenario_id]]
      
      purrr::map_dfr(
        n_vec_ci,
        function(n) {
          
          set.seed(base_seed + mode_offset + rep_id * 1000 + scenario_id * 100 + n)
          
          message("CI simulation: scenario = ", scenario_id,
                  ", n = ", n,
                  ", rep = ", rep_id)
          
          run_one_sim_ci(
            scenario = scenario,
            n = n,
            rep_id = rep_id,
            t_vec = t_vec_ci,
            y = y,
            m = m,
            B = B,
            BB = BB,
            alpha = alpha
          )
        }
      )
    }
    
  )
  
  out_file <- file.path(
    project_dir,
    "results", "raw", "ci",
    paste0("ci_rep", rep_id, ".csv")
  )
  
  readr::write_csv(out, out_file)
}

if (mode == "point") {
  run_point_round(rep_id, project_dir)
} else if (mode == "ci") {
  run_ci_round(rep_id, project_dir)
} else if (mode == "both") {
  run_point_round(rep_id, project_dir)
  run_ci_round(rep_id, project_dir)
} else {
  stop("mode must be 'point', 'ci', or 'both'")
}
