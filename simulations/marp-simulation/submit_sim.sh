#!/usr/bin/env bash
# Submit MARP simulation replicates as a SLURM array, then combine results once
# all array jobs finish successfully.
#
# Usage from the project directory:
#   bash submit_sim.sh
#
# Optional overrides:
#   MODE=point N_REPS=30 MAX_PARALLEL=10 bash submit_sim.sh
#   PROJECT_DIR=/path/to/marp-simulation MODE=ci bash submit_sim.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_DIR="${PROJECT_DIR:-$SCRIPT_DIR}"
MODE="${MODE:-both}"              # point, ci, or both
N_REPS="${N_REPS:-50}"
MAX_PARALLEL="${MAX_PARALLEL:-25}"
TIME_SIM="${TIME_SIM:-24:00:00}"
MEM_SIM="${MEM_SIM:-8G}"
CPUS_SIM="${CPUS_SIM:-1}"
TIME_COMBINE="${TIME_COMBINE:-01:00:00}"
MEM_COMBINE="${MEM_COMBINE:-8G}"
CPUS_COMBINE="${CPUS_COMBINE:-1}"
R_MODULE="${R_MODULE:-R}"
LOCAL="${LOCAL:-0}"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

load_r_module() {
  if command -v module >/dev/null 2>&1; then
    module load "$R_MODULE"
  fi
}

validate_config() {
  case "$MODE" in
    point|ci|both) ;;
    *) die "MODE must be 'point', 'ci', or 'both'." ;;
  esac

  [[ "$N_REPS" =~ ^[1-9][0-9]*$ ]] || die "N_REPS must be a positive integer."
  [[ "$MAX_PARALLEL" =~ ^[1-9][0-9]*$ ]] || die "MAX_PARALLEL must be a positive integer."
  [[ -d "$PROJECT_DIR/scripts" ]] || die "PROJECT_DIR does not look like the simulation project: $PROJECT_DIR"
}

run_sim_job() {
  local rep_id="${SLURM_ARRAY_TASK_ID:-}"
  [[ -n "$rep_id" ]] || die "SLURM_ARRAY_TASK_ID is not set for simulation job."

  echo "Running simulation replicate ${rep_id}"
  echo "MODE=${MODE}"
  echo "PROJECT_DIR=${PROJECT_DIR}"
  echo "HOSTNAME=$(hostname)"
  echo "START=$(date)"

  cd "$PROJECT_DIR"
  load_r_module
  Rscript scripts/run_sim_round.R "$rep_id" "$MODE" "$PROJECT_DIR"

  echo "END=$(date)"
}

run_combine_job() {
  echo "Combining results and creating figures"
  echo "PROJECT_DIR=${PROJECT_DIR}"
  echo "HOSTNAME=$(hostname)"
  echo "START=$(date)"

  cd "$PROJECT_DIR"
  load_r_module
  Rscript scripts/combine_results.R "$PROJECT_DIR"

  echo "END=$(date)"
}

run_local_jobs() {
  cd "$PROJECT_DIR"
  mkdir -p logs results/raw/point results/raw/ci results/summary results/figures

  echo "Running local simulation loop..."
  echo "MODE=${MODE}"
  echo "N_REPS=${N_REPS}"

  for rep_id in $(seq 1 "$N_REPS"); do
    echo "Running local replicate ${rep_id}"
    load_r_module
    Rscript scripts/run_sim_round.R "$rep_id" "$MODE" "$PROJECT_DIR"
  done

  echo "Combining local results"
  run_combine_job
}

submit_jobs() {
  if command -v sbatch >/dev/null 2>&1; then
    :
  elif [[ "$LOCAL" == "1" ]]; then
    run_local_jobs
    return
  else
    die "sbatch was not found. Run this on a SLURM login node or set LOCAL=1 to run locally."
  fi

  cd "$PROJECT_DIR"
  mkdir -p logs results/raw/point results/raw/ci results/summary results/figures

  echo "Submitting simulation array..."
  echo "PROJECT_DIR=${PROJECT_DIR}"
  echo "MODE=${MODE}"
  echo "N_REPS=${N_REPS}"
  echo "MAX_PARALLEL=${MAX_PARALLEL}"

  local sim_job_id
  sim_job_id="$(
    sbatch --parsable \
      --job-name=marp_sim \
      --output=logs/marp_%A_%a.out \
      --error=logs/marp_%A_%a.err \
      --time="$TIME_SIM" \
      --mem="$MEM_SIM" \
      --cpus-per-task="$CPUS_SIM" \
      --array="1-${N_REPS}%${MAX_PARALLEL}" \
      --export=RUN_TYPE=sim,MODE="$MODE",PROJECT_DIR="$PROJECT_DIR",R_MODULE="$R_MODULE" \
      "$0"
  )"

  echo "Submitted simulation array job: ${sim_job_id}"

  echo "Submitting combine job with dependency..."
  local combine_job_id
  combine_job_id="$(
    sbatch --parsable \
      --job-name=marp_combine \
      --output=logs/combine_%j.out \
      --error=logs/combine_%j.err \
      --time="$TIME_COMBINE" \
      --mem="$MEM_COMBINE" \
      --cpus-per-task="$CPUS_COMBINE" \
      --dependency="afterok:${sim_job_id}" \
      --export=RUN_TYPE=combine,PROJECT_DIR="$PROJECT_DIR",R_MODULE="$R_MODULE" \
      "$0"
  )"

  echo "Submitted combine job: ${combine_job_id}"
  echo "The combine job will run only after every simulation array task finishes successfully."
}

validate_config

case "${RUN_TYPE:-submit}" in
  sim)
    run_sim_job
    ;;
  combine)
    run_combine_job
    ;;
  submit)
    submit_jobs
    ;;
  *)
    die "RUN_TYPE must be unset, 'sim', or 'combine'."
    ;;
esac
