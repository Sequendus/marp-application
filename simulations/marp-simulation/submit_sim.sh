#!/bin/bash
# ============================================================
# submit_sim.sh
# One-file workflow:

# 1. Submit simulation array jobs

# 2. After all simulation jobs finish successfully

# 3. Run combine_results.R once to create summaries and figures

# Run this file using:
# bash submit_sim.sh
# ============================================================

PROJECT_DIR="$HOME/Desktop/marp-application/simulations/marp-simulation"

MODE="both"     # Options: "point", "ci", or "both"
N_REPS=50
MAX_PARALLEL=25  

if [ "$RUN_TYPE" = "sim" ]; then

echo "Running simulation replicate ${SLURM_ARRAY_TASK_ID}"
echo "MODE=${MODE}"
echo "PROJECT_DIR=${PROJECT_DIR}"

cd "$PROJECT_DIR" || exit 1

module load R

Rscript scripts/run_sim_round.R "$SLURM_ARRAY_TASK_ID" "$MODE" "$PROJECT_DIR"

exit 0
fi

if [ "$RUN_TYPE" = "combine" ]; then

echo "Combining results and creating figures"
echo "PROJECT_DIR=${PROJECT_DIR}"

cd "$PROJECT_DIR" || exit 1

module load R

Rscript scripts/combine_results.R "$PROJECT_DIR"

exit 0
fi

cd "$PROJECT_DIR" || exit 1

echo "Submitting simulation array..."
echo "MODE=${MODE}"
echo "N_REPS=${N_REPS}"
echo "MAX_PARALLEL=${MAX_PARALLEL}"

SIM_JOB_ID=$(sbatch --parsable 
--job-name=marp_sim 
--output=logs/marp_%A_%a.out 
--error=logs/marp_%A_%a.err 
--time=24:00:00 
--mem=8G 
--cpus-per-task=1 
--array=1-${N_REPS}%${MAX_PARALLEL} 
--export=ALL,RUN_TYPE=sim,MODE="${MODE}",PROJECT_DIR="${PROJECT_DIR}" 
"$0")

echo "Submitted simulation array job: ${SIM_JOB_ID}"

echo "Submitting combine job with dependency..."
COMBINE_JOB_ID=$(sbatch --parsable 
--job-name=marp_combine 
--output=logs/combine_%j.out 
--error=logs/combine_%j.err 
--time=01:00:00 
--mem=8G 
--cpus-per-task=1 
--dependency=afterok:${SIM_JOB_ID} 
--export=ALL,RUN_TYPE=combine,PROJECT_DIR="${PROJECT_DIR}" 
"$0")

echo "Submitted combine job: ${COMBINE_JOB_ID}"
echo "The combine job will run only after all simulation array jobs finish successfully."

