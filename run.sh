#!/bin/bash

#SBATCH --partition=mlhiwidlc_gpu-rtx2080
#SBATCH --job-name=smorf-data
#SBATCH --output=run_logs/%x-%A-output.out
#SBATCH --error=run_logs/%x-%A-error.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssurajj08@gmail.com
#SBATCH --mem 128GB
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00	

echo "Working directory: $PWD"
echo "Started at $(date)"
echo "Running job $SLURM_JOB_NAME using $SLURM_JOB_CPUS_PER_NODE CPUs with job ID $SLURM_JOB_ID"

# Activate conda environment
source ~/miniconda3/bin/activate
conda activate smorfenv

# Set PYTHONPATH
export PYTHONPATH=$PWD:$PYTHONPATH

# Set input/output paths
export FASTA_DIR="data/large_fasta/fna_files"
export GFF_DIR="data/large_gff/gff_files"
export OUTPUT_DIR="data/large_processed"

# Function to monitor progress
monitor_progress() {
    while true; do
        if [ -f "${OUTPUT_DIR}/progress.json" ]; then
            echo "=== Progress Update $(date) ==="
            python -c "
import json
with open('${OUTPUT_DIR}/progress.json') as f:
    p = json.load(f)
    print(f'Progress: {p[\"completion_percentage\"]:.1f}% complete')
    print(f'Processed: {p[\"processed_files\"]}/{p[\"total_files\"]} files')
    print(f'Batch: {p[\"current_batch\"]}/{p[\"total_batches\"]}')
    print(f'Total ORFs: {p[\"total_orfs\"]:,}')
    print(f'Memory Usage: {p[\"memory_usage\"]:.2f} GB')
    print(f'CPU Usage: {p[\"cpu_usage\"]:.1f}%')
"
        fi
        sleep 60  # Update every minute
    done
}

# Start progress monitoring in background
monitor_progress &
MONITOR_PID=$!

# Record start time
start=`date +%s`

# Run the pipeline
make run-all

# Stop progress monitoring
kill $MONITOR_PID

# Calculate runtime
end=`date +%s`
runtime=$((end-start))

echo "Job execution complete."
echo "Runtime: $runtime seconds"
echo "Finished at $(date)"
