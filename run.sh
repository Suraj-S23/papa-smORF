#!/bin/bash

#SBATCH --partition=mlhiwidlc_gpu-rtx2080
#SBATCH --job-name=smorf-data
#SBATCH --output=run_logs/proc2%x-%A-output.out
#SBATCH --error=run_logs/proc2%x-%A-error.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssurajj08@gmail.com
#SBATCH --mem 264GB
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00	

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
export OUTPUT_DIR="data/large_processed_2"

# Start progress monitoring in background using the dedicated script
echo "Monitoring Process"
python3 src/tools/monitor.py "${OUTPUT_DIR}" &
MONITOR_PID=$!

# Record start time
start=`date +%s`

# Run the pipeline
make process-large

# Stop progress monitoring
kill $MONITOR_PID

# Calculate runtime
end=`date +%s`
runtime=$((end-start))

echo "Job execution complete."
echo "Runtime: $runtime seconds"
echo "Finished at $(date)"
