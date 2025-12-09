#!/bin/bash 
## Submit job to our class's allocation
#SBATCH --partition=standard
#SBATCH --account=hpc4astro_crch_fall2025

## Time requested: 0 hours, 30 minutes, 0 seconds
#SBATCH --time=0:30:00 

## Ask for one core on one node for each element of the task array
#SBATCH --nodes=1 
#SBATCH --ntasks=6 

## Promise that each CPU core will use no more than 8GB of RAM
#SBATCH --mem-per-cpu=8GB

## Save STDOUT and STDERR into one file (%A_%a will expand to become the SLURM job id of the first job in the array and then the task id within the array)
#SBATCH --output=benchmarking_iloop_weak_scaling_%A_%a.log
## Optionally could uncomment line below to write STDERR to a separate file
##SBATCH --error=benchmarking_iloop_weak_scaling_%A_%a.stderr  

## Specificy job name, so easy to find using squeue 
#SBATCH --job-name=benchmarking_iloop_weak_scaling

## Submit a job array with 10 tasks, but limit to 5 tasks running at once
#SBATCH --array=0-9%5

## Uncomment next two lines (by removing one of #'s in each line) and replace with your email if you want to be notifed when jobs start and stop
##SBATCH --mail-user=vql5180@psu.edu
## Ask for emails when jobs begins, ends or fails (options are ALL, NONE, BEGIN, END, FAIL)
#SBATCH --mail-type=ALL

export JULIA_PROJECT=`pwd`
export JULIA_NUM_THREADS=$SLURM_NTASKS

echo "Starting job $SLURM_JOB_NAME"
echo "Job id: $SLURM_JOB_ID"
date

PARAMS=(16 32 64 128 256 512 1024 2048 4096 8192)
STEP_LIMIT=${PARAMS[$SLURM_ARRAY_TASK_ID]}

echo "Launching Julia with STEP_LIMIT = $STEP_LIMIT"

echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR            # Change into directory where job was submitted from

echo "About to start Julia"
julia --project=. -- benchmarking_weak_scaling.jl $STEP_LIMIT 
echo "Julia exited"
date