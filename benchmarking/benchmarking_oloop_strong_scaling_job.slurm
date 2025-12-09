#!/bin/bash 
## Submit job to our class's allocation
#SBATCH --partition=standard
#SBATCH --account=hpc4astro_crch_fall2025

## Time requested: 0 hours, 30 minutes, 0 seconds
#SBATCH --time=0:30:00 

#SBATCH --nodes=10 
#SBATCH --ntasks-per-node=3

echo "We are using $NODES nodes, each with $SLURM_NTASKS_PER_NODE threads."

## Promise that each node will use no more than 16GB of RAM
#SBATCH --mem=16GB

## Save STDOUT and STDERR into one file (%j will expand to become the SLURM job id)
#SBATCH --output=benchmarking_oloop_strong_scaling_%j.log
## Optionally could uncomment line below to write STDERR to a separate file
##SBATCH --error=benchmarking_oloop_strong_scaling_%j.stderr  

## Specify job name, so easy to find using squeue –u
#SBATCH --job-name=benchmarking_oloop_strong_scaling

## Uncomment next two lines (by removing one of #'s in each line) and replace with your email if you want to be notifed when jobs start and stop
##SBATCH --mail-user=vql5180@psu.edu
## Ask for emails when jobs begins, ends or fails (options are ALL, NONE, BEGIN, END, FAIL)
#SBATCH --mail-type=ALL

echo "Starting job $SLURM_JOB_NAME"
echo "Job id: $SLURM_JOB_ID"
date

echo "This job was assigned the following nodes"
echo $SLURM_NODELIST

#env |grep SLURM                 # Uncomment to print slurm environment variables

echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR            # Change into directory where job was submitted from

export JULIA_PROJECT=`pwd`
export JULIA_NUM_THREADS=$SLURM_NTASKS_PER_NODE

# Repeat nodes if setting ntasks, since likely uneven distribution of tasks across nodes
#srun -l /usr/bin/hostname | sort | awk '{print $2}' > machinefile.$SLURM_JOB_ID   
# Don't repeat nodes if setting ntasks-per-node and using multiple threads per node
scontrol show hostnames | sed -e s/^/$SLURM_NTASKS_PER_NODE*/ > machinefile.$SLURM_JOB_ID   

date
echo "About to start Julia, using ssh to connect to each CPU core allocated by slurm as a separate process"

julia --machine-file machinefile.$SLURM_JOB_ID -e 'using Distributed; @everywhere using Pkg; @everywhere Pkg.activate("."); include("benchmarking_oloop_strong_scaling.jl") ' 

echo "Julia exited"
date

# Cleanup temporary file (or leave it to see what happened)
rm machinefile.$SLURM_JOB_ID