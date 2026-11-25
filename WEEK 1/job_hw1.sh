#!/bin/bash
#SBATCH --job-name=homework1-
#SBATCH --output=hw1-%j.out
#SBATCH --error=hw1-%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00   
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=nsfql
#SBATCH --gres=gpu:1
#SBATCH --exclusive

echo "Job ID: $SLURM_JOB_ID"
echo "Number of tasks: $SLURM_NTASKS"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module load nvhpc/24.11
module load quantumespresso/7.5.1

mpirun -n $SLURM_NTASKS pw.x < si.in > si.out