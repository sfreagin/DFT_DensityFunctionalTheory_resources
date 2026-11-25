#!/bin/bash
#SBATCH --job-name=qe_%j
#SBATCH --output=qe-%j.out
#SBATCH --error=qe-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=gpuqs
#SBATCH --nodelist=g2,g3
##SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module load nvidia/nvhpc/24.11
module load quantumespresso/7.4.1

#mpirun -n $SLURM_NTASKS pw.x < nh3_relaxed.scf.in > nh3_relaxed.scf.out
mpirun -n $SLURM_NTASKS ph.x < nh3.ph.in > nh3.ph.out
