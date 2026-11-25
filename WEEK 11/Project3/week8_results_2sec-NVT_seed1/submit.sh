#!/bin/bash
#SBATCH --job-name=h2o_300k_lammps_64
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=8
#SBATCH --time=02-00:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=gpuqs
##SBATCH --nodelist=g2,g3
##SBATCH --exclusive

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_PLACES=cores
#export OMP_PROC_BIND=close

module load lammps/20250722
srun -n 1 lmp -in 512h2o_300k.lammps.in 
