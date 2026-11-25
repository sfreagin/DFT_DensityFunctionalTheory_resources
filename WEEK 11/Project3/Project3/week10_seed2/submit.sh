#!/bin/bash
#SBATCH --job-name=h2o_300k_lammps_64
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=8
#SBATCH --time=02-00:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=gpuqs
##SBATCH --nodelist=condo10,condo11
##SBATCH --exclusive

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#export OMP_PLACES=cores
#export OMP_PROC_BIND=close

module load lammps/deepmd-kit
srun -n 1 lmp -in 512h2o_300k_deepmd.lammps.in
