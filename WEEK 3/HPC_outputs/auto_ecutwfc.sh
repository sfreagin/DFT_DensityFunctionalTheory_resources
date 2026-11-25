#!/bin/bash
#SBATCH --job-name=qe_%j          
#SBATCH --output=qe-%j.out
##SBATCH --error=qe-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-01:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --partition=gpuqs
#SBATCH --nodelist=g2,g3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMPI_MCA_btl_openib_allow_ib="true"

module load nvidia/nvhpc/24.11 
module load quantumespresso/7.4.1

k=8
f=4

for ec in 16 20 24 28 32
do
er=$((ec*f))
        cat > si_scf_test.in << EOF
&control
  calculation = 'scf',
  prefix = 'si',
  pseudo_dir = '/scratch/mate250-fa26/pseudos/'
  outdir = './Si_out/',
/
&system
  ibrav = 2, celldm(1) = 10.2, nat = 2, ntyp = 1,
  ecutwfc = ${ec}, ecutrho =${er}
/
&electrons
  conv_thr = 1.0d-8,
  electron_maxstep = 500
  mixing_beta = 0.3D0
 !startingpot='file'
 !startingwfc='file'
 !mixing_mode='local-TF'
 !mixing_fixed_ns=3
/
&ions
 !upscale=20
/
&cell
 !cell_dofree='xyz'
/
ATOMIC_SPECIES
 Si  28.086 Si.pbe-n-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS
 Si 0.00 0.00 0.00
 Si 0.25 0.25 0.25
K_POINTS automatic
 4 4 4 0 0 0
EOF

file=si_${ec}_${k}.out
mpirun -n 1 pw.x < si_scf_test.in > ${file}
energy="$(grep -e '!    total energy' ${file} |awk '{print $5}')"
echo ${ec} ${energy} >> ecut.dat
done

