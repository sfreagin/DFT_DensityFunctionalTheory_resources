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
f=8

for ec in 24 28 32 36 40 44 48 52 56 60
do
er=$((ec*f))
        cat > ge_scf_test.in << EOF
&control
  calculation = 'scf',
  prefix = 'ge',
  pseudo_dir = '/fs/atipa/home/019119984/mate250/week4/germanium/pseudopotentials/' !define the path to the pseduopotential file
  outdir = './ge_out/'
  etot_conv_thr=1D-4 !this is the command to set the convergence threshold on total energy for ionic minimization
  forc_conv_thr=1D-4 !this is the command to set the convergence threshold on total force for ionic minimization
/
&system
  ibrav = 0,
  nat = 8, ntyp = 1,
  ecutwfc = ${ec}, ecutrho=${er}
  ! nbnd=20
/
&electrons
  conv_thr = 1.0d-8,
/
ATOMIC_SPECIES
 Ge 72.63 ge_pbe_v1.4.uspp.F.UPF
CELL_PARAMETERS {angstrom}
   5.6748542148543777    0.0000000000000000    0.0000000000000003
   0.0000000000000009    5.6748542148543777    0.0000000000000003
   0.0000000000000000    0.0000000000000000    5.6748542148543777
ATOMIC_POSITIONS {crystal}
 Ge  0.0000000000000000    0.0000000000000000    0.5000000000000000 
 Ge  0.2500000000000000    0.2500000000000000    0.7500000000000000 
 Ge  0.0000000000000000    0.5000000000000000    0.0000000000000000 
 Ge  0.2500000000000000    0.7500000000000000    0.2500000000000000 
 Ge  0.5000000000000000    0.0000000000000000    0.0000000000000000 
 Ge  0.7500000000000000    0.2500000000000000    0.2500000000000000 
 Ge  0.5000000000000000    0.5000000000000000    0.5000000000000000 
 Ge  0.7500000000000000    0.7500000000000000    0.7500000000000000 
K_POINTS automatic
 4 4 4 0 0 0
EOF

file=ge_redo_${ec}_${k}.out
mpirun -n 1 pw.x < ge_scf_test.in > ${file}
energy="$(grep -e '!    total energy' ${file} |awk '{print $5}')"
echo ${ec} ${energy} >> redo_ecut.dat
done

