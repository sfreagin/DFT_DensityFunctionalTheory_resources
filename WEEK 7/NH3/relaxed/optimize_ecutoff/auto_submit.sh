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

for ec in 60 62 64 66 68 70 72 74 76 78 80
do
er=$((ec*f))
        cat > nh3_scf.in << EOF
&control
  calculation = 'scf',
  prefix = 'nh3',
  pseudo_dir = '/scratch/mate250-fa26/pseudos/SSSP_1.3.0_PBE_efficiency/'
  outdir = './nh3_out/'
  forc_conv_thr=1.0D-5
/
&system
  ibrav = 0,
  nat = 4, ntyp = 2,
  ecutwfc = ${ec}, ecutrho=${er}
/
&electrons
  conv_thr = 1.0d-8,
/
&ions
/
ATOMIC_SPECIES
 H  1.008 H.pbe-rrkjus_psl.1.0.0.UPF
 N  14.007 N.pbe-n-radius_5.UPF 
CELL_PARAMETERS {angstrom}
  15.0   0.000000000    0.000000000000000
  0.0000000    15.0    0.0000000000000000
   0.0000000000000000    0.0000000000000000    15.0
ATOMIC_POSITIONS {crystal} 
 N 0.5 0.5   0.5
 H 0.5 0.437486667 0.47456
 H 0.55414 0.53126 0.47456
 H 0.44586 0.53126 0.47456
 ! NH3 values from NIST: https://cccbdb.nist.gov/expgeom2x.asp?casno=7664417&charge=0
K_POINTS automatic
 1 1 1 0 0 0
EOF

file=nh3_${ec}_${k}.out
mpirun -n 1 pw.x < nh3_scf.in > ${file}
energy="$(grep -e '!    total energy' ${file} |awk '{print $5}')"
echo ${ec} ${energy} >> ecut.dat
done

