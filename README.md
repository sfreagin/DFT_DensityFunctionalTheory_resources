# Density Functional Theory and Computational Materials work

This repository keeps the inputs, analysis scripts, plots, structures, and writeups for multiple density functional theory (DFT) models and other computational materials simulations. The original source is my own coursework from MATE 250, a graduate course in Advanced Computational Materials Engineering offered at San Jose State University (Fall 2026). The calculations use Quantum ESPRESSO for density-functional theory and LAMMPS for classical and machine-learned molecular dynamics.

## Start here

| Project | What it shows | Writeup |
| --- | --- | --- |
| Germanium electronic structure | Cutoff and k-point convergence, variable-cell relaxation, and a PBE band gap with ultrasoft (SSSP) and norm-conserving (ONCV) pseudopotentials, plus a SCAN single point | [`Link`](WEEK%204/germanium_report_combined/Week4_Germanium_Combined_Report.md) |
| Ammonia vibrations and IR | PBE relaxation, Phonopy finite-displacement modes, and DFPT Born charges and dielectric tensor for the IR intensities | [`Link`](WEEK%207/nh3_report/Week7_NH3_IR_Report.md) |
| Water self-diffusion | TIP4P/2005, a SCAN-trained DeePMD potential, and SCAN AIMD, compared through the oxygen Einstein slope | [`Link`](WEEK%2011/Paper/Project3_report.md) |

The week folders are the working directories from the course, including repeated AIMD attempts under `WEEK 11/week11/week11-round*`. The table above is the path through the finished analysis.

## What else is here

- Weeks 1–3: silicon SCF and variable-cell relaxation inputs, and the homework writeups.
- Week 6: infrared and Raman homework on oscillator strengths.
- Weeks 8–10: TIP4P/2005 timestep and ensemble tests, then DeePMD water, which the diffusion project builds on. The comparison plots are in `WEEK 9/comparison_plots/`.
- Week 13: silicon conventional cell, supercell, and a defective supercell (`Si.vasp`, `Si_supercell.vasp`, `Si_supercell_defect.vasp`).

## Re-running a calculation

Each production job has a Quantum ESPRESSO `*.in` file or a LAMMPS input next to a `submit.sh`. The reports give the code versions (Quantum ESPRESSO 7.4.1 or 7.5, Phonopy 2.21.2) and the converged settings. Pseudopotentials come from the [SSSP library](https://www.materialscloud.org/discover/sssp), [ONCVPSP](http://www.quantum-simulation.org/potentials/sg15_oncv/), and SG15, not from this repository.
