# Classical and DeepMD water diffusion

Methods and main results for the six finished trajectories. AIMD is not included. The diffusion coefficient is the oxygen Einstein slope. In one dimension, long-time diffusion gives \(\langle x^2 \rangle = 2Dt\) [source?]. The recorded MSD sums the three Cartesian directions after the center-of-mass motion is removed, so \(\mathrm{MSD}(t) = 6Dt\) [source?]. That relation applies once the short ballistic rise is over, and the slope here is fit from 50 to 200 ps. \(D\) in m²/s is \((\mathrm{slope\ in\ Å^2/ps})/6 \times 10^{-8}\). The reported uncertainty is the sample standard deviation of the three seeds, \(s = \sqrt{\sum (D_i - \bar D)^2 / 2}\).

Both methods use 512 water molecules (1536 atoms) in the same cubic cell, \(-6.2\) to \(18.6\) Å (24.8 Å on a side, volume 15252.992 Å³). The NVT density is 1.004 g/cm³. The target temperature is 300 K.

## Classical MD (TIP4P/2005)

The production input is `512h2o_300k.lammps.in` in each `week8_results_2sec-NVT*` folder. It does not build the 512-water box. `read_restart water512_new_equilibration.restart` loads the positions, velocities, cell, and TIP4P pair style from an earlier run. `reset_timestep 0` restarts the clock, and the job then runs 100,000 steps.

That earlier run is the equilibration. `WEEK 8/week8/512h2o_300k.lammps.in` reads `water512_equilibrated.data`, minimizes, assigns velocities with seed 4928459, runs 5,000 steps (10 ps) of NVT at 300 K, and writes `water512_new_equilibration.restart`. The three production jobs start from that file. They do not repeat the 10 ps segment. The production logs sit on 300 K from the first samples.

| Setting | Value |
|---|---|
| Units | `real` |
| Ensemble | NVT, Nosé–Hoover, 300 K, damping 100 fs |
| Timestep | 2.0 fs |
| Constraints | SHAKE on bonds and angles, tolerance \(10^{-5}\) |
| Electrostatics | `lj/cut/tip4p/long`, O type 1, H type 2, M-site offset 0.1546 Å, LJ cutoff 10.0 Å; `pppm/tip4p` accuracy \(10^{-5}\) |
| LJ (O–O) | \(\varepsilon = 0.1852\) kcal/mol, \(\sigma = 3.1589\) Å |
| Geometry | \(r_{\mathrm{OH}} = 0.9572\) Å, angle 104.52°, rigid |
| Charges | O \(-1.1128\), H \(+0.5564\) (M site carries the oxygen charge in the TIP4P evaluation) |
| Production | 100,000 steps = 200 ps |
| MSD | oxygen, center-of-mass corrected (`c_msdO[4]`), written every 100 steps (0.2 ps) to `msd_NVT.dat` |
| Processors | 2 |

Seed 1 keeps the velocities stored in the restart. Seeds 2 and 3 replace them with a new 300 K draw.

| Seed | Velocity seed | Wall time | Mean \(T\) (K) | Final MSD (Å²) | Slope (Å²/ps) | \(D\) (10⁻⁹ m²/s) | \(R^2\) |
|---|---|---|---|---|---|---|---|
| 1 | restart | 31 min 23 s | 300.1 ± 7.8 | 259.42 | 1.298 | 2.163 | 0.9963 |
| 2 | 111111 | 30 min 37 s | 300.0 ± 7.6 | 276.55 | 1.416 | 2.361 | 0.9959 |
| 3 | 234561 | 21 min 51 s | 299.6 ± 7.6 | 284.59 | 1.418 | 2.364 | 0.9977 |

The temperature column is the mean and the standard deviation of the instantaneous temperature over the 200 ps. That standard deviation is the thermal fluctuation, and it scales as \(1/\sqrt{N_\mathrm{dof}}\):

\[
\sigma_T / T \approx \sqrt{2/N_\mathrm{dof}}.
\]

Rigid TIP4P water has about 6 degrees of freedom per molecule, so \(N_\mathrm{dof} \approx 6 \times 512 = 3072\) and \(\sigma_T \approx 7.7\) K, matching the 7.6–7.8 K in the table. Flexible DeepMD water has 9 Cartesian degrees of freedom per molecule, \(N_\mathrm{dof} \approx 3 \times 1536 = 4608\), and \(\sigma_T \approx 6.2\) K, matching the DeepMD table below. The two models have the same number of molecules. DeepMD fluctuates less because each molecule has more degrees of freedom. This \(\pm\) is not the uncertainty of the mean temperature, which remains 300 K, and it is not the \(\pm 0.12 \times 10^{-9}\) on \(D\). The volume does not change.

\[
D_{\mathrm{TIP4P/2005}} = 2.30 \pm 0.12 \times 10^{-9}\ \mathrm{m^2/s}.
\]

The precise three-seed values are \(\bar D = 2.296 \times 10^{-9}\) and \(s = 0.115 \times 10^{-9}\). Relative to the course experimental reference, \(2.299 \times 10^{-9}\) m²/s, the mean is 0.1% lower. The individual differences are −5.9%, +2.7%, and +2.8%.

These \(D\) and \(R^2\) values are the 50–200 ps fits already stored in each seed’s `diffusion_summary.txt` (751 points, timestep flag 2 fs). The \(\pm 0.12 \times 10^{-9}\) above is the seed-to-seed standard deviation. The “std. dev. across windows” inside those summary files is about \(0.36 \times 10^{-9}\) and is a sensitivity check on 30 ps slices, not the error bar used here. Fit figures: `msd_fit.png` in each classical seed folder.

## DeepMD

The production input is `512h2o_300k_deepmd.lammps.in` in `week10_seed1`, `week10_seed2`, and `week10_seed3`. Each job reads the same `water512_dp_equilibrated.data` and the same `graph.pb` model.

| Setting | Value |
|---|---|
| Units | `metal` |
| Model | `pair_style deepmd graph.pb`, types mapped O then H |
| Ensemble | NVT, Nosé–Hoover, 300 K, damping 0.10 ps (100 fs) |
| Timestep | 0.0005 ps = 0.5 fs |
| Constraints | none; the molecules are flexible |
| Equilibration | 20,000 steps = 10 ps, before the MSD fix is defined |
| Production | 400,000 steps = 200 ps |
| MSD | oxygen, center-of-mass corrected, written every 200 steps (0.1 ps) to `msd_O.dat` |
| Processors | 1 |

Velocities are drawn at 300 K at the start of each job. The 10 ps equilibration is already near 300 K (means 300.2, 300.2, and 299.5 K).

| Seed | Velocity seed | Equilibration | Production | Total wall time | Production mean \(T\) (K) |
|---|---|---|---|---|---|
| 1 | 23456789 | 9 min 22 s | 3 h 7 min | 3 h 16 min 39 s | 299.8 ± 6.3 |
| 2 | 12345678 | 9 min 22 s | 3 h 7 min | 3 h 16 min 50 s | 300.0 ± 6.2 |
| 3 | 11111111 | 19 min 0 s | 6 h 17 min | 6 h 36 min 16 s | 300.0 ± 6.2 |

Seed 3 completed the same number of steps on a slower node. The density stays 1.004 g/cm³ in all three logs.

The `diffusion_summary.txt` files currently in the DeepMD folders were generated with `--dt-fs 2`. That stretches the time axis by a factor of four, so those saved \(D\) values are not used here. The table below refits `msd_O.dat` from 50 to 200 ps of production time at the input timestep of 0.5 fs (1501 points).

| Seed | Final MSD (Å²) | Slope (Å²/ps) | \(D\) (10⁻⁹ m²/s) | \(R^2\) |
|---|---|---|---|---|
| 1 | 82.76 | 0.4044 | 0.674 | 0.9963 |
| 2 | 87.16 | 0.4289 | 0.715 | 0.9868 |
| 3 | 62.69 | 0.3072 | 0.512 | 0.9898 |

\[
D_{\mathrm{DeepMD}} = 0.63 \pm 0.11 \times 10^{-9}\ \mathrm{m^2/s}.
\]

The precise three-seed values are \(\bar D = 0.634 \times 10^{-9}\) and \(s = 0.107 \times 10^{-9}\). The mean is 72% below \(2.299 \times 10^{-9}\) m²/s.

## Comparison

The two production runs match in molecule count, cell, density, thermostat target, and length. The fit window is the same. The potentials and timesteps are not: TIP4P/2005 is rigid and uses a 2 fs step; DeepMD is a flexible pretrained potential and uses a 0.5 fs step. Potential energies are not compared, because the two models do not share an energy zero.

TIP4P/2005, \(2.30 \pm 0.12 \times 10^{-9}\) m²/s, agrees with the experimental reference \(2.299 \times 10^{-9}\) m²/s. DeepMD, \(0.63 \pm 0.11 \times 10^{-9}\) m²/s, is lower by a factor of about 3.6. That gap is several times either error bar. Both sets of MSD curves are linear over 50–200 ps, and both stay at 300 K and 1.004 g/cm³, so the difference is in the potential rather than in an unfinished or unstable run.

A lower diffusion coefficient is not required by the DeepMD method. The model follows the first-principles data it was trained on. TIP4P/2005 was parameterized so that liquid-water diffusion lies near experiment. A potential trained on a density functional that over-structures water, as PBE and related GGAs do, comes out slow. The factor of 3.6 found here is a property of this `graph.pb` model. It is not a general DeepMD bias, and the 0.5 fs timestep is not the cause.

Published TIP4P/2005 self-diffusion at 298 K and 1 bar is about \(2.08 \times 10^{-9}\) m²/s (Abascal and Vega, J. Chem. Phys. 123, 234505, 2005). The result here is consistent with that value within the three-seed spread.

## All six runs

Every row is a 200 ps production trajectory of 512 waters at 300 K and 1.004 g/cm³. The slope, \(D\), and \(R^2\) are the 50–200 ps oxygen MSD fit. The percent column is \((D - 2.299\times10^{-9})/2.299\times10^{-9}\). Mean rows use the three seeds of that method. The \(\pm\) on \(D\) in a mean row is the sample standard deviation. The \(\pm\) on \(T\) in a seed row is the instantaneous fluctuation.

| Run | Timestep | Wall time | Mean \(T\) (K) | Final MSD (Å²) | Slope (Å²/ps) | \(D\) (10⁻⁹ m²/s) | \(R^2\) | vs \(2.299\times10^{-9}\) |
|---|---|---|---|---|---|---|---|---|
| TIP4P/2005 seed 1 | 2 fs | 31 min 23 s | 300.1 ± 7.8 | 259.42 | 1.298 | 2.163 | 0.9963 | −5.9% |
| TIP4P/2005 seed 2 | 2 fs | 30 min 37 s | 300.0 ± 7.6 | 276.55 | 1.416 | 2.361 | 0.9959 | +2.7% |
| TIP4P/2005 seed 3 | 2 fs | 21 min 51 s | 299.6 ± 7.6 | 284.59 | 1.418 | 2.364 | 0.9977 | +2.8% |
| TIP4P/2005 mean | 2 fs | | 299.9 | 273.5 | 1.377 | 2.296 ± 0.115 | | −0.1% |
| DeepMD seed 1 | 0.5 fs | 3 h 16 min 39 s | 299.8 ± 6.3 | 82.76 | 0.4044 | 0.674 | 0.9963 | −70.7% |
| DeepMD seed 2 | 0.5 fs | 3 h 16 min 50 s | 300.0 ± 6.2 | 87.16 | 0.4289 | 0.715 | 0.9868 | −68.9% |
| DeepMD seed 3 | 0.5 fs | 6 h 36 min 16 s | 300.0 ± 6.2 | 62.69 | 0.3072 | 0.512 | 0.9898 | −77.7% |
| DeepMD mean | 0.5 fs | | 299.9 | 77.5 | 0.380 | 0.634 ± 0.107 | | −72.4% |
