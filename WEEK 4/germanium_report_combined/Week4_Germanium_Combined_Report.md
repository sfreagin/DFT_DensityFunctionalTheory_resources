# Week 4 Project: Electronic structure of germanium with ultrasoft and norm-conserving pseudopotentials

**MATE 250**  
**Stephen Reagin**

*AI disclaimer: I did not use AI tools to design or perform any experimental runs in the HPC server, in part because access is gated through a secure shell, and in part because I am still learning the basics of SLURM scripting on my own. However, I did use Cursor AI to help grep values from output files, design plots, and organize the contents of this report. I also took advantage of the AI Research Assistant offered by the SJSU library to find references, using prompts such as: “I want to know what pseudo potentials and other parameters have been used to calculate the ground state and band structure of Germanium crystals in a DFT context.” Any misstatements of fact, whether hallucinated by AI or not, are ultimately my responsibility.*

## Summary

The band gap of diamond-structure germanium was computed with Quantum ESPRESSO using two pseudopotentials: the SSSP ultrasoft potential assigned for the project and a norm-conserving ONCV potential. Both were carried through cutoff and k-point convergence, a PBE variable-cell relaxation, and a final SCF. The ONCV set was also used for a SCAN meta-GGA calculation on the relaxed structure.

| Quantity | Ultrasoft, PBE | ONCV, PBE | ONCV, SCAN |
| --- | ---: | ---: | ---: |
| Converged `ecutwfc` | 56 Ry | 80 Ry | 80 Ry |
| Relaxed lattice parameter $a$ | 5.758803 &Aring; | 5.762454 &Aring; | (ONCV PBE cell used) |
| Band gap | 0.00 eV | 0.00 eV | **-0.097 eV** |

The two PBE relaxations agree to 0.0037 &Aring; (0.06%), and both are 1.8% larger than the experimental lattice parameter. Neither PBE calculation opens a gap. The SCAN calculation gives a negative gap of -0.097 eV: its lowest unoccupied level lies below its highest occupied level, so the valence and conduction bands overlap. None of the calculations reproduces the experimental gap of 0.66 to 0.74 eV.

## 1. Why two pseudopotentials

The project specifies the SSSP efficiency library, whose germanium entry is the ultrasoft file `ge_pbe_v1.4.uspp.F.UPF`. It is soft, so it converges at a low cutoff, and it gives an efficient PBE workflow. It does not support the functionals the project offers for the band gap.

- **SCAN.** A SCAN run with a PAW potential stopped immediately with `Meta-GGA not implemented with USPP/PAW`. In Quantum ESPRESSO, SCAN requires a norm-conserving pseudopotential.
- **HSE.** With the ultrasoft file, the hybrid-functional run warned that `ecutfock` is ignored for US/PAW potentials and that a nonlinear core correction is not consistent with hybrid exchange.

A norm-conserving ONCV potential removes both restrictions. Running the full PBE workflow again with that potential also tests how much the structural and electronic results depend on the pseudopotential, which is part of what a gap calculation for germanium has to establish.

## 2. Calculation settings

| Setting | Ultrasoft | ONCV (PBE) | ONCV (SCAN) |
| --- | --- | --- | --- |
| Exchange-correlation functional | PBE | PBE | SCAN (`input_dft = 'scan'`) |
| Code | QE 7.4.1 | QE 7.4.1 | QE 7.5 |
| Pseudopotential | `ge_pbe_v1.4.uspp.F.UPF` (SSSP 1.3.0 efficiency) | `Ge.upf` (ONCVPSP 3.3.0) | `Ge_ONCV_PBE-1.2.upf` (SG15) |
| Type | Ultrasoft, nonlinear core correction | Norm-conserving, nonlinear core correction | Norm-conserving, no core correction |
| Valence electrons | 14 (3d, 4s, 4p) | 14 (3d, 4s, 4p) | 14 (3d, 4s, 4p) |
| `ecutwfc` | 56 Ry | 80 Ry | 80 Ry |
| `ecutrho` | 448 Ry ($8 \times$) | 640 Ry ($8 \times$) | 320 Ry ($4 \times$, default) |
| k-point mesh, relaxation / final SCF | $8 \times 8 \times 8$ | $14 \times 14 \times 14$ | $8 \times 8 \times 8$ (SCF only) |
| `conv_thr` | `1.0d-8` Ry | `1.0d-8` Ry | `1.0d-8` Ry |
| `etot_conv_thr` | `1.0d-4` Ry | `1.0d-4` Ry | -- |
| `forc_conv_thr` | `1.0d-4` Ry/Bohr | `1.0d-4` Ry/Bohr | -- |
| Cell-pressure threshold | 0.5 kbar | 0.5 kbar | -- |
| Relaxation | `vc-relax`, `cell_dofree = 'xyz'` | `vc-relax`, `cell_dofree = 'xyz'` | Single point, $a$ = 5.762454 &Aring; (`ge.scf_final.out`) |

All calculations use the 8-atom conventional cube of diamond Ge. The starting structure is from the Materials Project, with $a$ = 5.674854 &Aring;. Both relaxations met the energy, force, and pressure thresholds in 4 BFGS steps.

An ultrasoft potential needs `ecutrho` $= 8 \times$ `ecutwfc` because the augmentation charges are much harder than the wave functions. For a norm-conserving potential $4 \times$ is sufficient. The ONCV PBE runs kept the $8 \times$ factor. That is a finer density grid than necessary, applied consistently to every cutoff test, the k-point test, and the relaxation. The SCAN run used the default $4 \times$.

Total energies from the two pseudopotentials are not comparable with each other, and neither are absolute Kohn-Sham eigenvalues. The reported highest occupied levels, 7.78 eV (ultrasoft PBE), 8.79 eV (ONCV PBE), and 7.64 eV (ONCV SCAN), differ because a periodic calculation has no absolute energy reference. Only energy differences within one calculation are meaningful.

## 3. Convergence

Cutoff tests were static SCF calculations on the starting cell with `ecutrho` $= 8 \times$ `ecutwfc`. The ultrasoft test used a $4 \times 4 \times 4$ mesh; the ONCV test used $8 \times 8 \times 8$. The k-point tests used the selected cutoff for each potential. Figure 1 shows both.

![Panel (a): energy per atom versus wave-function cutoff for the ultrasoft and ONCV germanium pseudopotentials, each relative to its highest tested cutoff, on a logarithmic scale. Panel (b): energy per atom versus Monkhorst-Pack mesh for both pseudopotentials, relative to a 14 x 14 x 14 mesh.](ge_convergence_us_vs_oncv.png)

**Figure 1.** (a) Energy above the highest tested cutoff (60 Ry ultrasoft, 92 Ry ONCV), in meV/atom. Dashed lines mark the selected cutoffs. (b) Energy above the $14 \times 14 \times 14$ result at the selected cutoff. The dotted line is 1 meV/atom.

| `ecutwfc` (Ry) | Ultrasoft, above 60 Ry (meV/atom) | ONCV, above 92 Ry (meV/atom) |
| ---: | ---: | ---: |
| 36 | 10.7 | 60472 |
| 44 | 4.97 | 17297 |
| 48 | 2.01 | 8167 |
| 52 | 0.27 | 3500 |
| **56** | **0.010** | 1346 |
| 60 | 0 | 455 |
| 68 | -- | 30.8 |
| 72 | -- | 5.29 |
| 76 | -- | 0.60 |
| **80** | -- | **0.14** |
| 88 | -- | 0.042 |

I treated a cutoff as converged once a further increase changed the energy by less than 1 meV/atom. The ultrasoft potential crosses that line between 48 and 52 Ry, and 56 Ry leaves a residual of 0.010 meV/atom. The ONCV potential crosses it between 72 and 76 Ry, and 80 Ry leaves 0.14 meV/atom. The ONCV energy is still 455 meV/atom from convergence at 60 Ry, where the ultrasoft energy has already stopped changing. That is the expected cost of norm conservation: without augmentation charges, the pseudo-wave functions of the localized Ge 3d states must be represented directly in plane waves.

| Mesh | Ultrasoft, 56 Ry (meV/atom) | ONCV, 80 Ry (meV/atom) |
| --- | ---: | ---: |
| $2 \times 2 \times 2$ | 144.5 | -- |
| $4 \times 4 \times 4$ | 6.96 | 6.98 |
| $6 \times 6 \times 6$ | 0.85 | -- |
| $8 \times 8 \times 8$ | 0.174 | 0.175 |
| $10 \times 10 \times 10$ | 0.047 | 0.048 |
| $12 \times 12 \times 12$ | 0.013 | 0.013 |

The k-point errors from the two pseudopotentials agree to within 0.02 meV/atom at every common mesh. k-point convergence is controlled by the band structure, which both potentials describe the same way, so this agreement is a consistency check on the two setups. An $8 \times 8 \times 8$ mesh is converged to 0.17 meV/atom. The ultrasoft relaxation used that mesh, and the ONCV relaxation used $14 \times 14 \times 14$.

## 4. Optimized lattice parameters

Both relaxed cells remain cubic, and the atoms stay on the diamond sites to within $3 \times 10^{-6}$ (ultrasoft) and $9 \times 10^{-7}$ (ONCV) in fractional coordinates.

| Quantity | Ultrasoft, PBE | ONCV, PBE |
| --- | ---: | ---: |
| $a = b = c$ | 5.758803 &Aring; | 5.762454 &Aring; |
| $\alpha = \beta = \gamma$ | 90&deg; | 90&deg; |
| Cell volume | 190.984 &Aring;<sup>3</sup> | 191.347 &Aring;<sup>3</sup> |
| Density | 5.052 g/cm<sup>3</sup> | 5.042 g/cm<sup>3</sup> |

### Literature comparison

**Computational references.** R&ouml;dl, Furthm&uuml;ller, Suckert, Armuzza, Bechstedt, and Botti, *Phys. Rev. Materials* **3**, 034602 (2019), obtained $a$ = 5.760 &Aring; with PBE for cubic Ge, using VASP with PAW potentials, a 500 eV cutoff, a $12 \times 12 \times 12$ $\Gamma$-centered mesh on the primitive cell, and the Ge 3d electrons in the valence. Haxhijaj, Riemelmoser, and Pasquarello, arXiv:2601.11177 (2026), obtained 5.765 &Aring; with PBE and 5.678 &Aring; with SCAN from Birch-Murnaghan fits, and list earlier SCAN values of 5.658 to 5.673 &Aring;.

**Experimental reference.** Baker and Hart, *Acta Crystallogr. Sect. A* **31**, 364 (1975), measured bulk single-crystal Ge by multiple-beam X-ray diffractometry, comparing the Ge 800 and Si 355 spacings in symmetric transmission with Mo K$\alpha_1$ radiation. At 25 &deg;C,

$$a = 5.6579060 \pm 0.0000092\ \text{\AA}.$$

R&ouml;dl et al. also compile 5.652 &Aring; from X-ray diffraction at 10 K.

| Source | $a$ (&Aring;) | Ultrasoft minus source | ONCV minus source |
| --- | ---: | ---: | ---: |
| R&ouml;dl et al., PBE, VASP PAW | 5.760 | -0.0012 &Aring; (-0.02%) | +0.0025 &Aring; (+0.04%) |
| Haxhijaj et al., PBE | 5.765 | -0.0062 &Aring; (-0.11%) | -0.0025 &Aring; (-0.04%) |
| Haxhijaj et al., SCAN | 5.678 | +0.081 &Aring; (+1.42%) | +0.084 &Aring; (+1.49%) |
| Baker and Hart, bulk crystal, XRD, 25 &deg;C | 5.657906 | +0.1009 &Aring; (+1.78%) | +0.1045 &Aring; (+1.85%) |
| XRD, 10 K (compiled by R&ouml;dl et al.) | 5.652 | +0.107 &Aring; (+1.89%) | +0.110 &Aring; (+1.95%) |

The four PBE values, from two pseudopotentials here and two independent PAW calculations, span only 0.006 &Aring;. The pseudopotential is not the source of the disagreement with experiment. The functional is: every PBE value is about 0.1 &Aring; too large, while SCAN gives 5.678 &Aring;. Haxhijaj et al. connect the size of the PBE error in Ge to the band structure: PBE predicts metallic Ge, metallic bonding is weaker, and the lattice expands. Section 5 shows that metallic state directly.

Temperature accounts for little of the discrepancy. The calculations are static lattices, and the measured lattice parameter changes by only 0.006 &Aring; between 10 K and room temperature. The experimental sample is a bulk single crystal measured by an absolute X-ray method, so no film-strain or nanocrystal correction applies.

## 5. Band gap

### Results

Each gap below is the difference between the lowest unoccupied and highest occupied Kohn-Sham levels printed at the end of the final SCF (64 bands, 56 occupied).

| Calculation | Highest occupied (eV) | Lowest unoccupied (eV) | $E_g$ (eV) |
| --- | ---: | ---: | ---: |
| Ultrasoft, PBE, relaxed cell | 7.7844 | 7.7844 | 0.00 |
| ONCV, PBE, relaxed cell | 8.7892 | 8.7892 | 0.00 |
| ONCV, SCAN (`ge.scf_final.out`) | 7.6353 | 7.5386 | **-0.097** |

For the SCAN calculation,

$$E_g(\text{SCAN}) = 7.5386 - 7.6353 = -0.097\ \text{eV}.$$

A negative gap means the lowest conduction state lies below the highest valence state: the bands overlap by 0.097 eV and the calculation describes a semimetal rather than a semiconductor. This is the same sense in which R&ouml;dl et al. describe the PBEsol band gap of Ge as negative. The SCAN SCF converged in 36 iterations with an $8 \times 8 \times 8$ mesh (65 irreducible k-points), `ecutwfc = 80` Ry, and `ecutrho = 320` Ry. In that input the atomic positions were entered as Cartesian values under the `{crystal}` tag, and Quantum ESPRESSO found 6 symmetry operations for the resulting cell rather than the 24 of the diamond structure. The result is reproducible: an independent run of the same input (`ge.scf_redo.out`) converged in the same 36 iterations to the same band edges, 7.6353 and 7.5386 eV, and the same total energy, -894.08455640 Ry. An earlier run with the atomic mass set to 28.086 (`ge.scf.out`) also gave identical eigenvalues, as expected, because a static SCF does not depend on the mass.

Two further SCAN single points were run with the atoms on the diamond sites, both with the SG15 potential. Each finds all 24 symmetry operations and has the same Ewald energy as the corresponding PBE diamond cell.

| Additional SCAN run | Cell | `ecutwfc` | Highest occupied (eV) | Lowest unoccupied (eV) | $E_g$ (eV) |
| --- | --- | ---: | ---: | ---: | ---: |
| `ge.scf_updated_lattice.out` | ONCV PBE, 5.762454 &Aring; | 80 Ry | 9.8406 | 9.8406 | 0.00 |
| `germanium/scan/scan/ge.scf.out` | Ultrasoft PBE, 5.758803 &Aring; | 56 Ry | 9.8706 | 9.8706 | 0.00 |

### Where the band edges are

In the 8-atom cubic cell, the fcc Brillouin zone is folded into a smaller cubic zone. The three inequivalent $X$ points of fcc fold onto $\Gamma$, and the four inequivalent $L$ points fold onto the cubic zone corner $R = (\tfrac12, \tfrac12, \tfrac12)\,2\pi/a$. The unshifted even meshes used here contain $R$, so the $\Gamma \to L$ gap is sampled. The table lists the eigenvalues of bands 54 to 57 relative to the valence-band maximum (VBM) for the diamond-cell calculations, where the high-symmetry states can be identified.

| State | Ultrasoft, PBE | ONCV, SCAN (`ge.scf_updated_lattice.out`) |
| --- | ---: | ---: |
| s-like state at $\Gamma$ (band 54) | -0.444 eV | -0.411 eV |
| p-like VBM at $\Gamma$, bands 55 to 57, triply degenerate | 0 | 0 |
| Conduction states from $X$, folded to $\Gamma$ (3-fold) | +0.817 eV | +1.069 eV |
| Conduction minimum at $L$, folded to $R$ (4-fold) | +0.027 eV | +0.167 eV |

The ONCV PBE final SCF used 120 irreducible k-points, and Quantum ESPRESSO prints eigenvalues only for fewer than 100 unless `verbosity = 'high'` is set. Its band-resolved values are therefore not available.

The zero gap has a specific origin. In real germanium, the p-like valence-band maximum at $\Gamma$ is fully occupied, the s-like antibonding state at $\Gamma$ lies 0.80 to 0.90 eV above it, and the lowest conduction state is at $L$. In both calculations the s-like state has dropped *below* the p-like triplet. It fills with two electrons, leaving the triply degenerate p level with only four of its six electrons. The highest occupied and lowest unoccupied levels are then members of the same degenerate level, and the system is a semimetal. The $\Gamma \to L$ separation itself stays positive: +0.027 eV in PBE and +0.167 eV in SCAN.

R&ouml;dl et al. report the same ordering with PBEsol and spin-orbit coupling at the PBEsol lattice parameter: the $\Gamma$ conduction state at -0.174 eV, $L$ at +0.005 eV, and the $X$ conduction state at +0.585 eV.

### Comparison with reference values

| Reference | Gap (eV) | PBE (0.00 eV) minus reference | SCAN (-0.097 eV) minus reference |
| --- | ---: | ---: | ---: |
| Experiment, indirect gap at 300 K | 0.66 | -0.66 | -0.76 |
| Experiment, indirect gap at low temperature | 0.744 | -0.74 | -0.84 |
| R&ouml;dl et al., HSE06 ($\Gamma \to L$) | 0.675 | -0.68 | -0.77 |
| Yao and Kanai, SCAN with SCAN pseudopotential | 0.57 | -0.57 | -0.67 |
| Yao and Kanai, SCAN with PBE pseudopotential | 0.19 | -0.19 | -0.29 |
| Yao and Kanai, PBE | 0 | 0.00 | -0.10 |

The 300 K value is the reference given in the project. The low-temperature indirect gap is 0.744 eV; Varshni's fit gives $E_g(0\ \text{K}) = 0.7412$ eV. Yao and Kanai, *J. Chem. Phys.* **146**, 224105 (2017), computed SCAN with Troullier-Martins norm-conserving pseudopotentials. They found that the gap drops from 0.57 eV to 0.19 eV when the pseudopotential is generated with PBE rather than SCAN.

### Discussion

**Functional.** A PBE eigenvalue difference lacks the derivative discontinuity and underestimates gaps. The gap of Ge is small enough that the error closes it and inverts the $\Gamma$ states. Both pseudopotentials give this result, and it matches the literature PBE and PBEsol calculations. It is a property of the functional, not of the setup.

**The SCAN result.** SCAN is expected to improve on PBE for Ge: Yao and Kanai obtain 0.19 to 0.57 eV depending on the pseudopotential. The SCAN value here, -0.097 eV, is further from experiment than the PBE value in raw difference, by 0.10 eV. Two features of the calculation account for that. First, it was run at the PBE lattice parameter, whose effect is estimated below. Second, its cell has lower symmetry than diamond, so the band overlap reflects the structure as entered as well as the functional. The SCAN runs with the atoms on the diamond sites give 0.00 eV at both PBE lattice parameters: the band inversion at $\Gamma$ persists, but the bands no longer overlap.

**Pseudopotential.** Changing from the ultrasoft to the ONCV potential moves the PBE lattice parameter by 0.06% and leaves the gap at zero. For PBE the pseudopotential does not matter at this level. For SCAN it does. The SG15 potential used here was generated with PBE and has no core correction, which is the situation in which Yao and Kanai found only 0.19 eV rather than 0.57 eV.

**Lattice parameter.** This is the largest effect that can be quantified from these results. All of the SCAN calculations were run at a PBE lattice parameter; the reported one used 5.762454 &Aring;, which is 1.85% too large, a volume expansion of 5.65%. The direct gap of Ge is strongly volume dependent; R&ouml;dl et al. quote a volume deformation potential of -9.0 eV for it. A linear estimate gives

$$\Delta E_{\Gamma} \approx -9.0\ \text{eV} \times 0.0565 \approx -0.51\ \text{eV}.$$

The PBE cell therefore lowers the $\Gamma$ conduction state by about half an eV relative to the experimental lattice. The inverted s-like state in the SCAN diamond-cell calculation is only 0.41 eV below the VBM. Running SCAN on a PBE-relaxed cell builds in an error comparable to the inversion itself. The $\Gamma \to L$ separation is less sensitive, which is consistent with SCAN raising $L$ to +0.167 eV even in the expanded cell.

**Other contributions.** Spin-orbit coupling, not included, splits the valence band by 0.30 eV and raises the VBM by about 0.1 eV, lowering all gaps by about that amount. Electron-phonon renormalization lowers the gap by about 0.08 eV between 0 K and 300 K. Both are an order of magnitude smaller than the 0.74 eV discrepancy.

## 6. Runs that did not produce results

- **SCAN with a PAW potential** (ultrasoft workflow): stopped with `Meta-GGA not implemented with USPP/PAW`. No energy or gap.
- **HSE, ultrasoft**: completed the initial PBE SCF and stopped during the exact-exchange stage without `JOB DONE`. The eigenvalues printed in that file are PBE values.
- **HSE, "ONCV" folder**: every submitted run read `ge.scf.in`, which points to the ultrasoft SSSP file; the outputs report `Pseudo is Ultrasoft + core correction`. They stopped at the start of exact exchange, with estimated memory of 96 to 189 GB for the 8-atom cell. The input that uses the SG15 ONCV potential, `ge-Copy1.scf.in`, was never submitted.

## 7. Proposed improvement

Relax the structure with SCAN, using the norm-conserving potential, and take the gap from the final SCF of that relaxation. Section 5 shows that the gap calculation is currently limited by the geometry as well as by the functional. The PBE lattice parameter alone is worth about 0.5 eV at $\Gamma$, which is larger than both the 0.097 eV band overlap in the reported SCAN result and the 0.41 eV inversion in the SCAN diamond-cell calculation. Haxhijaj et al. find a SCAN lattice parameter of 5.678 &Aring;, within 0.4% of experiment, so a SCAN relaxation would remove most of that error at meta-GGA cost. A norm-conserving potential generated for SCAN, rather than for PBE, would further close the 0.19 to 0.57 eV range identified by Yao and Kanai.

The more expensive alternative is an HSE06 SCF. It should use the ONCV input that was never submitted, and a 2-atom primitive cell instead of the 8-atom cube. The primitive cell has one quarter of the atoms, bands, and volume, and the cost and memory of exact exchange grow faster than linearly in all three; the estimated 96 to 189 GB comes largely from the 8-atom cell.

## References

1. J. P. Perdew, K. Burke, and M. Ernzerhof, "Generalized Gradient Approximation Made Simple," *Phys. Rev. Lett.* **77**, 3865 (1996).
2. J. Sun, A. Ruzsinszky, and J. P. Perdew, "Strongly Constrained and Appropriately Normed Semilocal Density Functional," *Phys. Rev. Lett.* **115**, 036402 (2015).
3. Y. Yao and Y. Kanai, "Plane-wave pseudopotential implementation and performance of SCAN meta-GGA exchange-correlation functional for extended systems," *J. Chem. Phys.* **146**, 224105 (2017).
4. C. R&ouml;dl, J. Furthm&uuml;ller, J. R. Suckert, V. Armuzza, F. Bechstedt, and S. Botti, "Accurate electronic and optical properties of hexagonal germanium for optoelectronic applications," *Phys. Rev. Materials* **3**, 034602 (2019).
5. A. Haxhijaj, S. Riemelmoser, and A. Pasquarello, "Lattice dynamics and structural phase stability of group-IV elemental solids with the r2SCAN functional," arXiv:2601.11177 (2026).
6. J. F. C. Baker and M. Hart, "An absolute measurement of the lattice parameter of germanium using multiple-beam X-ray diffractometry," *Acta Crystallogr. Sect. A* **31**, 364 (1975).
7. Y. P. Varshni, "Temperature dependence of the energy gap in semiconductors," *Physica* **34**, 149 (1967).
8. D. R. Hamann, "Optimized norm-conserving Vanderbilt pseudopotentials," *Phys. Rev. B* **88**, 085117 (2013).
9. M. Schlipf and F. Gygi, "Optimization algorithm for the generation of ONCV pseudopotential libraries for density functional theory," *Comput. Phys. Commun.* **196**, 36 (2015).
10. G. Prandini, A. Marrazzo, I. E. Castelli, N. Mounet, and N. Marzari, "Precision and efficiency in solid-state pseudopotential calculations," *npj Comput. Mater.* **4**, 72 (2018).
11. P. Giannozzi et al., "Advanced capabilities for materials modelling with Quantum ESPRESSO," *J. Phys.: Condens. Matter* **29**, 465901 (2017).
