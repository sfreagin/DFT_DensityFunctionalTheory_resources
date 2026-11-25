"""Energy-cutoff convergence plot for the Week 4 germanium project."""

import os
from pathlib import Path

import matplotlib.pyplot as plt

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

# Total energies from Runtime/ecut.dat (8-atom conventional cell, 4x4x4 mesh).
ecut = [24, 28, 32, 36, 40, 44, 48, 52, 56, 60]
energy = [
    -1707.72654357,
    -1708.88691527,
    -1709.15054545,
    -1709.19327673,
    -1709.19589391,
    -1709.19664908,
    -1709.19839129,
    -1709.19941308,
    -1709.19956510,
    -1709.19957118,
]
selected = 56
n_atoms = 8
ry_to_mev = 13605.693009  # 1 Ry in meV

e_ref = energy[-1]
delta_mev_per_atom = [(e - e_ref) * ry_to_mev / n_atoms for e in energy]

fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)

ax = axes[0]
ax.plot(ecut, energy, "o-", color="#1f4e79", ms=5.5, lw=1.4, label="SCF total energy")
ax.axvline(selected, color="#c45c26", ls="--", lw=1.2, label=f"Selected cutoff ({selected} Ry)")
ax.set_xlabel(r"Wave-function cutoff, $E_{\mathrm{cut}}^{\mathrm{wfc}}$ (Ry)")
ax.set_ylabel("Total energy (Ry)")
ax.set_title("(a) Total energy")
ax.legend(frameon=False, fontsize=8)

ax = axes[1]
ax.semilogy(ecut[:-1], delta_mev_per_atom[:-1], "o-", color="#1f4e79", ms=5.5, lw=1.4)
ax.axvline(selected, color="#c45c26", ls="--", lw=1.2)
ax.axhline(1.0, color="0.45", ls=":", lw=1.0, label="1 meV/atom")
ax.set_xlabel(r"Wave-function cutoff, $E_{\mathrm{cut}}^{\mathrm{wfc}}$ (Ry)")
ax.set_ylabel(r"$|\Delta E|$ per atom relative to 60 Ry (meV)")
ax.set_title("(b) Energy change per atom")
ax.set_ylim(3e-3, 8e3)
ax.legend(frameon=False, fontsize=8)

for axis in axes:
    axis.set_xlim(20, 64)
    axis.tick_params(direction="in", top=True, right=True)

out = Path(__file__).resolve().parent
fig.savefig(out / "ge_ecut_convergence.png", dpi=300)
fig.savefig(out / "ge_ecut_convergence.pdf")

print("cutoff  E (Ry)  dE from 60 Ry (meV/atom)  step to next (meV/atom)")
for i, (c, e, d) in enumerate(zip(ecut, energy, delta_mev_per_atom)):
    if i < len(energy) - 1:
        step = (energy[i + 1] - e) * ry_to_mev / n_atoms
        print(f"{c:6.0f}  {e:15.8f}  {d:10.4f}  {step:10.4f}")
    else:
        print(f"{c:6.0f}  {e:15.8f}  {d:10.4f}")
