"""Energy-cutoff convergence plot for the ONCV germanium calculations."""

import os
from pathlib import Path

import matplotlib.pyplot as plt

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

# Total energies from ecut.dat (8-atom cell, 8x8x8 mesh, ecutrho = 8*ecutwfc).
ecut = [24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92]
energy = [
    -1286.67633147,
    -1334.51892464,
    -1369.37983387,
    -1393.57972872,
    -1409.34550293,
    -1418.96623924,
    -1424.33461983,
    -1427.07893014,
    -1428.34519899,
    -1428.86930385,
    -1429.05957264,
    -1429.11863429,
    -1429.13361283,
    -1429.13636913,
    -1429.13663912,
    -1429.13665959,
    -1429.13669724,
    -1429.13672220,
]
selected = 80
n_atoms = 8
ry_to_mev = 13605.693009

e_ref = energy[-1]
delta = [(e - e_ref) * ry_to_mev / n_atoms for e in energy]

fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.8), constrained_layout=True)

ax = axes[0]
ax.plot(ecut, energy, "o-", color="#1f4e79", ms=4.5, lw=1.4, label="SCF total energy")
ax.axvline(selected, color="#c45c26", ls="--", lw=1.2, label=f"Selected cutoff ({selected} Ry)")
ax.set_xlabel(r"Wave-function cutoff, $E_{\mathrm{cut}}^{\mathrm{wfc}}$ (Ry)")
ax.set_ylabel("Total energy (Ry)")
ax.set_title("(a) Total energy")
ax.legend(frameon=False, fontsize=8)

ax = axes[1]
ax.semilogy(ecut[:-1], delta[:-1], "o-", color="#1f4e79", ms=4.5, lw=1.4)
ax.axvline(selected, color="#c45c26", ls="--", lw=1.2)
ax.axhline(1.0, color="0.45", ls=":", lw=1.0, label="1 meV/atom")
ax.set_xlabel(r"Wave-function cutoff, $E_{\mathrm{cut}}^{\mathrm{wfc}}$ (Ry)")
ax.set_ylabel(r"$\Delta E$ per atom relative to 92 Ry (meV)")
ax.set_title("(b) Energy change per atom")
ax.set_ylim(3e-3, 3e5)
ax.legend(frameon=False, fontsize=8)

for axis in axes:
    axis.set_xlim(20, 96)
    axis.tick_params(direction="in", top=True, right=True)

out = Path(__file__).resolve().parent
fig.savefig(out / "ge_oncv_ecut_convergence.png", dpi=300)
fig.savefig(out / "ge_oncv_ecut_convergence.pdf")

print("cutoff  E (Ry)  dE from 92 Ry (meV/atom)  step (meV/atom)")
for i, (c, e, d) in enumerate(zip(ecut, energy, delta)):
    step = (energy[i + 1] - e) * ry_to_mev / n_atoms if i < len(energy) - 1 else None
    print(f"{c:6.0f}  {e:15.8f}  {d:12.4f}  {step}")

print("--- k mesh at 80 Ry, relative to 14x14x14 ---")
ks = [4, 8, 10, 12, 14]
kes = [-1429.13263976, -1429.13663912, -1429.13671387, -1429.13673444, -1429.13674182]
kref = kes[-1]
for k, e in zip(ks, kes):
    print(f"{k:3d}  {e:.8f}  {(e - kref) * ry_to_mev / n_atoms:10.4f}")

a = 5.762454494
exp = 5.6579060
lit = 5.760
print("delta exp A", a - exp, "pct", (a - exp) / exp * 100)
print("delta lit A", a - lit, "pct", (a - lit) / lit * 100)
print("scan gap", 7.5386 - 7.6353)
