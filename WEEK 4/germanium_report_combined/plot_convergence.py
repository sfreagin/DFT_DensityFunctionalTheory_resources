"""Cutoff and k-point convergence for the ultrasoft and ONCV germanium calculations."""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

import matplotlib.pyplot as plt

RY_TO_MEV = 13605.693009
N_ATOMS = 8

# Ultrasoft ge_pbe_v1.4.uspp.F.UPF: Runtime/ecut.dat (4x4x4) and Runtime/kpoints.dat (56 Ry).
us_ecut = [24, 28, 32, 36, 40, 44, 48, 52, 56, 60]
us_e = [-1707.72654357, -1708.88691527, -1709.15054545, -1709.19327673, -1709.19589391,
        -1709.19664908, -1709.19839129, -1709.19941308, -1709.19956510, -1709.19957118]
us_k = [2, 4, 6, 8, 10, 12, 14]
us_ke = [-1709.11869389, -1709.19956510, -1709.20315909, -1709.20355376, -1709.20362832,
         -1709.20364874, -1709.20365611]

# ONCV Ge.upf: ecut.dat (8x8x8) and kpoints.dat (80 Ry).
nc_ecut = [24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92]
nc_e = [-1286.67633147, -1334.51892464, -1369.37983387, -1393.57972872, -1409.34550293,
        -1418.96623924, -1424.33461983, -1427.07893014, -1428.34519899, -1428.86930385,
        -1429.05957264, -1429.11863429, -1429.13361283, -1429.13636913, -1429.13663912,
        -1429.13665959, -1429.13669724, -1429.13672220]
nc_k = [4, 8, 10, 12, 14]
nc_ke = [-1429.13263976, -1429.13663912, -1429.13671387, -1429.13673444, -1429.13674182]


def rel(energies):
    ref = energies[-1]
    return [(e - ref) * RY_TO_MEV / N_ATOMS for e in energies[:-1]]


US_COLOR, NC_COLOR = "#1f4e79", "#b5462a"
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9), constrained_layout=True)

ax = axes[0]
ax.semilogy(us_ecut[:-1], rel(us_e), "o-", color=US_COLOR, ms=4.5, lw=1.4,
            label="Ultrasoft (vs. 60 Ry)")
ax.semilogy(nc_ecut[:-1], rel(nc_e), "s-", color=NC_COLOR, ms=4.5, lw=1.4,
            label="ONCV (vs. 92 Ry)")
ax.axvline(56, color=US_COLOR, ls="--", lw=1.0)
ax.axvline(80, color=NC_COLOR, ls="--", lw=1.0)
ax.axhline(1.0, color="0.45", ls=":", lw=1.0, label="1 meV/atom")
ax.text(56.8, 3e4, "56 Ry", color=US_COLOR, fontsize=8)
ax.text(80.8, 3e4, "80 Ry", color=NC_COLOR, fontsize=8)
ax.set_xlabel(r"Wave-function cutoff, $E_{\mathrm{cut}}^{\mathrm{wfc}}$ (Ry)")
ax.set_ylabel(r"$\Delta E$ per atom (meV)")
ax.set_title("(a) Cutoff convergence")
ax.set_xlim(20, 96)
ax.set_ylim(3e-3, 3e5)
ax.legend(frameon=False, fontsize=8, loc="lower left")

ax = axes[1]
ax.semilogy(us_k[:-1], rel(us_ke), "o-", color=US_COLOR, ms=4.5, lw=1.4,
            label="Ultrasoft, 56 Ry")
ax.semilogy(nc_k[:-1], rel(nc_ke), "s-", color=NC_COLOR, ms=4.5, lw=1.4,
            label="ONCV, 80 Ry")
ax.axhline(1.0, color="0.45", ls=":", lw=1.0, label="1 meV/atom")
ax.set_xlabel(r"Monkhorst-Pack mesh, $n \times n \times n$")
ax.set_ylabel(r"$\Delta E$ per atom relative to $14^3$ (meV)")
ax.set_title("(b) k-point convergence")
ax.set_xticks([2, 4, 6, 8, 10, 12])
ax.set_ylim(5e-3, 5e2)
ax.legend(frameon=False, fontsize=8)

for axis in axes:
    axis.tick_params(direction="in", top=True, right=True, which="both")

out = Path(__file__).resolve().parent
fig.savefig(out / "ge_convergence_us_vs_oncv.png", dpi=300)
fig.savefig(out / "ge_convergence_us_vs_oncv.pdf")

a_us, a_nc, a_exp = 5.758802973, 5.762454494, 5.6579060
for name, a in [("US", a_us), ("ONCV", a_nc)]:
    dv = (a / a_exp) ** 3 - 1
    print(f"{name}: da = {a - a_exp:+.4f} A ({(a - a_exp) / a_exp * 100:+.2f}%), "
          f"dV/V = {dv * 100:.2f}%, -9.0 eV * dV/V = {-9.0 * dv:+.3f} eV")
print(f"ONCV - US: {a_nc - a_us:+.4f} A ({(a_nc - a_us) / a_us * 100:+.3f}%)")
