"""Kohn-Sham eigenvalues near the Fermi level at Gamma and R for each Ge calculation.

In the 8-atom conventional cubic cell, the fcc X points fold onto Gamma and the
fcc L points fold onto the cubic zone corner R = (1/2, 1/2, 1/2) 2pi/a.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "week4"
N_OCC = 56  # 112 valence electrons, spin-unpolarized

RUNS = {
    "Ultrasoft PBE": ROOT / "germanium/germanium/Runtime/ge_scf_bandgap.out",
    "ONCV PBE": ROOT / "germanium-ONCV/germanium-ONCV/ge.relax.out",
    "ONCV SCAN (corrected)": ROOT / "germanium-ONCV/scan/scan/ge.scf_updated_lattice.out",
}

NUM = re.compile(r"-?\d+\.\d+")
KLINE = re.compile(r"^\s*k =\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)")


def last_band_block(path):
    lines = path.read_text(errors="replace").splitlines()
    end = max(i for i, l in enumerate(lines) if "highest occupied" in l)
    start = max(i for i, l in enumerate(lines[:end]) if "End of self-consistent" in l or "End of band structure" in l)
    bands = {}
    kpt = None
    for line in lines[start:end]:
        m = KLINE.match(line)
        if m:
            kpt = tuple(round(abs(float(x)), 4) for x in m.groups())
            bands[kpt] = []
            continue
        if kpt is not None and line.strip() and "occupation" not in line:
            bands[kpt].extend(float(x) for x in NUM.findall(line))
    edges = [float(x) for x in NUM.findall(lines[end])]
    return bands, edges


for label, path in RUNS.items():
    bands, edges = last_band_block(path)
    vbm = edges[0]
    print(f"== {label}: highest occupied / lowest unoccupied = {edges}")
    if not bands:
        print("  eigenvalues not printed (more than 100 k-points without verbosity='high')")
        continue
    for name, k in [("Gamma", (0.0, 0.0, 0.0)), ("R (fcc L)", (0.5, 0.5, 0.5))]:
        e = bands[k]
        window = e[N_OCC - 4:N_OCC + 4]
        rel = " ".join(f"{x - vbm:+.3f}" for x in window)
        print(f"  {name:10s} bands {N_OCC - 3}-{N_OCC + 4} relative to VBM: {rel}")
    lumo = min(e[N_OCC] for e in bands.values())
    homo = max(e[N_OCC - 1] for e in bands.values())
    print(f"  max band {N_OCC} = {homo:.4f}, min band {N_OCC + 1} = {lumo:.4f}, gap = {lumo - homo:+.4f} eV")
