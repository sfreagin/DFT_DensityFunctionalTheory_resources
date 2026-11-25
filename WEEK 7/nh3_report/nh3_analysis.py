"""Figures and tables for the Week 7 NH3 IR report."""

import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplconfig"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np

WEEK7 = HERE.parent
ECUT_DAT = WEEK7 / "NH3/relaxed/optimize_ecutoff/ecut.dat"
NH3_IN = WEEK7 / "NH3/phonon/nh3.in"
QPOINTS = WEEK7 / "NH3/phonon/qpoints.yaml"
OSC = WEEK7 / "NH3/combined_files/OscillatorStrengths.txt"
PH_OUT = WEEK7 / "NH3/dielectric/nh3.ph.out"

RY_TO_MEV = 13605.693009
N_ATOMS = 4
BOX = 15.0
SELECTED_ECUT = 70
MASSES = np.array([14.0067, 1.00794, 1.00794, 1.00794])
# |dmu/dQ|^2 in e^2/amu -> IR intensity in km/mol: (4.803204 D/(e A))^2 * 42.2561
E2_AMU_TO_KM_MOL = 4.803204**2 * 42.2561

# Shimanouchi (NSRDS-NBS 39) fundamentals; Koops, Visser, and Smit (1983) band intensities.
EXPERIMENT = [
    ("Umbrella (sym. deformation)", "A1", 950.0, 138.0),
    ("Bend (deg. deformation)", "E", 1627.0, 28.2),
    ("Symmetric N-H stretch", "A1", 3337.0, 7.1),
    ("Asymmetric N-H stretch", "E", 3444.0, 4.3),
]
EXP_RNH, EXP_HNH = 1.012, 106.67  # Herzberg (1966)
EXP_ALPHA = 2.103  # A^3, Olney et al. (1997)

COLOR, EXP_COLOR, ACCENT = "#1f4e79", "#8a8a8a", "#c45c26"


def read_ecut():
    energy = {}
    for line in ECUT_DAT.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            energy[int(parts[0])] = float(parts[1])
    ecut = sorted(energy)
    return np.array(ecut), np.array([energy[c] for c in ecut])


def read_positions():
    lines = NH3_IN.read_text().splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip().startswith("ATOMIC_POSITIONS"))
    rows = [l.split() for l in lines[start + 1 : start + 1 + N_ATOMS]]
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows]) * BOX


def read_modes():
    bands = re.split(r"\n  - # \d+\n", QPOINTS.read_text())[1:]
    freqs, vecs = [], []
    for band in bands:
        freqs.append(float(re.search(r"frequency:\s+(\S+)", band).group(1)))
        re_parts = [float(v) for v in re.findall(r"- \[\s*(-?[\d.]+),", band)][: 3 * N_ATOMS]
        vecs.append(np.array(re_parts).reshape(N_ATOMS, 3))
    return np.array(freqs), vecs


def read_oscillator():
    text = OSC.read_text()
    head, tail = text.split("Oscillator Strengths:")
    num = r"-\s+(-?[\d.]+(?:e-?\d+)?)"
    return (np.array([float(x) for x in re.findall(num, head)]),
            np.array([float(x) for x in re.findall(num, tail)]))


def read_dielectric():
    lines = PH_OUT.read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if "Dielectric constant in cartesian axis" in l)
    eps = np.array([[float(x) for x in re.findall(r"-?\d+\.\d+", lines[i + 2 + r])] for r in range(3)])
    j = next(k for k, l in enumerate(lines) if "with asr applied" in l)
    zstar = []
    for a in range(N_ATOMS):
        block = lines[j + 1 + 4 * a : j + 5 + 4 * a]
        zstar.append(np.array([[float(x) for x in re.findall(r"-?\d+\.\d+", row)] for row in block[1:]]))
    return eps, zstar


def geometry(x):
    bonds = np.array([np.linalg.norm(x[i] - x[0]) for i in (1, 2, 3)])
    angles = []
    for i, j in ((1, 2), (1, 3), (2, 3)):
        a, b = x[i] - x[0], x[j] - x[0]
        angles.append(np.degrees(np.arccos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))))
    return bonds, np.array(angles)


def plot_ecut(ecut, energy):
    delta = (energy[:-1] - energy[-1]) * RY_TO_MEV / N_ATOMS
    fig, ax = plt.subplots(figsize=(6.0, 3.8), constrained_layout=True)
    ax.semilogy(ecut[:-1], delta, "o-", color=COLOR, ms=4.5, lw=1.4)
    ax.axvline(SELECTED_ECUT, color=ACCENT, ls="--", lw=1.2, label=f"Selected cutoff ({SELECTED_ECUT} Ry)")
    ax.axhline(1.0, color="0.45", ls=":", lw=1.0, label="1 meV/atom")
    ax.set_xlabel(r"Wave-function cutoff, $E_{\mathrm{cut}}^{\mathrm{wfc}}$ (Ry)")
    ax.set_ylabel(f"Energy above {ecut[-1]} Ry (meV/atom)")
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(HERE / "nh3_ecut_convergence.png", dpi=300)
    fig.savefig(HERE / "nh3_ecut_convergence.pdf")
    plt.close(fig)


def project(x, elev, azim):
    e, a = np.radians(elev), np.radians(azim)
    right = np.array([-np.sin(a), np.cos(a), 0.0])
    up = np.array([-np.sin(e) * np.cos(a), -np.sin(e) * np.sin(a), np.cos(e)])
    return np.stack([x @ right, x @ up], axis=-1)


def plot_modes(symbols, x0, freqs, vecs, picks):
    fig, axes = plt.subplots(1, 4, figsize=(11.0, 3.4), constrained_layout=True)
    center = x0.mean(axis=0)
    for ax, (k, title, side) in zip(axes, picks):
        elev, azim = (15, -15) if side else (90, -90)
        u = vecs[k] / np.sqrt(MASSES)[:, None]
        u = 0.55 * u / np.abs(u).max()
        p = project(x0 - center, elev, azim)
        d = project(u, elev, azim)
        for i in (1, 2, 3):
            ax.plot(*zip(p[0], p[i]), color="0.55", lw=2.2, zorder=1)
        for i, s in enumerate(symbols):
            ax.scatter(*p[i], s=520 if s == "N" else 260, color="#3b6fb6" if s == "N" else "0.92",
                       edgecolor="0.25", lw=1.0, zorder=2)
            ax.annotate("", xy=p[i] + d[i], xytext=p[i],
                        arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=1.6), zorder=3)
        ax.set_title(f"{title}\n{freqs[k]:.0f} cm$^{{-1}}$", fontsize=9)
        ax.set_aspect("equal")
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.4, 1.4)
        ax.axis("off")
        ax.text(0.5, -0.02, "side view" if side else "view down the C$_3$ axis",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5, color="0.35")
    fig.savefig(HERE / "nh3_mode_sketches.png", dpi=300)
    fig.savefig(HERE / "nh3_mode_sketches.pdf")
    plt.close(fig)


def plot_ir(freqs, strengths, labels, exp_bands, sigma=20.0):
    grid = np.linspace(500, 4000, 7001)
    spectrum = sum(s * np.exp(-0.5 * ((grid - f) / sigma) ** 2) for f, s in zip(freqs, strengths))
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 6.2), constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[3, 2]))
    ax = axes[0]
    ax.plot(grid, spectrum, color=COLOR, lw=1.4)
    ax.vlines(freqs, 0, strengths, color=COLOR, lw=0.8, alpha=0.5)
    for f, s, text, ha in labels:
        dx = {"center": 0, "right": -3, "left": 3}[ha]
        ax.annotate(f"{text}\n{f:.0f} cm$^{{-1}}$", xy=(f, s), xytext=(dx, 8), textcoords="offset points",
                    ha=ha, va="bottom", fontsize=8)
    ax.set_xlim(4000, 500)
    ax.set_ylim(0, 1.32 * strengths.max())
    ax.set_xlabel(r"Wavenumber (cm$^{-1}$)")
    ax.set_ylabel(r"Oscillator strength (e$^2$/amu)")
    ax.set_title(f"(a) Calculated IR spectrum (Gaussian broadening, $\\sigma$ = {sigma:.0f} cm$^{{-1}}$)", fontsize=10)
    ax.tick_params(direction="in", top=True, right=True)

    ax = axes[1]
    calc_f = [b[0] for b in exp_bands]
    calc_i = [b[1] for b in exp_bands]
    exp_f = [b[2] for b in exp_bands]
    exp_i = [b[3] for b in exp_bands]
    ax.vlines(calc_f, 0, calc_i, color=COLOR, lw=3.0, label="This work (PBE, harmonic)")
    ax.vlines(exp_f, 0, exp_i, color=EXP_COLOR, lw=3.0, label="Experiment (gas)")
    ax.set_xlim(4000, 500)
    ax.set_yscale("log")
    ax.set_ylim(1, 400)
    ax.set_xlabel(r"Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Band intensity (km/mol)")
    ax.set_title("(b) Band positions and integrated intensities", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.tick_params(direction="in", top=True, right=True, which="both")
    fig.savefig(HERE / "nh3_ir_spectrum.png", dpi=300)
    fig.savefig(HERE / "nh3_ir_spectrum.pdf")
    plt.close(fig)


def main():
    ecut, energy = read_ecut()
    plot_ecut(ecut, energy)
    print("ecutwfc  E (Ry)  above max (meV/atom)")
    for c, e in zip(ecut, energy):
        print(f"{c:5d}  {e:.8f}  {(e - energy[-1]) * RY_TO_MEV / N_ATOMS:9.3f}")

    symbols, x0 = read_positions()
    bonds, angles = geometry(x0)
    print(f"\nrelaxed r(NH) = {bonds.mean():.4f} A ({(bonds.mean() / EXP_RNH - 1) * 100:+.2f}% vs {EXP_RNH})")
    print(f"relaxed HNH = {angles.mean():.2f} deg ({angles.mean() - EXP_HNH:+.2f} deg vs {EXP_HNH})")

    freqs, vecs = read_modes()
    osc_f, osc_s = read_oscillator()
    ratio = osc_f[6:] / freqs[6:]
    print(f"\nIR-code / qpoints frequency ratio: {ratio.min():.5f} to {ratio.max():.5f}")

    groups = [(6,), (7, 8), (9,), (10, 11)]
    peak_labels = [("Umbrella (A$_1$)", "center"), ("Bend (E)", "center"),
                   ("Sym. stretch (A$_1$)", "left"), ("Asym. stretch (E)", "right")]
    print("\nmode  calc (cm-1)  exp (cm-1)  error  strength (e2/amu)  calc (km/mol)  exp (km/mol)")
    exp_bands, labels = [], []
    for g, (name, _, f_exp, i_exp), (text, ha) in zip(groups, EXPERIMENT, peak_labels):
        f = freqs[list(g)].mean()
        s = osc_s[list(g)].sum()
        exp_bands.append((f, s * E2_AMU_TO_KM_MOL, f_exp, i_exp))
        labels.append((f, s, text, ha))
        print(f"{name:30s} {f:8.1f}  {f_exp:7.0f}  {(f / f_exp - 1) * 100:+5.1f}%  "
              f"{s:.5f}  {s * E2_AMU_TO_KM_MOL:7.1f}  {i_exp:6.1f}")

    picks = [(6, "Umbrella, A$_1$", True), (7, "Deg. deformation, E", False),
             (9, "Sym. stretch, A$_1$", True), (10, "Asym. stretch, E", False)]
    plot_modes(symbols, x0, freqs, vecs, picks)

    keep = freqs > 500
    plot_ir(freqs[keep], osc_s[keep], labels, exp_bands)

    eps, zstar = read_dielectric()
    volume = BOX**3
    alpha = 3 * volume / (4 * np.pi) * (np.diag(eps) - 1) / (np.diag(eps) + 2)
    print(f"\nepsilon diag = {np.round(np.diag(eps), 6)}")
    print(f"alpha (Clausius-Mossotti) = {np.round(alpha, 3)} A^3, mean {alpha.mean():.3f} "
          f"({(alpha.mean() / EXP_ALPHA - 1) * 100:+.1f}% vs {EXP_ALPHA})")
    for s, z in zip(symbols, zstar):
        print(f"Z*({s}) mean {np.trace(z) / 3:+.4f}  diag {np.round(np.diag(z), 4)}")
    print(f"sum of Z* (should be 0): {np.round(sum(zstar), 5).tolist()}")


if __name__ == "__main__":
    main()
