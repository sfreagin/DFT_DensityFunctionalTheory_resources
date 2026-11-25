#!/usr/bin/env python3
"""
Week 9 comparison figures for the report.

Part 1: time-step test (0.5, 1, 2, 4 fs, NVT, 20 ps production each)
Part 2: NVE / NVT / NPT branches (2 fs, 50 ps production each)

Run from the WEEK 9 folder (or anywhere; paths are relative to this script):

    python make_comparison_plots.py

Outputs (written to comparison_plots/):
    timestep_thermo.png     T, PE, Etot vs time for the four time steps
    timestep_energies.png   PE, KE, Etot vs time (raw + running average)
    timestep_cost.png       mean PE vs time step, and wall time per ps
    timestep_msd.png        oxygen MSD for the four time steps
    ensemble_msd.png        oxygen MSD + linear fits for NVE/NVT/NPT
    ensemble_thermo.png     T, density, Etot vs time for NVE/NVT/NPT
    timestep_summary.csv
    ensemble_summary.csv
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / "comparison_plots"

TIMESTEP_DIRS = {
    0.5: "week8_results_0.5sec",
    1.0: "week8_results_1sec",
    2.0: "week8_results_2sec",
    4.0: "week8_results_4sec",
}
ENSEMBLES = ["NVE", "NVT", "NPT"]
ENSEMBLE_DT_FS = 2.0

TIMESTEP_FIT_PS = (2.0, 20.0)
ENSEMBLE_FIT_PS = (2.0, 40.0)
N_WINDOWS = 5
N_BLOCKS = 5
RUNNING_AVG_PS = 2.0


# ------------------------------------------------------------
# Readers
# ------------------------------------------------------------

def read_production_thermo(log_path):
    """
    Return the last thermo block in a LAMMPS log (the production run)
    as {column_name: array}, plus that run's loop time in seconds.
    """
    blocks = []
    header, rows = None, []
    with open(log_path, errors="ignore") as f:
        for line in f:
            parts = line.split()
            if parts and parts[0] == "Step":
                header, rows = parts, []
            elif header and line.startswith("Loop time"):
                blocks.append((header, rows, float(parts[3])))
                header = None
            elif header and len(parts) == len(header):
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    pass

    if not blocks:
        raise RuntimeError(f"No completed thermo block found in {log_path}")

    header, rows, loop_time = blocks[-1]
    data = np.asarray(rows)
    thermo = {name: data[:, i] for i, name in enumerate(header)}
    thermo["t_ps"] = (thermo["Time"] - thermo["Time"][0]) / 1000.0
    return thermo, loop_time


def read_msd(path, dt_fs):
    step, msd = np.loadtxt(path, comments="#", unpack=True)
    t_ps = (step - step[0]) * dt_fs / 1000.0
    return t_ps, msd


def find_ensemble_msd(folder, ensemble):
    for name in (f"msd_{ensemble}.dat", "msd_NVT.dat"):
        path = folder / name
        if path.exists():
            return path
    raise FileNotFoundError(f"No MSD file for {ensemble} in {folder}")


# ------------------------------------------------------------
# Analysis helpers
# ------------------------------------------------------------

def fit_diffusion(t, msd, window):
    """Linear fit of MSD over window (ps). Returns D (m^2/s), slope, intercept."""
    mask = (t >= window[0]) & (t <= window[1])
    slope, intercept = np.polyfit(t[mask], msd[mask], 1)
    return slope / 6.0 * 1.0e-8, slope, intercept


def window_sem(t, msd, window, n_windows=N_WINDOWS):
    """Same sub-window spread used by analyze_msd.py."""
    edges = np.linspace(window[0], window[1], n_windows + 1)
    values = []
    for i in range(n_windows):
        lo, hi = edges[i], edges[i + 1]
        mask = (t >= lo) & ((t <= hi) if i == n_windows - 1 else (t < hi))
        if mask.sum() >= 3:
            values.append(np.polyfit(t[mask], msd[mask], 1)[0] / 6.0 * 1.0e-8)
    if len(values) < 2:
        return np.nan
    return np.std(values, ddof=1) / np.sqrt(len(values))


def block_sem(x, n_blocks=N_BLOCKS):
    """Standard error of the mean from n_blocks contiguous block averages."""
    blocks = np.array_split(np.asarray(x), n_blocks)
    means = np.array([b.mean() for b in blocks])
    return means.std(ddof=1) / np.sqrt(n_blocks)


def drift_per_ps(t, y):
    return np.polyfit(t, y, 1)[0]


# ------------------------------------------------------------
# Part 1: time-step comparison
# ------------------------------------------------------------

def timestep_comparison():
    results = []
    for dt, dirname in TIMESTEP_DIRS.items():
        folder = HERE / dirname
        thermo, loop_time = read_production_thermo(folder / "log.lammps")
        t_msd, msd = read_msd(folder / "msd_O_512.dat", dt)
        d, _, _ = fit_diffusion(t_msd, msd, TIMESTEP_FIT_PS)
        results.append({
            "dt": dt, "thermo": thermo, "loop_time": loop_time,
            "t_msd": t_msd, "msd": msd, "D": d,
            "D_sem": window_sem(t_msd, msd, TIMESTEP_FIT_PS),
        })

    # T, PE, Etot vs time
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    for r in results:
        th, label = r["thermo"], f"{r['dt']:g} fs"
        axes[0].plot(th["t_ps"], th["Temp"], lw=0.9, label=label)
        axes[1].plot(th["t_ps"], th["PotEng"], lw=0.9,
                     label=f"{label} (mean {th['PotEng'].mean():.0f})")
        axes[2].plot(th["t_ps"], th["TotEng"], lw=0.9, label=label)
    axes[0].set_ylabel("Temperature (K)")
    axes[1].set_ylabel("Potential energy\n(kcal/mol)")
    axes[2].set_ylabel("Total energy\n(kcal/mol)")
    axes[2].set_xlabel("Production time (ps)")
    axes[0].set_title("Time-step test (NVT, 300 K, 20 ps production)")
    for ax in axes:
        ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "timestep_thermo.png", dpi=220)
    plt.close(fig)

    # PE, KE, Etot vs time: raw (faint) + running average (bold)
    fig, axes = plt.subplots(3, 1, figsize=(8, 9.5), sharex=True)
    energy_cols = [("PotEng", "Potential energy"), ("KinEng", "Kinetic energy"),
                   ("TotEng", "Total energy")]
    for r in results:
        th = r["thermo"]
        n_avg = max(1, int(round(RUNNING_AVG_PS / (th["t_ps"][1] - th["t_ps"][0]))))
        for ax, (col, _) in zip(axes, energy_cols):
            line, = ax.plot(th["t_ps"], th[col], lw=0.6, alpha=0.35)
            smooth = np.convolve(th[col], np.ones(n_avg) / n_avg, mode="valid")
            t_smooth = th["t_ps"][(n_avg - 1) // 2:(n_avg - 1) // 2 + len(smooth)]
            ax.plot(t_smooth, smooth, lw=1.8, color=line.get_color(),
                    label=f"{r['dt']:g} fs (mean {th[col].mean():.0f}, std {th[col].std():.0f})")
    for ax, (_, name) in zip(axes, energy_cols):
        ax.set_ylabel(f"{name}\n(kcal/mol)")
        ax.legend(fontsize=8, ncol=2)
    axes[2].set_xlabel("Production time (ps)")
    axes[0].set_title(f"NVT energies vs time step (faint: raw, bold: {RUNNING_AVG_PS:g} ps running average)",
                      fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "timestep_energies.png", dpi=220)
    plt.close(fig)

    # Mean PE vs dt, and cost
    dts = np.array([r["dt"] for r in results])
    pe_mean = np.array([r["thermo"]["PotEng"].mean() for r in results])
    pe_sem = np.array([block_sem(r["thermo"]["PotEng"]) for r in results])
    prod_ps = np.array([r["thermo"]["t_ps"][-1] for r in results])
    wall_per_ps = np.array([r["loop_time"] for r in results]) / prod_ps

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    ax1.errorbar(dts, pe_mean, yerr=pe_sem, marker="o", capsize=4)
    ax1.set_xscale("log", base=2)
    ax1.set_xticks(dts)
    ax1.set_xticklabels([f"{d:g}" for d in dts])
    ax1.set_xlabel("Time step (fs)")
    ax1.set_ylabel("Mean potential energy (kcal/mol)")
    ax1.set_title(f"Mean PE ({N_BLOCKS}-block SEM)")
    ax2.bar([f"{d:g} fs" for d in dts], wall_per_ps)
    for i, v in enumerate(wall_per_ps):
        ax2.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=9)
    ax2.set_ylabel("Wall time per simulated ps (s)")
    ax2.set_title("Computational cost (2 procs)")
    fig.tight_layout()
    fig.savefig(OUT / "timestep_cost.png", dpi=220)
    plt.close(fig)

    # MSD overlay
    fig = plt.figure(figsize=(7.4, 4.8))
    for r in results:
        plt.plot(r["t_msd"], r["msd"], lw=1.2,
                 label=f"{r['dt']:g} fs: D = {r['D'] / 1e-9:.2f} +/- {r['D_sem'] / 1e-9:.2f}")
    plt.axvspan(*TIMESTEP_FIT_PS, color="0.92", zorder=0)
    plt.xlabel("Time (ps)")
    plt.ylabel(r"Oxygen MSD ($\AA^2$)")
    plt.title(f"Oxygen MSD vs time step (fit {TIMESTEP_FIT_PS[0]:g}-{TIMESTEP_FIT_PS[1]:g} ps, "
              r"D in $10^{-9}$ m$^2$/s)", fontsize=10)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(OUT / "timestep_msd.png", dpi=220)
    plt.close(fig)

    with open(OUT / "timestep_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dt_fs", "mean_T_K", "std_T_K", "mean_PE_kcal_mol", "PE_block_sem", "std_PE",
                    "mean_KE_kcal_mol", "std_KE", "mean_Etot_kcal_mol", "std_Etot",
                    "Etot_drift_kcal_mol_per_ps",
                    "prod_wall_s", "wall_s_per_ps", "D_m2_s", "D_window_sem_m2_s"])
        for r, sem, wps in zip(results, pe_sem, wall_per_ps):
            th = r["thermo"]
            w.writerow([r["dt"], f"{th['Temp'].mean():.2f}", f"{th['Temp'].std():.2f}",
                        f"{th['PotEng'].mean():.2f}", f"{sem:.2f}", f"{th['PotEng'].std():.1f}",
                        f"{th['KinEng'].mean():.2f}", f"{th['KinEng'].std():.1f}",
                        f"{th['TotEng'].mean():.2f}", f"{th['TotEng'].std():.1f}",
                        f"{drift_per_ps(th['t_ps'], th['TotEng']):.3f}",
                        f"{r['loop_time']:.1f}", f"{wps:.2f}",
                        f"{r['D']:.4e}", f"{r['D_sem']:.2e}"])

    return results


# ------------------------------------------------------------
# Part 2: ensemble comparison
# ------------------------------------------------------------

def ensemble_comparison():
    results = []
    for ens in ENSEMBLES:
        folder = HERE / f"week8_results_2sec-{ens}"
        thermo, loop_time = read_production_thermo(folder / "log.lammps")
        t_msd, msd = read_msd(find_ensemble_msd(folder, ens), ENSEMBLE_DT_FS)
        d, slope, intercept = fit_diffusion(t_msd, msd, ENSEMBLE_FIT_PS)
        results.append({
            "ens": ens, "thermo": thermo, "loop_time": loop_time,
            "t_msd": t_msd, "msd": msd, "D": d, "slope": slope, "intercept": intercept,
            "D_sem": window_sem(t_msd, msd, ENSEMBLE_FIT_PS),
        })

    # MSD + fits
    fig = plt.figure(figsize=(7.4, 4.8))
    for r in results:
        line, = plt.plot(r["t_msd"], r["msd"], lw=1.2,
                         label=f"{r['ens']}: D = {r['D'] / 1e-9:.2f} +/- {r['D_sem'] / 1e-9:.2f}")
        tf = np.linspace(*ENSEMBLE_FIT_PS, 50)
        plt.plot(tf, r["slope"] * tf + r["intercept"], "--", lw=1.0, color=line.get_color())
    plt.axvspan(*ENSEMBLE_FIT_PS, color="0.94", zorder=0)
    plt.xlabel("Time (ps)")
    plt.ylabel(r"Oxygen MSD ($\AA^2$)")
    plt.title(f"Oxygen MSD by ensemble, 2 fs (fit {ENSEMBLE_FIT_PS[0]:g}-{ENSEMBLE_FIT_PS[1]:g} ps, "
              r"D in $10^{-9}$ m$^2$/s)", fontsize=10)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(OUT / "ensemble_msd.png", dpi=220)
    plt.close(fig)

    # T, density, Etot vs time
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    for r in results:
        th, ens = r["thermo"], r["ens"]
        axes[0].plot(th["t_ps"], th["Temp"], lw=0.9,
                     label=f"{ens} (mean {th['Temp'].mean():.1f} K)")
        axes[1].plot(th["t_ps"], th["Density"], lw=0.9,
                     label=f"{ens} (mean {th['Density'].mean():.4f})")
        axes[2].plot(th["t_ps"], th["TotEng"], lw=0.9,
                     label=f"{ens} (std {th['TotEng'].std():.2f})")
    axes[0].axhline(300.0, color="k", lw=0.6, ls=":")
    axes[0].set_ylabel("Temperature (K)")
    axes[1].set_ylabel(r"Density (g/cm$^3$)")
    axes[2].set_ylabel("Total energy\n(kcal/mol)")
    axes[2].set_xlabel("Production time (ps)")
    axes[0].set_title("Ensemble comparison (2 fs, 50 ps production)")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "ensemble_thermo.png", dpi=220)
    plt.close(fig)

    with open(OUT / "ensemble_summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ensemble", "mean_T_K", "std_T_K", "mean_density_g_cm3", "std_density",
                    "mean_press_atm", "mean_Etot_kcal_mol", "std_Etot", "Etot_drift_kcal_mol_per_ps",
                    "fit_start_ps", "fit_end_ps", "D_m2_s", "D_window_sem_m2_s"])
        for r in results:
            th = r["thermo"]
            w.writerow([r["ens"], f"{th['Temp'].mean():.2f}", f"{th['Temp'].std():.2f}",
                        f"{th['Density'].mean():.4f}", f"{th['Density'].std():.4f}",
                        f"{th['Press'].mean():.0f}",
                        f"{th['TotEng'].mean():.2f}", f"{th['TotEng'].std():.2f}",
                        f"{drift_per_ps(th['t_ps'], th['TotEng']):.3f}",
                        ENSEMBLE_FIT_PS[0], ENSEMBLE_FIT_PS[1],
                        f"{r['D']:.4e}", f"{r['D_sem']:.2e}"])

    return results


def main():
    OUT.mkdir(exist_ok=True)
    timestep_comparison()
    ensemble_comparison()
    for name in ("timestep_summary.csv", "ensemble_summary.csv"):
        print(f"\n{name}")
        print((OUT / name).read_text())
    print(f"Figures written to {OUT}")


if __name__ == "__main__":
    main()
