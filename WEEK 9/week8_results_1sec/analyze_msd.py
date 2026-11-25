#!/usr/bin/env python3
"""
Analyze oxygen mean-squared displacement (MSD) from a LAMMPS water simulation.

This script is designed for the output produced by:

    group oxygen type 1
    compute msdO oxygen msd com yes
    fix msdout all ave/time 100 1 100 c_msdO[4] file msd_O.dat

The file is expected to contain:
    timestep   MSD(A^2)

The diffusion coefficient is obtained from the Einstein relation:

    MSD(t) = 6 D t

for three-dimensional diffusion.

Example
-------
For a 2 fs timestep and a fit from 50 to 500 ps:

    python analyze_msd.py msd_O.dat --dt-fs 2.0 \
        --fit-start 50 --fit-end 500

Outputs
-------
msd_vs_time.png
msd_fit.png
diffusion_summary.txt
diffusion_windows.csv

Important statistical note
--------------------------
MSD values at neighboring times are strongly correlated. Therefore, the
ordinary linear-regression standard error is NOT a reliable uncertainty in D.

This script estimates a more useful uncertainty by repeating the diffusion
fit over several non-overlapping time windows within the selected diffusive
regime. The spread among those window-based D values is reported as a
"window-to-window" uncertainty.

For rigorous transport-property uncertainties, use a trajectory analysis with
multiple time origins and/or independent trajectories.
"""

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def load_msd(path):
    """Read numeric lines from a LAMMPS fix ave/time output."""
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                pass

    if not rows:
        raise RuntimeError(f"No numeric MSD data found in {path}")

    arr = np.asarray(rows, dtype=float)
    return arr[:, 0], arr[:, 1]


def linear_fit(t, msd):
    """Return slope, intercept, R^2."""
    slope, intercept = np.polyfit(t, msd, 1)
    pred = slope * t + intercept

    ss_res = np.sum((msd - pred) ** 2)
    ss_tot = np.sum((msd - np.mean(msd)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return slope, intercept, r2


def diffusion_from_slope(slope_A2_per_ps):
    """
    MSD = 6Dt in 3D.

    D in A^2/ps:
        D = slope/6

    Conversion:
        1 A^2/ps = 1e-8 m^2/s
    """
    d_a2_ps = slope_A2_per_ps / 6.0
    d_m2_s = d_a2_ps * 1.0e-8
    return d_a2_ps, d_m2_s


def window_uncertainty(t, msd, fit_start, fit_end, n_windows):
    """
    Divide the selected fit interval into non-overlapping windows and fit
    each window independently.

    This is intended as a pedagogical estimate of sensitivity to fit window,
    not as a rigorous statistical confidence interval.
    """
    edges = np.linspace(fit_start, fit_end, n_windows + 1)
    results = []

    for i in range(n_windows):
        lo, hi = edges[i], edges[i + 1]

        if i == n_windows - 1:
            mask = (t >= lo) & (t <= hi)
        else:
            mask = (t >= lo) & (t < hi)

        if mask.sum() < 3:
            continue

        slope, intercept, r2 = linear_fit(t[mask], msd[mask])
        d_a2_ps, d_m2_s = diffusion_from_slope(slope)

        results.append({
            "window": i + 1,
            "start_ps": lo,
            "end_ps": hi,
            "n_points": int(mask.sum()),
            "slope_A2_per_ps": slope,
            "D_A2_per_ps": d_a2_ps,
            "D_m2_per_s": d_m2_s,
            "R2": r2,
        })

    return results


def main():
    p = argparse.ArgumentParser(
        description="Fit water MSD and estimate the self-diffusion coefficient."
    )
    p.add_argument("msd_file", help="LAMMPS MSD file, e.g. msd_O.dat")
    p.add_argument(
        "--dt-fs", type=float, default=2.0,
        help="LAMMPS MD timestep in fs (default: 2.0)"
    )
    p.add_argument(
        "--fit-start", type=float, required=True,
        help="start of the linear MSD fitting region in ps"
    )
    p.add_argument(
        "--fit-end", type=float, required=True,
        help="end of the linear MSD fitting region in ps"
    )
    p.add_argument(
        "--windows", type=int, default=5,
        help="number of non-overlapping subwindows used to estimate fit-window uncertainty (default: 5)"
    )
    p.add_argument(
        "--experimental-D", type=float, default=2.299e-9,
        help="reference diffusion coefficient in m^2/s (default: 2.299e-9)"
    )
    args = p.parse_args()

    step, msd = load_msd(args.msd_file)

    # Convert step number to physical time.
    time_ps = step * args.dt_fs / 1000.0

    # Shift time so the first stored MSD point is t=0 if necessary.
    # This is useful if the production run started after reset_timestep 0
    # but the first saved point is not exactly zero.
    time_ps = time_ps - time_ps[0]

    # ----------------------------------------------------------
    # Overall MSD plot
    # ----------------------------------------------------------
    fig = plt.figure(figsize=(7.4, 4.8))
    plt.plot(time_ps, msd, linewidth=1.2)
    plt.xlabel("Time (ps)")
    plt.ylabel(r"Oxygen MSD ($\AA^2$)")
    plt.tight_layout()
    plt.savefig("msd_vs_time.png", dpi=220)
    plt.close(fig)

    # ----------------------------------------------------------
    # Main diffusion fit
    # ----------------------------------------------------------
    mask = (time_ps >= args.fit_start) & (time_ps <= args.fit_end)

    if mask.sum() < 3:
        raise RuntimeError(
            "The requested fit interval contains fewer than three MSD points."
        )

    tfit = time_ps[mask]
    mfit = msd[mask]

    slope, intercept, r2 = linear_fit(tfit, mfit)
    d_a2_ps, d_m2_s = diffusion_from_slope(slope)

    # ----------------------------------------------------------
    # Window-to-window analysis
    # ----------------------------------------------------------
    window_results = window_uncertainty(
        time_ps, msd, args.fit_start, args.fit_end, args.windows
    )

    if len(window_results) >= 2:
        d_windows = np.array([x["D_m2_per_s"] for x in window_results])
        d_window_mean = np.mean(d_windows)
        d_window_std = np.std(d_windows, ddof=1)
        d_window_sem = d_window_std / math.sqrt(len(d_windows))
    else:
        d_window_mean = np.nan
        d_window_std = np.nan
        d_window_sem = np.nan

    # ----------------------------------------------------------
    # Plot MSD with fitted line
    # ----------------------------------------------------------
    fig = plt.figure(figsize=(7.4, 4.8))
    plt.plot(time_ps, msd, label="Oxygen MSD", linewidth=1.2)
    plt.plot(
        tfit,
        slope * tfit + intercept,
        "--",
        linewidth=1.7,
        label=f"Linear fit: {args.fit_start:g}-{args.fit_end:g} ps"
    )
    plt.axvline(args.fit_start, linestyle=":", linewidth=1.0)
    plt.axvline(args.fit_end, linestyle=":", linewidth=1.0)
    plt.xlabel("Time (ps)")
    plt.ylabel(r"Oxygen MSD ($\AA^2$)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("msd_fit.png", dpi=220)
    plt.close(fig)

    # ----------------------------------------------------------
    # Save subwindow results
    # ----------------------------------------------------------
    with open("diffusion_windows.csv", "w", newline="") as f:
        fields = [
            "window", "start_ps", "end_ps", "n_points",
            "slope_A2_per_ps", "D_A2_per_ps", "D_m2_per_s", "R2"
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(window_results)

    # ----------------------------------------------------------
    # Comparison with reference experiment
    # ----------------------------------------------------------
    ref = args.experimental_D
    percent_error = 100.0 * (d_m2_s - ref) / ref
    abs_percent_error = abs(percent_error)

    lines = []
    lines.append("MSD diffusion analysis")
    lines.append("=" * 60)
    lines.append(f"Input file                  : {args.msd_file}")
    lines.append(f"MD timestep                 : {args.dt_fs:g} fs")
    lines.append(f"Fit interval                : {args.fit_start:g} to {args.fit_end:g} ps")
    lines.append(f"Number of fit points        : {len(tfit)}")
    lines.append("")
    lines.append(f"MSD slope                   : {slope:.8g} A^2/ps")
    lines.append(f"Fit R^2                     : {r2:.6f}")
    lines.append(f"D                           : {d_a2_ps:.8g} A^2/ps")
    lines.append(f"D                           : {d_m2_s:.8g} m^2/s")
    lines.append(f"D                           : {d_m2_s/1e-9:.5f} x 10^-9 m^2/s")
    lines.append("")
    lines.append(f"Reference D                 : {ref:.8g} m^2/s")
    lines.append(f"Signed percent difference   : {percent_error:.3f} %")
    lines.append(f"Absolute percent difference : {abs_percent_error:.3f} %")
    lines.append("")
    lines.append("Window/bin sensitivity analysis")
    lines.append("-" * 60)
    lines.append(f"Requested windows           : {args.windows}")
    lines.append(f"Usable windows              : {len(window_results)}")

    if len(window_results) >= 2:
        lines.append(
            f"Mean D across windows       : {d_window_mean/1e-9:.5f} x 10^-9 m^2/s"
        )
        lines.append(
            f"Std. dev. across windows    : {d_window_std/1e-9:.5f} x 10^-9 m^2/s"
        )
        lines.append(
            f"SEM across windows          : {d_window_sem/1e-9:.5f} x 10^-9 m^2/s"
        )

    lines.append("")
    lines.append("Interpretation")
    lines.append("-" * 60)
    lines.append(
        "Use only a time region where MSD is approximately linear. "
        "Do not fit the very early ballistic regime."
    )
    lines.append(
        "Changing the fit interval is a useful way to see whether the reported "
        "diffusion coefficient is stable."
    )
    lines.append(
        "The subwindow spread is a pedagogical uncertainty estimate. "
        "MSD points are strongly correlated, so ordinary least-squares fit "
        "errors should not be interpreted as rigorous statistical uncertainty."
    )
    lines.append(
        "For publication-quality uncertainty, use multiple time origins and/or "
        "several statistically independent trajectories."
    )

    summary = "\n".join(lines)

    with open("diffusion_summary.txt", "w") as f:
        f.write(summary + "\n")

    print(summary)
    print("\nCreated:")
    print("  msd_vs_time.png")
    print("  msd_fit.png")
    print("  diffusion_windows.csv")
    print("  diffusion_summary.txt")


if __name__ == "__main__":
    main()
