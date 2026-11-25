#!/usr/bin/env python3
"""
Analyze LAMMPS PE, KE, and total energy versus MD step.
Also perform block/bin-size analysis and estimate correlation-aware uncertainty.

Expected thermo output, e.g.:
    thermo_style custom step time temp pe ke etotal press density

Examples:
    python analyze_lammps_energy.py log.lammps
    python analyze_lammps_energy.py log.lammps --start-step 25000

Outputs:
    energy_vs_step.png
    blocking_analysis.png
    energy_statistics.csv
"""

import argparse
import csv
import math
from statistics import NormalDist

import numpy as np
import matplotlib.pyplot as plt

ALIASES = {
    "Step": {"Step", "step"},
    "Time": {"Time", "time"},
    "Temp": {"Temp", "temp"},
    "PE": {"PotEng", "PE", "pe"},
    "KE": {"KinEng", "KE", "ke"},
    "Etot": {"TotEng", "Etot", "ETotal", "etotal"},
}


def canonical_name(token):
    for key, names in ALIASES.items():
        if token in names:
            return key
    return None


def read_lammps_log(path):
    """Read all LAMMPS thermo blocks containing Step, PE, KE, and Etot."""
    rows = []
    header = None
    idx = None

    with open(path, "r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            parts = s.split()
            mapped = [canonical_name(x) for x in parts]

            if "Step" in mapped and "PE" in mapped and "KE" in mapped and "Etot" in mapped:
                header = parts
                idx = {}
                for i, tok in enumerate(parts):
                    name = canonical_name(tok)
                    if name is not None and name not in idx:
                        idx[name] = i
                continue

            if header is None:
                continue

            if s.startswith(("Loop time", "ERROR", "WARNING")):
                header = None
                idx = None
                continue

            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue

            if len(vals) != len(header):
                continue

            rows.append({
                "Step": vals[idx["Step"]],
                "Time": vals[idx["Time"]] if "Time" in idx else np.nan,
                "Temp": vals[idx["Temp"]] if "Temp" in idx else np.nan,
                "PE": vals[idx["PE"]],
                "KE": vals[idx["KE"]],
                "Etot": vals[idx["Etot"]],
            })

    if not rows:
        raise RuntimeError(
            "No thermo table found. Expected Step, PotEng, KinEng, and TotEng columns."
        )

    return {
        key: np.asarray([r[key] for r in rows], dtype=float)
        for key in ("Step", "Time", "Temp", "PE", "KE", "Etot")
    }


def autocorrelation_fft(x):
    """Normalized autocorrelation function computed with FFT."""
    x = np.asarray(x, dtype=float)
    x = x - np.mean(x)
    n = len(x)
    nfft = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x, n=nfft)
    acf = np.fft.irfft(f * np.conjugate(f), n=nfft)[:n]
    acf /= np.arange(n, 0, -1)
    if acf[0] == 0:
        return np.ones(n)
    return acf / acf[0]


def integrated_autocorrelation_time(x):
    """
    Estimate integrated autocorrelation time in stored-sample units.
    Sum the ACF until its first non-positive value.
    """
    acf = autocorrelation_fft(x)
    total = 0.0
    for rho in acf[1:]:
        if not np.isfinite(rho) or rho <= 0.0:
            break
        total += rho
    return max(0.5, 0.5 + total)


def block_sem(x, block_size):
    """SEM calculated from contiguous, non-overlapping block means."""
    x = np.asarray(x, dtype=float)
    b = int(block_size)
    nblocks = len(x) // b
    if nblocks < 2:
        return np.nan, nblocks

    y = x[:nblocks * b].reshape(nblocks, b).mean(axis=1)
    return y.std(ddof=1) / math.sqrt(nblocks), nblocks


def blocking_curve(x, min_blocks=8):
    sizes, sems, nblocks = [], [], []
    b = 1
    while len(x) // b >= min_blocks:
        sem, nb = block_sem(x, b)
        sizes.append(b)
        sems.append(sem)
        nblocks.append(nb)
        b *= 2
    return np.asarray(sizes), np.asarray(sems), np.asarray(nblocks)


def choose_block_size(x, tau_int, min_blocks=10):
    """Choose a power-of-two block length of about 2*tau_int samples."""
    desired = max(1, int(math.ceil(2.0 * tau_int)))
    b = 1
    while b < desired:
        b *= 2

    max_b = max(1, len(x) // min_blocks)
    if b > max_b:
        b = 1
        while 2 * b <= max_b:
            b *= 2
    return max(1, b)


def summarize(x):
    x = np.asarray(x, dtype=float)
    n = len(x)
    mean = x.mean()
    std = x.std(ddof=1)
    naive_sem = std / math.sqrt(n)

    tau = integrated_autocorrelation_time(x)
    g = max(1.0, 2.0 * tau)
    neff = n / g
    corr_sem = std / math.sqrt(neff)

    b = choose_block_size(x, tau)
    bsem, nblocks = block_sem(x, b)
    z95 = NormalDist().inv_cdf(0.975)
    ci95 = z95 * bsem

    return {
        "N": n,
        "mean": mean,
        "std": std,
        "naive_sem": naive_sem,
        "tau_int_samples": tau,
        "stat_inefficiency_g": g,
        "N_eff": neff,
        "corr_sem": corr_sem,
        "block_size": b,
        "n_blocks": nblocks,
        "block_sem": bsem,
        "CI95_halfwidth": ci95,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("logfile", help="LAMMPS log file, e.g. log.lammps")
    parser.add_argument("--start-step", type=float, default=None,
                        help="discard data before this MD step")
    parser.add_argument("--end-step", type=float, default=None,
                        help="discard data after this MD step")
    args = parser.parse_args()

    data = read_lammps_log(args.logfile)
    mask = np.ones(len(data["Step"]), dtype=bool)

    if args.start_step is not None:
        mask &= data["Step"] >= args.start_step
    if args.end_step is not None:
        mask &= data["Step"] <= args.end_step

    for key in data:
        data[key] = data[key][mask]

    if len(data["Step"]) < 20:
        raise RuntimeError("Too few thermo samples remain after filtering.")

    step = data["Step"]
    series = {
        "Temperature": data["Temp"],
        "Potential energy": data["PE"],
        "Kinetic energy": data["KE"],
        "Total energy": data["Etot"],
    }

    # Plot energy versus MD step.
    
    for name, values in series.items():
        fig = plt.figure(figsize=(8.2, 5.0))
        plt.plot(step, values, linewidth=1.0, label=name)
        plt.xlabel("MD step")
        plt.ylabel(name + " (kcal/mol)" if name != "Temperature" else "Temperature (K)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(name + "vs_step.png", dpi=220)
        plt.close(fig)

    # Statistical analysis.
    stats = {name: summarize(values) for name, values in series.items()}

    # Plot estimated SEM as a function of block/bin size.
    fig = plt.figure(figsize=(8.2, 5.0))
    for name, values in series.items():
        sizes, sems, _ = blocking_curve(values)
        plt.plot(sizes, sems, marker="o", label=name)
    plt.xscale("log", base=2)
    plt.xlabel("Block/bin size (stored thermo samples)")
    plt.ylabel("Estimated SEM of the mean (kcal/mol)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("blocking_analysis.png", dpi=220)
    plt.close(fig)

    fields = [
        "quantity", "N", "mean", "std", "naive_sem",
        "tau_int_samples", "stat_inefficiency_g", "N_eff",
        "corr_sem", "block_size", "n_blocks", "block_sem",
        "CI95_halfwidth"
    ]

    with open("energy_statistics.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for name, s in stats.items():
            row = {"quantity": name}
            row.update(s)
            writer.writerow(row)

    print("\nEnergy statistics")
    print("=" * 72)
    print(f"Samples analyzed: {len(step)}")
    print(f"Step range: {step[0]:.0f} to {step[-1]:.0f}\n")

    for name, s in stats.items():
        print(name)
        print("-" * len(name))
        print(f"Mean                  = {s['mean']:.8g} kcal/mol")
        print(f"Raw std. deviation    = {s['std']:.8g} kcal/mol")
        print(f"Naive SEM             = {s['naive_sem']:.8g} kcal/mol")
        print(f"tau_int               = {s['tau_int_samples']:.3f} stored samples")
        print(f"Effective sample size = {s['N_eff']:.1f} of {s['N']}")
        print(f"Correlation-aware SEM = {s['corr_sem']:.8g} kcal/mol")
        print(f"Chosen block size     = {s['block_size']} samples")
        print(f"Number of blocks      = {s['n_blocks']}")
        print(f"Block SEM             = {s['block_sem']:.8g} kcal/mol")
        print(f"Mean +/- 95% CI       = {s['mean']:.8g} +/- {s['CI95_halfwidth']:.8g} kcal/mol")
        print()

    print("Files written:")
    print("  energy_vs_step.png")
    print("  blocking_analysis.png")
    print("  energy_statistics.csv")
    print()
    print("Blocking interpretation:")
    print("  MD samples are correlated, so std/sqrt(N) usually underestimates")
    print("  the uncertainty of the mean. As block size increases, the block")
    print("  SEM should rise and then approach a plateau once a block is longer")
    print("  than the relevant correlation time.")


if __name__ == "__main__":
    main()
