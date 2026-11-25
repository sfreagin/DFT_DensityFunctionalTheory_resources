#!/usr/bin/env python3
"""Recompute Project 3 diffusion and AIMD equilibration numbers.

Reads only:
  WEEK 11/Project3/Project3/week8_results_2sec-NVT*   (msd_NVT.dat, log.lammps)
  WEEK 11/Project3/Project3/week10_seed*             (msd_O.dat, log.lammps)
  WEEK 11/week11/week11 and week11-round2..round8    (qe_thermo.txt, 32h2o_qe.out)

Writes figures and project3_analysis.json in this directory.
"""

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/Users/stephenreagin/Desktop/SJSU/MATE 250/WEEK 11")
PAPER = Path("/Users/stephenreagin/Desktop/SJSU/MATE 250/WEEK 11/Paper")
CLASSICAL = [
    ("seed 1", ROOT / "Project3/Project3/week8_results_2sec-NVT_seed1"),
    ("seed 2", ROOT / "Project3/Project3/week8_results_2sec-NVT-seed2"),
    ("seed 3", ROOT / "Project3/Project3/week8_results_2sec-NVT-seed3"),
]
DEEPMD = [
    ("seed 1", ROOT / "Project3/Project3/week10_seed1"),
    ("seed 2", ROOT / "Project3/Project3/week10_seed2"),
    ("seed 3", ROOT / "Project3/Project3/week10_seed3"),
]
AIMD_ROOT = ROOT / "week11"
AIMD_ROUNDS = ["week11"] + [f"week11-round{i}" for i in range(2, 9)]
BOX_AIMD = 9.841886522
DT_AIMD_FS = 0.5
EXP_D = 2.299e-9


def load_msd(path):
    rows = []
    with open(path, errors="ignore") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 2 or parts[0].startswith("#"):
                continue
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    arr = np.asarray(rows, dtype=float)
    return arr[:, 0], arr[:, 1]


def linear_fit(t, y):
    slope, intercept = np.polyfit(t, y, 1)
    pred = slope * t + intercept
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(slope), float(intercept), float(r2)


def d_from_slope(slope_a2_per_ps):
    return slope_a2_per_ps / 6.0 * 1.0e-8


def fit_window(t, msd, start, end):
    mask = (t >= start) & (t <= end)
    slope, intercept, r2 = linear_fit(t[mask], msd[mask])
    return {
        "n_points": int(np.count_nonzero(mask)),
        "slope_A2_per_ps": slope,
        "intercept_A2": intercept,
        "R2": r2,
        "D_m2_s": d_from_slope(slope),
        "final_msd_A2": float(msd[-1]),
        "t_final_ps": float(t[-1]),
    }


def window_spread(t, msd, start, end, n_windows):
    edges = np.linspace(start, end, n_windows + 1)
    values = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (t >= lo) & (t <= hi)
        if np.count_nonzero(mask) < 5:
            continue
        slope, _, r2 = linear_fit(t[mask], msd[mask])
        values.append({"lo": float(lo), "hi": float(hi), "D_m2_s": d_from_slope(slope), "R2": r2})
    ds = np.array([row["D_m2_s"] for row in values])
    return values, float(ds.std(ddof=1)) if len(ds) > 1 else None


def parse_lammps_log(path):
    """Return thermo rows for each run section and the wall-time string."""
    header = None
    sections = []
    rows = []
    wall = None
    loop = None
    with open(path, errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("Step ") and "Temp" in stripped:
                if rows:
                    sections.append((header, np.asarray(rows, dtype=float)))
                    rows = []
                header = stripped.split()
                continue
            if line.startswith("Loop time"):
                loop = line.strip()
                if rows:
                    sections.append((header, np.asarray(rows, dtype=float)))
                    rows = []
                header = None
                continue
            if line.startswith("Total wall time:"):
                wall = line.split(":", 1)[1].strip()
                continue
            if header is None:
                continue
            parts = line.split()
            if len(parts) != len(header):
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    if rows:
        sections.append((header, np.asarray(rows, dtype=float)))
    return sections, wall, loop


def column_stats(sections, name, section_index=-1):
    header, data = sections[section_index]
    idx = header.index(name)
    values = data[:, idx]
    return {
        "mean": float(values.mean()),
        "std": float(values.std(ddof=1)),
        "first": float(values[0]),
        "last": float(values[-1]),
        "n": int(len(values)),
    }


def lammps_thermo(sections, time_scale):
    """Temperature and energies versus time, in ps. Concatenate every thermo block."""
    chunks = []
    for header, data in sections:
        chunks.append({
            "t_ps": data[:, header.index("Time")] * time_scale,
            "temp": data[:, header.index("Temp")],
            "pe": data[:, header.index("PotEng")],
            "ke": data[:, header.index("KinEng")],
            "etotal": data[:, header.index("TotEng")],
        })
    if len(chunks) == 1:
        return chunks[0]
    return {key: np.concatenate([chunk[key] for chunk in chunks]) for key in chunks[0]}


def analyze_classical():
    results = []
    curves = []
    thermo = []
    for label, folder in CLASSICAL:
        steps, msd = load_msd(folder / "msd_NVT.dat")
        t = steps * 2.0 / 1000.0
        fit = fit_window(t, msd, 50.0, 200.0)
        windows, window_std = window_spread(t, msd, 50.0, 200.0, 5)
        sections, wall, loop = parse_lammps_log(folder / "log.lammps")
        temp = column_stats(sections, "Temp")
        density = column_stats(sections, "Density")
        results.append({
            "label": label,
            "folder": str(folder.relative_to(ROOT)),
            "wall_time": wall,
            "loop": loop,
            "temperature_K": temp,
            "density_g_cm3": density,
            "fit_50_200ps": fit,
            "window_D_std_m2_s": window_std,
            "windows": windows,
        })
        curves.append((label, t, msd, fit))
        # Real-unit time is printed in fs for these 2 fs runs. MSD time uses the same conversion.
        thermo.append((label, lammps_thermo(sections, 0.001)))
    return results, curves, thermo


def analyze_deepmd():
    results = []
    curves = []
    thermo = []
    for label, folder in DEEPMD:
        steps, msd = load_msd(folder / "msd_O.dat")
        t = (steps - steps[0]) * 0.0005
        fit = fit_window(t, msd, 50.0, 200.0)
        windows, window_std = window_spread(t, msd, 50.0, 200.0, 5)
        sections, wall, loop = parse_lammps_log(folder / "log.lammps")
        # First run is the 10 ps equilibration; last run is 200 ps production.
        eq = sections[0]
        prod = sections[-1]
        eq_temp = {
            "mean": float(eq[1][:, eq[0].index("Temp")].mean()),
            "std": float(eq[1][:, eq[0].index("Temp")].std(ddof=1)),
        }
        temp = column_stats(sections, "Temp", -1)
        density = column_stats(sections, "Density", -1)
        results.append({
            "label": label,
            "folder": str(folder.relative_to(ROOT)),
            "wall_time": wall,
            "production_loop": loop,
            "equilibration_temperature_K": eq_temp,
            "temperature_K": temp,
            "density_g_cm3": density,
            "fit_50_200ps": fit,
            "window_D_std_m2_s": window_std,
            "windows": windows,
            "msd_step0": float(steps[0]),
            "msd_step_final": float(steps[-1]),
            "n_msd": int(len(steps)),
        })
        curves.append((label, t, msd, fit))
        # Metal units: the Time column is already in ps and continues through the 10 ps warmup.
        thermo.append((label, lammps_thermo(sections, 1.0)))
    return results, curves, thermo


STEP_RE = re.compile(r"Entering\s+Dynamics\s*:\s*iteration\s*=\s*(\d+)", re.I)
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
EKIN_RE = re.compile(r"kinetic\s+energy\s*\(Ekin\)\s*=\s*(" + NUMBER + r")", re.I)
TEMP_RE = re.compile(r"^\s*temperature\s*=\s*(" + NUMBER + r")\s*K\b", re.I)
SUM_RE = re.compile(r"Ekin\s*\+\s*Etot\s*\(const\)\s*=\s*(" + NUMBER + r")", re.I)
NOSE_RE = re.compile(r"Ekin\s*\+\s*Etot\s*\+\s*Ions\s+Nose\s*=\s*(" + NUMBER + r")", re.I)
POS_RE = re.compile(r"^\s*ATOMIC_POSITIONS\s*(?:\(\s*([^\)]+)\s*\))?", re.I)


def fnum(text):
    return float(text.replace("D", "E").replace("d", "e"))


def parse_thermo(path):
    records = []
    step = None
    pending = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        match = STEP_RE.search(line)
        if match:
            step = int(match.group(1))
        match = EKIN_RE.search(line)
        if match:
            pending["ekin"] = fnum(match.group(1))
        match = TEMP_RE.search(line)
        if match:
            pending["temp"] = fnum(match.group(1))
        match = SUM_RE.search(line)
        if match:
            pending["sum"] = fnum(match.group(1))
        match = NOSE_RE.search(line)
        if match:
            pending["nose"] = fnum(match.group(1))
            if step is not None and all(k in pending for k in ("ekin", "temp", "sum", "nose")):
                ekin = pending["ekin"]
                records.append((step, pending["temp"], pending["sum"] - ekin, ekin, pending["sum"], pending["nose"]))
            pending = {}
    return records


def parse_oxygen(path):
    """Oxygen coordinates printed after each 'Entering Dynamics' line."""
    frames = {}
    step = None
    expect = False
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    i = 0
    while i < len(lines):
        match = STEP_RE.search(lines[i])
        if match:
            step = int(match.group(1))
            expect = True
            i += 1
            continue
        match = POS_RE.match(lines[i])
        if match and expect and step is not None:
            unit = (match.group(1) or "").strip().lower()
            oxygen = []
            j = i + 1
            while j < len(lines) and len(oxygen) < 32:
                parts = lines[j].split()
                if len(parts) < 4:
                    break
                if re.fullmatch(r"O\d*", parts[0], re.I):
                    oxygen.append([fnum(v) for v in parts[1:4]])
                j += 1
            if len(oxygen) == 32 and unit.startswith("angstrom"):
                frames[step] = oxygen
            expect = False
            i = j
            continue
        i += 1
    return frames


def unwrap(coords, box):
    wrapped = np.mod(coords, box)
    jumps = (wrapped[1:] - wrapped[:-1]) / box
    jumps -= np.rint(jumps)
    out = np.empty_like(wrapped)
    out[0] = wrapped[0]
    out[1:] = wrapped[0] + np.cumsum(jumps * box, axis=0)
    com = out.mean(axis=1, keepdims=True)
    return out - com


def multi_origin_msd(unwrapped, max_lag):
    n, n_ox, _ = unwrapped.shape
    coords = unwrapped.reshape(n, -1)
    sq = np.einsum("ij,ij->i", coords, coords)
    prefix = np.r_[0.0, np.cumsum(sq)]
    lag = np.arange(max_lag + 1)
    sum_sq = prefix[n] - prefix[lag] + prefix[n - lag]
    acf = np.zeros(max_lag + 1)
    fft_n = 2 * n
    for col in range(coords.shape[1]):
        ft = np.fft.rfft(coords[:, col], n=fft_n)
        acf += np.fft.irfft(ft * np.conj(ft), n=fft_n)[: max_lag + 1]
    msd = (sum_sq - 2.0 * acf) / ((n - lag) * n_ox)
    msd[0] = 0.0
    return np.maximum(msd, 0.0)


def analyze_aimd():
    thermo = []
    frames = {}
    rounds = []
    for name in AIMD_ROUNDS:
        folder = AIMD_ROOT / name
        rec = parse_thermo(folder / "qe_thermo.txt")
        pos = parse_oxygen(folder / "32h2o_qe.out")
        rounds.append({
            "name": name,
            "thermo_first": rec[0][0],
            "thermo_last": rec[-1][0],
            "n_thermo": len(rec),
            "n_position_frames": len(pos),
            "pos_first": min(pos) if pos else None,
            "pos_last": max(pos) if pos else None,
        })
        thermo.extend(rec)
        frames.update(pos)
        print(f"  {name}: thermo {rec[0][0]}-{rec[-1][0]} ({len(rec)}), positions {len(pos)}")

    thermo = sorted(set(thermo))
    steps = np.array([row[0] for row in thermo], dtype=int)
    if len(steps) != len(set(steps.tolist())) or np.any(np.diff(steps) != 1):
        raise RuntimeError("AIMD thermo iterations are not a continuous unit sequence")
    temp = np.array([row[1] for row in thermo])
    etot = np.array([row[2] for row in thermo])
    ekin = np.array([row[3] for row in thermo])
    conserved = np.array([row[4] for row in thermo])
    nose = np.array([row[5] for row in thermo])
    t_ps = (steps - 1) * DT_AIMD_FS / 1000.0

    order = np.array(sorted(frames))
    if not np.array_equal(order, steps):
        missing = np.setdiff1d(steps, order)
        extra = np.setdiff1d(order, steps)
        raise RuntimeError(f"Position frames do not match thermo steps. missing={missing[:8]} extra={extra[:8]}")
    coords = np.asarray([frames[int(s)] for s in order], dtype=float)
    unwrapped = unwrap(coords, BOX_AIMD)
    # Subsample every 10 steps (5 fs) for the all-origin MSD. Oxygen does not
    # cross half the 9.84 Å box in 5 fs, and 6.9 ps at every step is unnecessary.
    stride = 10
    unwrapped = unwrapped[::stride]
    frame_ps = DT_AIMD_FS * stride / 1000.0
    max_lag = len(unwrapped) // 2
    msd = multi_origin_msd(unwrapped, max_lag)
    lag_ps = np.arange(max_lag + 1) * frame_ps
    single = np.mean(np.sum((unwrapped - unwrapped[0]) ** 2, axis=2), axis=1)
    single_ps = np.arange(len(unwrapped)) * frame_ps

    windows = []
    for lo, hi in [(0.2, 1.0), (0.5, 1.5), (1.0, 2.0), (1.5, 2.5), (1.0, 3.0), (2.0, 3.4)]:
        if hi > lag_ps[-1]:
            continue
        windows.append({"lo": lo, "hi": hi, **fit_window(lag_ps, msd, lo, hi)})
    ds = np.array([row["D_m2_s"] for row in windows])

    def segment(mask):
        return {
            "mean_T_K": float(temp[mask].mean()),
            "std_T_K": float(temp[mask].std(ddof=1)),
            "mean_Etot_Ry": float(etot[mask].mean()),
            "std_Etot_Ry": float(etot[mask].std(ddof=1)),
            "mean_nose_Ry": float(nose[mask].mean()),
            "drift_nose_mRy_per_ps": float(np.polyfit(t_ps[mask], nose[mask], 1)[0] * 1000.0),
            "drift_Etot_mRy_per_ps": float(np.polyfit(t_ps[mask], etot[mask], 1)[0] * 1000.0),
        }

    volume = BOX_AIMD ** 3
    mass_g = 32.0 * 18.01528 / 6.02214076e23
    density = mass_g / (volume * 1.0e-24)
    kB, xi, eta = 1.380649e-23, 2.837297, 8.90e-4
    yeh = {}
    for name, length in (("AIMD_32", BOX_AIMD * 1e-10), ("classical_512", 24.8 * 1e-10)):
        yeh[name] = kB * 300.0 * xi / (6.0 * np.pi * eta * length)

    return {
        "rounds": rounds,
        "n_steps": int(len(steps)),
        "duration_ps": float(t_ps[-1]),
        "temperature_K": {"mean": float(temp.mean()), "std": float(temp.std(ddof=1))},
        "first_0.5ps": segment(t_ps <= 0.5),
        "after_0.5ps": segment(t_ps > 0.5),
        "last_2ps": segment(t_ps >= t_ps[-1] - 2.0),
        "density_g_cm3": float(density),
        "box_A": BOX_AIMD,
        "msd_lag_ps": lag_ps,
        "msd": msd,
        "single_ps": single_ps,
        "single_msd": single,
        "windows": windows,
        "window_D_min": float(ds.min()),
        "window_D_max": float(ds.max()),
        "yeh_hummer_experimental_eta_m2_s": yeh,
        "series": {
            "t_ps": t_ps,
            "temp": temp,
            "etot": etot,
            "ekin": ekin,
            "conserved": conserved,
            "nose": nose,
        },
    }


def mean_std(values):
    arr = np.asarray(values, dtype=float)
    return float(arr.mean()), float(arr.std(ddof=1))


def plot_msd(classical_curves, deep_curves):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    colors = ["#1f4e79", "#2e75b6", "#9dc3e6"]
    for ax, curves, title in (
        (axes[0], classical_curves, "TIP4P/2005, 2 fs"),
        (axes[1], deep_curves, "DeepMD, 0.5 fs"),
    ):
        for (label, t, msd, fit), color in zip(curves, colors):
            ax.plot(t, msd, color=color, lw=1.3, label=label)
            mask = (t >= 50) & (t <= 200)
            line = fit["slope_A2_per_ps"] * t[mask] + fit["intercept_A2"]
            ax.plot(t[mask], line, color=color, lw=1.0, ls="--")
        ax.axvspan(50, 200, color="#f4e4c4", alpha=0.35, zorder=0)
        ax.set_xlabel("Production time (ps)")
        ax.set_ylabel(r"Oxygen MSD ($\mathrm{\AA}^2$)")
        ax.set_title(title)
        ax.legend(frameon=False)
        ax.grid(alpha=0.25)
    fig.savefig(PAPER / "msd_classical_deepmd.png", dpi=180)
    plt.close(fig)


def plot_diffusion(classical, deep, aimd):
    labels = ["Classical MD", "DeepMD", "AIMD"]
    means = []
    err_low = []
    err_high = []
    for rows in (classical, deep):
        ds = [row["fit_50_200ps"]["D_m2_s"] for row in rows]
        mean, std = mean_std(ds)
        means.append(mean * 1e9)
        err_low.append(std * 1e9)
        err_high.append(std * 1e9)
    t = aimd["single_ps"]
    msd = aimd["single_msd"]
    whole = fit_window(t, msd, 0.5, float(t[-1]))
    pieces = [
        fit_window(t, msd, lo, hi)["D_m2_s"] * 1e9
        for lo, hi in ((0.5, 2.5), (2.5, 4.5), (4.5, float(t[-1])))
    ]
    aimd_d = whole["D_m2_s"] * 1e9
    means.append(aimd_d)
    err_low.append(aimd_d - min(pieces))
    err_high.append(max(pieces) - aimd_d)
    fig, ax = plt.subplots(figsize=(6.6, 4.2), constrained_layout=True)
    x = np.arange(3)
    ax.bar(
        x, means, yerr=[err_low, err_high],
        color=["#1f4e79", "#c45911", "#548235"], capsize=5, width=0.62, zorder=2,
    )
    ax.axhline(EXP_D * 1e9, color="#833c0c", ls="--", lw=1.2, label=r"Experiment, 25 °C, $2.299\times10^{-9}$")
    ax.axhline(2.08, color="#833c0c", ls=":", lw=1.2, label=r"TIP4P/2005 literature, 298 K, $2.08\times10^{-9}$")
    ax.set_xticks(x, labels)
    ax.set_ylabel(r"$D$ ($10^{-9}\ \mathrm{m}^2/\mathrm{s}$)")
    ax.set_title("Self-diffusion coefficient")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(PAPER / "diffusion_comparison_aimd.png", dpi=180)
    plt.close(fig)


def centered_average(values, times, window_ps):
    """Centered average over a fixed time window. Ends use whatever samples fall inside the window."""
    values = np.asarray(values, dtype=float)
    times = np.asarray(times, dtype=float)
    half = 0.5 * window_ps
    cumulative = np.concatenate([[0.0], np.cumsum(values)])
    lower = np.searchsorted(times, times - half, side="left")
    upper = np.searchsorted(times, times + half, side="right")
    return (cumulative[upper] - cumulative[lower]) / np.maximum(upper - lower, 1)


def plot_lammps_thermo(series_list, path, title, energy_unit, xlabel, production_mark=None, average_ps=5.0):
    panels = (
        ("temp", "T (K)"),
        ("pe", rf"$E_\mathrm{{pot}}$ ({energy_unit})"),
        ("ke", rf"$E_\mathrm{{kin}}$ ({energy_unit})"),
        ("etotal", rf"$E_\mathrm{{tot}}$ ({energy_unit})"),
    )
    colors = ["#1f4e79", "#2e75b6", "#9dc3e6"]
    fig, axes = plt.subplots(4, 1, figsize=(8.2, 9.2), sharex=True, constrained_layout=True)
    for ax, (key, ylabel) in zip(axes, panels):
        for (label, series), color in zip(series_list, colors):
            t = series["t_ps"]
            y = series[key]
            ax.plot(t, y, color="#9a9a9a", lw=0.5, alpha=0.35, zorder=1)
            ax.plot(
                t, centered_average(y, t, average_ps),
                color=color, lw=1.6, label=label, zorder=3,
            )
        if key == "temp":
            ax.axhline(300, color="#833c0c", ls="--", lw=1, zorder=4)
        if production_mark is not None:
            ax.axvline(production_mark, color="0.45", ls=":", lw=1, zorder=4)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
    axes[0].legend(frameon=True, framealpha=0.92, edgecolor="none", ncol=3, loc="upper right", fontsize=8)
    axes[0].set_title(title)
    axes[-1].set_xlabel(xlabel)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_aimd(aimd, classical_curves, deep_curves):
    s = aimd["series"]
    t = s["t_ps"]
    panels = (
        (s["temp"], "T (K)", "#1f4e79", True),
        (s["etot"], r"DFT $E_\mathrm{tot}$ (Ry)", "#548235", False),
        (s["ekin"], r"Ionic $E_\mathrm{kin}$ (Ry)", "#2e75b6", False),
        (s["conserved"], r"$E_\mathrm{kin}+E_\mathrm{tot}$ (Ry)", "#c45911", False),
    )
    fig, axes = plt.subplots(4, 1, figsize=(8.2, 9.2), sharex=True, constrained_layout=True)
    for ax, (values, ylabel, color, show_target) in zip(axes, panels):
        ax.plot(t, values, color="#9a9a9a", lw=0.45, alpha=0.35, zorder=1)
        ax.plot(t, centered_average(values, t, 0.5), color=color, lw=1.6, zorder=3)
        if show_target:
            ax.axhline(300, color="#833c0c", ls="--", lw=1, zorder=4)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
    axes[0].set_title("AIMD, 32 waters, 6.86 ps")
    axes[-1].set_xlabel("Time (ps)")
    fig.savefig(PAPER / "aimd_thermo_avg.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    for curves, style, name in (
        (classical_curves, "-", "TIP4P/2005"),
        (deep_curves, "--", "DeepMD"),
    ):
        first = True
        for _label, t, msd, _fit in curves:
            mask = t <= 7.0
            ax.plot(
                t[mask], msd[mask], color="#7f7f7f", lw=1.15, ls=style, alpha=0.55,
                label=name if first else None, zorder=1,
            )
            first = False
    ax.plot(
        aimd["single_ps"], aimd["single_msd"], color="#1f4e79", lw=1.8,
        label="AIMD", zorder=3,
    )
    ax.set_xlim(0, 7.0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Production time (ps)")
    ax.set_ylabel(r"Oxygen MSD ($\mathrm{\AA}^2$)")
    ax.set_title("Oxygen MSD, first 7 ps")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.savefig(PAPER / "aimd_msd_early.png", dpi=180)
    plt.close(fig)


def slim(result):
    out = dict(result)
    out.pop("windows", None)
    return out


def main():
    print("Classical and DeepMD")
    classical, c_curves, c_thermo = analyze_classical()
    deep, d_curves, d_thermo = analyze_deepmd()
    print("AIMD")
    aimd = analyze_aimd()
    plot_msd(c_curves, d_curves)
    plot_diffusion(classical, deep, aimd)
    plot_lammps_thermo(
        c_thermo, PAPER / "classical_thermo_avg.png", "TIP4P/2005",
        "kcal/mol", "Production time (ps)",
    )
    plot_lammps_thermo(
        d_thermo, PAPER / "deepmd_thermo_avg.png", "DeepMD",
        "eV", "Time (ps)", production_mark=10.0,
    )
    plot_aimd(aimd, c_curves, d_curves)

    def pack(rows):
        ds = [row["fit_50_200ps"]["D_m2_s"] for row in rows]
        mean, std = mean_std(ds)
        return {"runs": [slim(row) | {"windows": row["windows"]} for row in rows],
                "D_mean_m2_s": mean, "D_sample_std_m2_s": std,
                "D_sem_m2_s": std / np.sqrt(len(ds))}

    payload = {
        "classical": pack(classical),
        "deepmd": pack(deep),
        "aimd": {
            "rounds": aimd["rounds"],
            "n_steps": aimd["n_steps"],
            "duration_ps": aimd["duration_ps"],
            "temperature_K": aimd["temperature_K"],
            "first_0.5ps": aimd["first_0.5ps"],
            "after_0.5ps": aimd["after_0.5ps"],
            "last_2ps": aimd["last_2ps"],
            "density_g_cm3": aimd["density_g_cm3"],
            "box_A": aimd["box_A"],
            "windows": [{k: v for k, v in row.items()} for row in aimd["windows"]],
            "window_D_min": aimd["window_D_min"],
            "window_D_max": aimd["window_D_max"],
            "yeh_hummer_experimental_eta_m2_s": aimd["yeh_hummer_experimental_eta_m2_s"],
        },
    }
    (PAPER / "project3_analysis.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
