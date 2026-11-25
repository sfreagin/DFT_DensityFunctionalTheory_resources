#!/usr/bin/env python3
"""Review figures for section 7. Writes png files in this directory only."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import analyze_project3 as ap

PAPER = ap.PAPER
BOX_512 = 24.8
CLASSICAL_TRJ = ap.CLASSICAL[0][1] / "water512_NVT.lammpstrj"
DEEPMD_TRJ = ap.DEEPMD[0][1] / "deepmd_water.lammpstrj"


def read_lammpstrj(path, oxygen_type=1, has_mol=False):
    frames = []
    steps = []
    with open(path) as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(handle.readline())
            handle.readline()  # ITEM: NUMBER OF ATOMS
            n_atoms = int(handle.readline())
            handle.readline()  # ITEM: BOX BOUNDS
            for _ in range(3):
                handle.readline()
            header = handle.readline()
            oxygen = []
            for _ in range(n_atoms):
                parts = handle.readline().split()
                atom_type = int(parts[2] if has_mol else parts[1])
                if atom_type == oxygen_type:
                    xyz = parts[3:6] if has_mol else parts[2:5]
                    oxygen.append([float(v) for v in xyz])
            if len(oxygen) == 0:
                raise RuntimeError(f"No oxygen atoms in frame {step} of {path}")
            steps.append(step)
            frames.append(oxygen)
    return np.asarray(steps, dtype=int), np.asarray(frames, dtype=float)


def goo(frames, box, r_max, dr=0.05):
    edges = np.arange(dr, r_max + dr * 0.5, dr)
    centers = edges[:-1] + 0.5 * dr
    hist = np.zeros(len(centers))
    n_frames, n_ox, _ = frames.shape
    for frame in frames:
        delta = frame[:, None, :] - frame[None, :, :]
        delta -= box * np.round(delta / box)
        dist = np.sqrt(np.sum(delta * delta, axis=-1))
        iu = np.triu_indices(n_ox, k=1)
        hist += np.histogram(dist[iu], bins=edges)[0]
    rho = n_ox / box**3
    shell = 4.0 * np.pi * centers**2 * dr
    g = (2.0 * hist) / (n_frames * n_ox * rho * shell)
    return centers, g


def running_slope(t, y, half):
    slopes = np.full(len(t), np.nan)
    for i, ti in enumerate(t):
        mask = (t >= ti - half) & (t <= ti + half)
        if np.count_nonzero(mask) < 5:
            continue
        if t[mask][-1] - t[mask][0] < half:
            continue
        slopes[i] = np.polyfit(t[mask], y[mask], 1)[0]
    return slopes


def oxygen_vacf(unwrapped, dt_ps, max_lag):
    """Central-difference oxygen velocities, COM already removed. Returns lag and C(t)."""
    vel = (unwrapped[2:] - unwrapped[:-2]) / (2.0 * dt_ps)
    n = len(vel)
    max_lag = min(max_lag, n // 2)
    # Average v(t)·v(t+lag) over atoms and origins.
    dots = np.einsum("taf,taf->ta", vel, vel)
    c0 = float(dots.mean())
    corr = np.empty(max_lag + 1)
    corr[0] = c0
    for lag in range(1, max_lag + 1):
        corr[lag] = np.mean(np.einsum("taf,taf->ta", vel[:-lag], vel[lag:]))
    lag_ps = np.arange(max_lag + 1) * dt_ps
    # D(t) = (1/3) int_0^t <v·v> dt, in Å^2/ps, then to m^2/s.
    integral = np.cumsum(corr) * dt_ps
    integral[1:] -= 0.5 * (corr[1:] + corr[0]) * dt_ps  # trapezoid correction vs left Riemann
    # Recalculate trapezoid properly.
    trap = np.zeros_like(corr)
    trap[1:] = np.cumsum(0.5 * (corr[:-1] + corr[1:]) * dt_ps)
    d_m2_s = trap / 3.0 * 1.0e-8
    return lag_ps, corr / c0, d_m2_s


def plot_goo(classical, deep, aimd):
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax.plot(*goo(classical, BOX_512, 8.0), color="#1f4e79", lw=1.6, label="Classical MD, seed 1")
    ax.plot(*goo(deep, BOX_512, 8.0), color="#c45911", lw=1.6, label="DeepMD, seed 1")
    r_aimd, g_aimd = goo(aimd, ap.BOX_AIMD, ap.BOX_AIMD / 2.0)
    ax.plot(r_aimd, g_aimd, color="#548235", lw=1.6, label="AIMD")
    ax.axhline(1.0, color="0.6", lw=0.8)
    ax.set_xlim(0, 8)
    ax.set_xlabel(r"$r$ (Å)")
    ax.set_ylabel(r"$g_{\mathrm{OO}}(r)$")
    ax.set_title("Oxygen–oxygen radial distribution")
    ax.legend(frameon=False)
    ax.grid(alpha=0.25)
    fig.savefig(PAPER / "goo_comparison.png", dpi=180)
    plt.close(fig)


def plot_running_slope(classical_curves, deep_curves, aimd):
    fig, ax = plt.subplots(figsize=(7.4, 4.4), constrained_layout=True)
    colors = ["#1f4e79", "#2e75b6", "#9dc3e6"]
    for (label, t, msd, _fit), color in zip(classical_curves, colors):
        slope = running_slope(t, msd, 2.0)
        mask = (t <= 15) & np.isfinite(slope)
        ax.plot(t[mask], slope[mask], color=color, lw=1.3, label=f"Classical {label}")
    deep_colors = ["#c45911", "#e07a3d", "#f0b27a"]
    for (label, t, msd, _fit), color in zip(deep_curves, deep_colors):
        slope = running_slope(t, msd, 2.0)
        mask = (t <= 15) & np.isfinite(slope)
        ax.plot(t[mask], slope[mask], color=color, lw=1.3, label=f"DeepMD {label}")
    t = aimd["single_ps"]
    slope = running_slope(t, aimd["single_msd"], 0.5)
    mask = np.isfinite(slope)
    ax.plot(t[mask], slope[mask], color="#548235", lw=1.5, label="AIMD")
    ax.axhline(2.299e-9 / 1e-8 * 6.0, color="#833c0c", ls="--", lw=1, label="Experiment, 6D")
    ax.set_xlim(0, 15)
    ax.set_ylim(-0.5, 4.0)
    ax.set_xlabel("Production time (ps)")
    ax.set_ylabel(r"Running slope of MSD (Å$^2$/ps)")
    ax.set_title("Running slope of the oxygen MSD")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.grid(alpha=0.25)
    fig.savefig(PAPER / "running_slope.png", dpi=180)
    plt.close(fig)


def plot_green_kubo(aimd_pos, deep_pos, deep_dt):
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True, constrained_layout=True)
    series = [("AIMD", aimd_pos, 0.005, "#548235")]
    if deep_dt is not None and deep_dt <= 0.05:
        series.append(("DeepMD, seed 1", deep_pos, deep_dt, "#c45911"))
    for label, pos, dt, color in series:
        lag, vacf, d_run = oxygen_vacf(pos, dt, max_lag=int(round(2.0 / dt)))
        axes[0].plot(lag, vacf, color=color, lw=1.4, label=label)
        axes[1].plot(lag, d_run * 1e9, color=color, lw=1.4, label=label)
    axes[0].axhline(0, color="0.6", lw=0.8)
    axes[0].set_ylabel(r"Normalized $\langle v(0)\cdot v(t)\rangle$")
    axes[0].set_title("Oxygen velocity autocorrelation")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)
    axes[1].axhline(2.299, color="#833c0c", ls="--", lw=1, label="Experiment")
    axes[1].set_xlim(0, 2.0)
    axes[1].set_xlabel("Lag (ps)")
    axes[1].set_ylabel(r"Running $D$ ($10^{-9}$ m$^2$/s)")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.25)
    fig.savefig(PAPER / "green_kubo.png", dpi=180)
    plt.close(fig)
    return deep_dt is not None and deep_dt <= 0.05


def plot_blocking():
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2), sharey=False, constrained_layout=True)
    windows = [(50, 80), (80, 110), (110, 140), (140, 170), (170, 200)]
    mids = [0.5 * (lo + hi) for lo, hi in windows]
    colors = ["#1f4e79", "#2e75b6", "#9dc3e6"]
    for ax, specs, scale, title in (
        (axes[0], ap.CLASSICAL, lambda steps: steps * 2.0 / 1000.0, "Classical MD"),
        (axes[1], ap.DEEPMD, lambda steps: (steps - steps[0]) * 0.0005, "DeepMD"),
    ):
        seed_means = []
        for (label, folder), color in zip(specs, colors):
            filename = "msd_NVT.dat" if "week8" in folder.name else "msd_O.dat"
            steps, msd = ap.load_msd(folder / filename)
            t = scale(steps)
            ds = [ap.fit_window(t, msd, lo, hi)["D_m2_s"] * 1e9 for lo, hi in windows]
            seed_means.append(float(np.mean(ds)))
            ax.plot(mids, ds, marker="o", color=color, lw=1.3, label=label)
        mean, std = ap.mean_std([
            ap.fit_window(
                scale(ap.load_msd(folder / ("msd_NVT.dat" if "week8" in folder.name else "msd_O.dat"))[0]),
                ap.load_msd(folder / ("msd_NVT.dat" if "week8" in folder.name else "msd_O.dat"))[1],
                50, 200,
            )["D_m2_s"] * 1e9
            for _label, folder in specs
        ])
        ax.axhspan(mean - std, mean + std, color="#1f4e79", alpha=0.12, label=f"Seed spread, ±{std:.2f}")
        ax.axhline(mean, color="#1f4e79", lw=1.0)
        ax.set_title(title)
        ax.set_xlabel("Center of the 30 ps window (ps)")
        ax.set_ylabel(r"$D$ ($10^{-9}$ m$^2$/s)")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.25)
    fig.savefig(PAPER / "diffusion_blocking.png", dpi=180)
    plt.close(fig)


def main():
    print("Reading trajectories")
    _c_steps, classical = read_lammpstrj(CLASSICAL_TRJ, has_mol=True)
    d_steps, deep = read_lammpstrj(DEEPMD_TRJ, has_mol=False)
    deep_dt = float(np.median(np.diff(d_steps))) * 0.0005
    print(f"  classical frames {len(classical)}, deep frames {len(deep)}, deep dt {deep_dt:.4f} ps")

    print("AIMD positions")
    frames = {}
    for name in ap.AIMD_ROUNDS:
        frames.update(ap.parse_oxygen(ap.AIMD_ROOT / name / "32h2o_qe.out"))
    order = np.array(sorted(frames))
    aimd = np.asarray([frames[int(s)] for s in order], dtype=float)
    aimd_rdf = aimd[::20]
    aimd_unwrapped = ap.unwrap(aimd[::10], ap.BOX_AIMD)

    print("gOO")
    plot_goo(classical[::2], deep[::2], aimd_rdf)
    print("Running slope")
    _c_res, c_curves, _c_th = ap.analyze_classical()
    _d_res, d_curves, _d_th = ap.analyze_deepmd()
    aimd_analysis = ap.analyze_aimd()
    plot_running_slope(c_curves, d_curves, aimd_analysis)
    print("Green–Kubo")
    deep_used = plot_green_kubo(aimd_unwrapped, ap.unwrap(deep, BOX_512), deep_dt)
    print("  deepmd included" if deep_used else "  deepmd frames are too far apart for a velocity autocorrelation")
    print("Blocking")
    plot_blocking()
    print("done")


if __name__ == "__main__":
    main()
