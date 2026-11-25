#!/usr/bin/env python3
"""Fourier content of the AIMD thermodynamic traces.

The figure is aimd_thermo_fft.png. It is not part of the report.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_project3 as ap

# Nose–Hoover frequency set in 32h2o_qe.in (fnosep = 10 THz, period 100 fs).
# Sines slower than this are what the eye reads as the periodic pattern.
F_CUT_THZ = 10.0
ZOOM = (2.0, 3.2)
PAPER = Path(__file__).resolve().parent


def load_series():
    thermo = []
    for name in ap.AIMD_ROUNDS:
        thermo.extend(ap.parse_thermo(ap.AIMD_ROOT / name / "qe_thermo.txt"))
    thermo = sorted(set(thermo))
    steps = np.array([row[0] for row in thermo], dtype=int)
    if np.any(np.diff(steps) != 1):
        raise RuntimeError("AIMD thermo iterations are not a continuous unit sequence")
    t_ps = (steps - 1) * ap.DT_AIMD_FS / 1000.0
    columns = {
        "temp": np.array([row[1] for row in thermo]),
        "etot": np.array([row[2] for row in thermo]),
        "ekin": np.array([row[3] for row in thermo]),
        "conserved": np.array([row[4] for row in thermo]),
    }
    return t_ps, columns


def peak_indices(freq, amp, count=2, min_df=0.45):
    """Tallest local maxima, kept far enough apart to be different oscillations."""
    candidates = []
    for i in range(2, len(amp) - 1):
        if freq[i] < 0.4:
            continue
        if amp[i] >= amp[i - 1] and amp[i] > amp[i + 1]:
            candidates.append(i)
    candidates.sort(key=lambda i: amp[i], reverse=True)
    chosen = []
    for i in candidates:
        if all(abs(freq[i] - freq[j]) >= min_df for j in chosen):
            chosen.append(i)
        if len(chosen) == count:
            break
    if len(chosen) < count:
        raise RuntimeError("Fewer distinct Fourier peaks than requested")
    return chosen


def sine_on_trend(spec, trend, index):
    one = np.zeros_like(spec)
    one[index] = spec[index]
    return trend + np.fft.irfft(one, n=len(trend))


def fourier(t_ps, values):
    """Linear-detrended FFT, the slow sum, and the two strongest sines."""
    n = len(values)
    dt = float(np.median(np.diff(t_ps)))
    trend = np.polyval(np.polyfit(t_ps, values, 1), t_ps)
    resid = values - trend
    spec = np.fft.rfft(resid)
    freq = np.fft.rfftfreq(n, d=dt)

    keep = freq <= F_CUT_THZ
    band = np.zeros_like(spec)
    band[keep] = spec[keep]
    slow = trend + np.fft.irfft(band, n=n)

    window = np.hanning(n)
    shown = np.fft.rfft(resid * window)
    amp = np.abs(shown) / window.sum()
    amp[1:-1] *= 2.0
    indices = peak_indices(freq, amp)
    sines = [sine_on_trend(spec, trend, int(k)) for k in indices]
    frac = float(np.var(slow - trend) / np.var(resid))
    return {
        "freq": freq,
        "amp": amp,
        "slow": slow,
        "sines": sines,
        "peaks": [(float(freq[k]), float(amp[k])) for k in indices],
        "frac": frac,
    }


def mark_peaks(ax, peaks, color):
    """Mark the two sines. A filled dot is the stronger one."""
    for (freq, amp), filled in zip(peaks, (True, False)):
        ax.plot(
            freq, amp, "o", ms=4.2, zorder=4,
            color=color, mfc=color if filled else "white", mec=color, mew=1.15,
        )
    (f1, a1), (f2, a2) = peaks
    if abs(f1 - f2) < 1.8:
        fright, aright = max(peaks, key=lambda item: item[0])
        text = "\n".join(f"{freq:.2f} THz, {1000.0 / freq:.0f} fs" for freq, _amp in peaks)
        ax.annotate(
            text, (fright, aright),
            textcoords="offset points", xytext=(8, 4),
            ha="left", va="bottom", fontsize=7.5, color=color,
        )
        return
    for freq, amp in peaks:
        ax.annotate(
            f"{freq:.2f} THz\n{1000.0 / freq:.0f} fs",
            (freq, amp),
            textcoords="offset points", xytext=(6, 4),
            ha="left", va="bottom", fontsize=7.5, color=color,
        )


def main():
    t_ps, columns = load_series()
    panels = (
        ("temp", "T (K)", "#1f4e79"),
        ("etot", r"DFT $E_\mathrm{tot}$ (Ry)", "#548235"),
        ("ekin", r"Ionic $E_\mathrm{kin}$ (Ry)", "#2e75b6"),
        ("conserved", r"$E_\mathrm{kin}+E_\mathrm{tot}$ (Ry)", "#c45911"),
    )
    packed = [(ylabel, color, fourier(t_ps, columns[key])) for key, ylabel, color in panels]

    fig, axes = plt.subplots(
        4, 3, figsize=(12.6, 9.6), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.35, 1.05, 1.15]},
    )
    zoom = (t_ps >= ZOOM[0]) & (t_ps <= ZOOM[1])
    for row, (ylabel, color, spec) in enumerate(packed):
        values = columns[panels[row][0]]
        ax_full, ax_zoom, ax_spec = axes[row]

        for ax, mask in ((ax_full, slice(None)), (ax_zoom, zoom)):
            ax.plot(t_ps[mask], values[mask], color="#9a9a9a", lw=0.45, alpha=0.85, zorder=1)
            ax.plot(t_ps[mask], spec["slow"][mask], color=color, lw=1.15, zorder=3)
            ax.grid(alpha=0.25)
            ax.set_ylabel(ylabel)
        for sine, style in zip(spec["sines"], ("--", "-.")):
            ax_zoom.plot(t_ps[zoom], sine[zoom], color=color, lw=1.05, ls=style, alpha=0.95, zorder=2)
        ax_full.text(
            0.99, 0.96,
            f"{spec['frac']:.0%} of the fluctuation\nis slower than 10 THz",
            transform=ax_full.transAxes, ha="right", va="top", fontsize=7.5, color=color,
        )

        freq, amp = spec["freq"], spec["amp"]
        shown = freq <= 20.0
        ax_spec.fill_between(freq[shown], amp[shown], color=color, alpha=0.18, zorder=1)
        ax_spec.plot(freq[shown], amp[shown], color=color, lw=0.9, zorder=2)
        ax_spec.axvspan(0.0, F_CUT_THZ, color=color, alpha=0.06, zorder=0)
        ax_spec.axvline(F_CUT_THZ, color="#833c0c", ls="--", lw=1.0, zorder=3)
        mark_peaks(ax_spec, spec["peaks"], color)
        ax_spec.set_xlim(0.0, 20.0)
        ax_spec.set_ylim(0.0, float(amp[shown].max()) * 1.32)
        ax_spec.set_ylabel("Amplitude")
        ax_spec.grid(alpha=0.25)

        if row < 3:
            for ax in (ax_full, ax_zoom, ax_spec):
                ax.tick_params(labelbottom=False)
        peak_text = ", ".join(
            f"{freq:.3f} THz ({1000.0 / freq:.1f} fs)" for freq, _amp in spec["peaks"]
        )
        print(f"{panels[row][0]}: {peak_text}; fraction below 10 THz = {spec['frac']:.3f}")

    axes[0, 0].plot([], [], color="#9a9a9a", lw=1.2, label="recorded")
    axes[0, 0].plot([], [], color="#1f4e79", lw=1.4, label="sum of sines below 10 THz")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="upper left")
    axes[0, 1].plot([], [], color="#1f4e79", lw=1.05, ls="--", label="strongest sine")
    axes[0, 1].plot([], [], color="#1f4e79", lw=1.05, ls="-.", label="second sine")
    axes[0, 1].legend(frameon=False, fontsize=8, loc="upper right")
    axes[0, 2].plot([], [], color="#833c0c", ls="--", lw=1.0, label="Nose–Hoover, 10 THz")
    axes[0, 2].legend(frameon=False, fontsize=8, loc="upper right")

    axes[0, 0].set_title("Full trajectory")
    axes[0, 1].set_title("2.0–3.2 ps")
    axes[0, 2].set_title("Fourier amplitude")
    axes[-1, 0].set_xlabel("Time (ps)")
    axes[-1, 1].set_xlabel("Time (ps)")
    axes[-1, 2].set_xlabel("Frequency (THz)")
    axes[0, 0].set_xlim(t_ps[0], 7.0)
    axes[0, 1].set_xlim(*ZOOM)

    period_axis = axes[0, 2].secondary_xaxis(
        "top",
        functions=(lambda f: 1000.0 / np.maximum(f, 1e-6), lambda p: 1000.0 / np.maximum(p, 1e-6)),
    )
    period_axis.set_xticks([1000, 200, 100, 50])
    period_axis.set_xlabel("Period (fs)")

    fig.suptitle("AIMD, 32 waters: each thermodynamic trace as a sum of sine waves", fontsize=13)
    fig.savefig(PAPER / "aimd_thermo_fft.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
