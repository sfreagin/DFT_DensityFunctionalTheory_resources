#!/usr/bin/env python3
"""Oxygen MSD and self-diffusion from a Quantum ESPRESSO pw.x AIMD trajectory.

For Week 11's cubic 32-water cell (~9.8419 Angstrom per side):

  python analyze_qe_diffusion.py water32.out --box 9.8419 --dt-fs 0.5 \\
      --fit-start-ps 5 --fit-end-ps 20

Input may be (a) the FULL pw.x output containing repeated ATOMIC_POSITIONS
blocks, or (b) a standard multi-frame XYZ trajectory. Do NOT use the grep
thermodynamic output; it does not contain the coordinates.

Requires: numpy, matplotlib
Creates: qe_oxygen_msd.png and qe_oxygen_msd.csv (by default).

Algorithm:
  1. Read all O coordinates; QE output supports angstrom, bohr, alat,
     or crystal coordinates. The cell is assumed orthorhombic/constant.
  2. Unwrap crossing of periodic boundaries using minimum-image jumps.
  3. Calculate oxygen MSD averaged over ALL eligible time origins, via FFT.
  4. Fit MSD(t) = 6 D t + intercept in the chosen long-time window.

D = slope / 6. If slope is in Angstrom^2/ps, D in m^2/s is
(slope / 6) * 1e-8. Short trajectories yield very uncertain estimates.

Important: the separation between *saved frames* must be short enough that
no oxygen moves more than half a box length between consecutive snapshots.
Use an equilibrated trajectory, and select the fitting window by examining
whether the long-time MSD is approximately linear. Fixed-volume NVT only.
"""

import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BOHR_TO_ANG = 0.529177210903
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
STEP_RE = re.compile(r"Entering\s+Dynamics\s*:\s*iteration\s*=\s*(\d+)", re.I)
ALAT_RE = re.compile(r"lattice\s+parameter\s*\(alat\)\s*=\s*(" + NUMBER + r")\s*a\.u\.", re.I)
POS_RE = re.compile(r"^\s*ATOMIC_POSITIONS\s*(?:\(\s*([^\)]+)\s*\))?", re.I)


def qe_frames(path, nat, box, alat_override):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    frames = []
    steps = []
    step = None
    alat = alat_override
    i = 0
    while i < len(lines):
        line = lines[i]
        m = ALAT_RE.search(line)
        if m and alat_override is None:
            alat = float(m.group(1).replace("D", "E")) * BOHR_TO_ANG
        m = STEP_RE.search(line)
        if m:
            step = int(m.group(1))
        m = POS_RE.match(line)
        if m:
            unit = (m.group(1) or "alat").strip().lower().replace(" ", "")
            atoms = []
            for j in range(i + 1, min(i + 1 + nat, len(lines))):
                parts = lines[j].split()
                if len(parts) < 4:
                    break
                try:
                    xyz = [float(s.replace("D", "E").replace("d", "e")) for s in parts[1:4]]
                except ValueError:
                    break
                atoms.append((parts[0], xyz))
            if len(atoms) == nat:
                oxygen = np.asarray([xyz for element, xyz in atoms if re.fullmatch(r"O\d*", element, re.I)], dtype=float)
                if unit.startswith("angstrom"):
                    pass
                elif unit.startswith("bohr"):
                    oxygen *= BOHR_TO_ANG
                elif unit.startswith("alat"):
                    if alat is None:
                        raise ValueError("QE trajectory uses ATOMIC_POSITIONS (alat), but the "
                                         "lattice parameter was not found. Supply --alat-angstrom L.")
                    oxygen *= alat
                elif unit.startswith("crystal"):
                    oxygen *= box
                else:
                    raise ValueError(f"Unsupported QE coordinate unit: {unit!r}")
                frames.append(oxygen)
                steps.append(step)
                i += nat + 1
                continue
        i += 1
    if not frames:
        raise ValueError("Found no complete ATOMIC_POSITIONS blocks in pw.x output. "
                         "Ensure the trajectory is printed (e.g. iprint = 1).")
    return np.asarray(frames), steps


def xyz_frames(path, nat):
    frames = []
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].strip())
        except ValueError as exc:
            raise ValueError(f"Expected XYZ atom count at line {i + 1}") from exc
        if n != nat:
            raise ValueError(f"XYZ frame contains {n} atoms, but --nat is {nat}")
        if i + 2 + n > len(lines):
            raise ValueError("Incomplete final XYZ frame")
        atoms = []
        for line in lines[i + 2:i + 2 + n]:
            parts = line.split()
            if len(parts) < 4:
                raise ValueError("Invalid XYZ atomic-coordinate line")
            atoms.append((parts[0], [float(v) for v in parts[1:4]]))
        oxygen = np.asarray([xyz for element, xyz in atoms if re.fullmatch(r"O\d*", element, re.I)], dtype=float)
        frames.append(oxygen)
        i += n + 2
    if not frames:
        raise ValueError("No XYZ frames found")
    return np.asarray(frames), [None] * len(frames)


def unwrap_pbc(wrapped, box):
    wrapped = np.mod(wrapped, box)
    jumps = (wrapped[1:] - wrapped[:-1]) / box
    jumps -= np.rint(jumps)
    return np.concatenate([wrapped[:1], wrapped[:1] + np.cumsum(jumps * box, axis=0)], axis=0)


def multi_origin_msd(unwrapped, max_lag):
    """Unbiased all-time-origin MSD from FFT autocorrelations; result in A^2."""
    n, oxygen_count, _ = unwrapped.shape
    coords = unwrapped.reshape(n, -1)
    sq = np.einsum("ij,ij->i", coords, coords)
    prefix = np.r_[0.0, np.cumsum(sq)]
    lag = np.arange(max_lag + 1)
    sum_sq = prefix[n] - prefix[lag] + prefix[n - lag]
    acf = np.zeros(max_lag + 1, dtype=float)
    fft_length = 2 * n  # zero-padding prevents circular-correlation aliasing
    for col in range(coords.shape[1]):
        ft = np.fft.rfft(coords[:, col], n=fft_length)
        acf += np.fft.irfft(ft * np.conj(ft), n=fft_length)[:max_lag + 1]
    result = (sum_sq - 2.0 * acf) / ((n - lag) * oxygen_count)
    result[0] = 0.0
    return np.maximum(result, 0.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("trajectory", help="Full QE pw.x output or multi-frame XYZ file")
    parser.add_argument("--box", type=float, nargs="+", required=True, metavar="L",
                        help="Periodic box in Angstrom: one length for cubic, or Lx Ly Lz")
    parser.add_argument("--dt-fs", type=float, default=0.5,
                        help="QE ionic MD step in femtoseconds (default 0.5)")
    parser.add_argument("--frame-stride", type=int, default=1,
                        help="QE MD steps between frames if not encoded in trajectory (default 1)")
    parser.add_argument("--nat", type=int, default=96,
                        help="Number of atoms per snapshot (default 96 = 32 H2O)")
    parser.add_argument("--n-oxygen", type=int, default=32,
                        help="Expected oxygen atoms per frame (default 32)")
    parser.add_argument("--alat-angstrom", type=float, default=None,
                        help="Optional QE alat in Angstrom for ATOMIC_POSITIONS (alat)")
    parser.add_argument("--discard-ps", type=float, default=0.0,
                        help="Discard this many initial ps as equilibration (default 0)")
    parser.add_argument("--fit-start-ps", type=float, default=None,
                        help="Start of linear fit interval in ps (default 20%% of max lag)")
    parser.add_argument("--fit-end-ps", type=float, default=None,
                        help="End of linear fit interval in ps (default 80%% of max lag)")
    parser.add_argument("--max-lag-fraction", type=float, default=0.5,
                        help="Maximum MSD lag as fraction of saved trajectory (default 0.5)")
    parser.add_argument("--prefix", default="qe_oxygen_msd", help="Output file prefix")
    args = parser.parse_args()

    if len(args.box) not in (1, 3) or any(x <= 0 for x in args.box):
        parser.error("--box requires one positive length or three positive lengths")
    box = np.array(args.box * 3 if len(args.box) == 1 else args.box, dtype=float)
    if args.dt_fs <= 0 or args.frame_stride < 1 or args.nat < 1:
        parser.error("--dt-fs, --frame-stride, and --nat must be positive")
    if not (0 < args.max_lag_fraction <= 0.9):
        parser.error("--max-lag-fraction must be between 0 and 0.9")

    first_line = next((l.strip() for l in Path(args.trajectory).read_text(encoding="utf-8", errors="replace").splitlines() if l.strip()), "")
    is_xyz = first_line.isdigit()
    if is_xyz:
        frames, indices = xyz_frames(args.trajectory, args.nat)
    else:
        frames, indices = qe_frames(args.trajectory, args.nat, box, args.alat_angstrom)

    if frames.shape[1:] != (args.n_oxygen, 3):
        raise ValueError(f"Expected {args.n_oxygen} oxygen coordinates per frame, found {frames.shape[1]}")
    if len(frames) < 12:
        raise ValueError("Need at least 12 saved frames to calculate a meaningful MSD")

    if all(k is not None for k in indices) and len(indices) > 1:
        step_changes = np.diff(indices)
        if np.any(step_changes <= 0) or not np.all(step_changes == step_changes[0]):
            raise ValueError("QE step indices are not strictly increasing at a uniform interval. "
                             "Analyze each continuous trajectory separately or export a regular XYZ.")
        frames_per_step = int(step_changes[0])
    else:
        frames_per_step = args.frame_stride
        print(f"No consistent QE iteration data: assuming {frames_per_step} MD steps per saved frame.")
    frame_ps = args.dt_fs * frames_per_step / 1000.0

    discard_n = int(np.ceil(args.discard_ps / frame_ps))
    frames = frames[discard_n:]
    if len(frames) < 12:
        raise ValueError("Too few frames remain after discarding equilibration")
    n = len(frames)
    full_ps = (n - 1) * frame_ps
    if full_ps < 10:
        print("WARNING: production trajectory is shorter than 10 ps; diffusion estimate is not reliable.")

    unwrapped = unwrap_pbc(frames, box)
    max_lag = min(int(n * args.max_lag_fraction), n - 1)
    msd = multi_origin_msd(unwrapped, max_lag)
    lag = np.arange(max_lag + 1)
    time_ps = lag * frame_ps

    fit_start = args.fit_start_ps if args.fit_start_ps is not None else 0.2 * time_ps[-1]
    fit_end = args.fit_end_ps if args.fit_end_ps is not None else 0.8 * time_ps[-1]
    fit_mask = (time_ps >= fit_start) & (time_ps <= fit_end)
    if not (0 <= fit_start < fit_end <= time_ps[-1]) or np.count_nonzero(fit_mask) < 5:
        raise ValueError(f"Invalid fit window: choose 0 <= start < end <= {time_ps[-1]:.3f} ps, with >= 5 points")

    slope, intercept = np.polyfit(time_ps[fit_mask], msd[fit_mask], 1)
    fitted = slope * time_ps + intercept
    D_m2_s = slope * 1e-8 / 6.0
    if D_m2_s <= 0:
        print("WARNING: nonpositive diffusion slope; select a different fit window or acquire more data.")

    with open(f"{args.prefix}.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["lag_frames", "lag_time_ps", "oxygen_MSD_angstrom2"])
        writer.writerows(zip(lag, time_ps, msd))

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.plot(time_ps, msd, label="Oxygen MSD (all time origins)")
    ax.plot(time_ps[fit_mask], fitted[fit_mask], linestyle="--", linewidth=2,
            label=f"Linear fit: D = {D_m2_s:.3e} m$^2$/s")
    ax.set(xlabel="Lag time (ps)", ylabel="Oxygen MSD ($\\mathrm{\\AA}^2$)",
           title="Quantum ESPRESSO AIMD: Oxygen self-diffusion")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(f"{args.prefix}.png", dpi=200)
    plt.close(fig)

    print(f"Frames: {n}; oxygen atoms: {args.n_oxygen}; saved-frame interval: {frame_ps:.5f} ps")
    print(f"Production duration after discard: {full_ps:.3f} ps")
    print(f"Fit interval: {fit_start:.3f} to {fit_end:.3f} ps")
    print(f"MSD slope: {slope:.6g} Angstrom^2/ps")
    print(f"Self-diffusion coefficient D: {D_m2_s:.5e} m^2/s")
    print(f"Wrote {args.prefix}.png and {args.prefix}.csv")
    print("Check the MSD fit window and long-trajectory convergence before reporting D.")


if __name__ == "__main__":
    main()
