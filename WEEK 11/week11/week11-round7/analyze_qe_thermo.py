#!/usr/bin/env python3
"""Plot temperature and energies from Quantum ESPRESSO pw.x AIMD output.

Designed for the Week 11 (32-water, SCAN, Nose thermostat) hands-on.

First extract the relevant lines with the command in the tutorial:

  grep -E "Entering Dynamics|!    total energy|kinetic energy \\(Ekin\\)|temperature|Ekin \\+ Etot" water32.out > qe_thermo.txt

Then run:

  python analyze_qe_thermo.py qe_thermo.txt --dt-fs 0.5

You can also pass the FULL pw.x output directly in place of qe_thermo.txt.

Requires: numpy, matplotlib
Creates: qe_thermo.png and qe_thermo.csv (by default).

The QE output labels the Born-Oppenheimer DFT potential energy as Etot,
ionic kinetic energy as Ekin, and their sum as 'Ekin + Etot (const)'.
Under Nose NVT dynamics, the sum need not be constant because the thermostat
exchanges energy with the ions; the extended-system energy is a different
quantity. We reconstruct Etot = (Ekin + Etot) - Ekin using values from the
SAME ionic step, avoiding possible one-step offsets in SCF energy printouts.
All energy quantities below are in Rydberg (Ry).
"""

import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
STEP_RE = re.compile(r"Entering\s+Dynamics\s*:\s*iteration\s*=\s*(\d+)", re.I)
EKIN_RE = re.compile(r"kinetic\s+energy\s*\(Ekin\)\s*=\s*(" + NUMBER + r")", re.I)
TEMP_RE = re.compile(r"^\s*temperature\s*=\s*(" + NUMBER + r")\s*K\b", re.I)
SUM_RE = re.compile(r"Ekin\s*\+\s*Etot\s*\(const\)\s*=\s*(" + NUMBER + r")", re.I)


def parse_thermo(filename):
    records = []
    step = None
    pending = {}
    for line in Path(filename).read_text(encoding="utf-8", errors="replace").splitlines():
        match = STEP_RE.search(line)
        if match:
            step = int(match.group(1))
        match = EKIN_RE.search(line)
        if match:
            pending["ekin"] = float(match.group(1).replace("D", "E").replace("d", "e"))
        match = TEMP_RE.search(line)
        if match:
            pending["temp"] = float(match.group(1).replace("D", "E").replace("d", "e"))
        match = SUM_RE.search(line)
        if match:
            pending["sum"] = float(match.group(1).replace("D", "E").replace("d", "e"))
            if all(k in pending for k in ("ekin", "temp", "sum")):
                kinetic = pending["ekin"]
                total = pending["sum"]
                records.append((step, pending["temp"], total - kinetic, kinetic, total))
            pending.clear()

    if not records:
        raise ValueError(
            "No complete MD thermodynamic records found. Use the exact grep "
            "command in the script header, or supply the full pw.x output."
        )

    explicit = [row[0] is not None for row in records]
    if all(explicit):
        steps = np.asarray([row[0] for row in records], dtype=int)
        if np.any(np.diff(steps) <= 0):
            print("WARNING: QE iteration numbers repeat/decrease (perhaps joined restarts). "
                  "Using sequential step numbers for plotting.")
            steps = np.arange(1, len(records) + 1)
    else:
        # Some QE builds or grep choices omit the 'Entering Dynamics' lines.
        print("WARNING: not every frame has a QE iteration marker; using sequential steps.")
        steps = np.arange(1, len(records) + 1)
    thermo = np.asarray([row[1:] for row in records], dtype=float)
    return steps, thermo


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file", help="Output from the tutorial's grep command (or full pw.x output)")
    parser.add_argument("--dt-fs", type=float, default=0.5,
                        help="Ionic time step in femtoseconds (default: 0.5 fs, from dt=10.34)")
    parser.add_argument("--start-step", type=int, default=0,
                        help="Optional first step to plot (e.g. discard equilibration)")
    parser.add_argument("--prefix", default="qe_thermo", help="Output file prefix (default: qe_thermo)")
    args = parser.parse_args()
    if args.dt_fs <= 0:
        parser.error("--dt-fs must be positive")

    steps, data = parse_thermo(args.file)
    mask = steps >= args.start_step
    steps, data = steps[mask], data[mask]
    if len(steps) == 0:
        parser.error("No data remain after --start-step")

    temp, epot, ekin, esum = data.T
    time_ps = steps * args.dt_fs / 1000.0
    csv_name = f"{args.prefix}.csv"
    with open(csv_name, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["step", "time_ps", "temperature_K", "Etot_DFT_Ry", "Ekin_ionic_Ry", "Ekin_plus_Etot_Ry"])
        for values in zip(steps, time_ps, temp, epot, ekin, esum):
            writer.writerow([f"{v:.10g}" if isinstance(v, (float, np.floating)) else int(v) for v in values])

    fig, axes = plt.subplots(4, 1, figsize=(9, 10), sharex=True, constrained_layout=True)
    for ax, y, ylabel in zip(
        axes,
        (temp, epot, ekin, esum),
        ("Temperature (K)", "DFT Etot (Ry)", "Ionic Ekin (Ry)", "Ekin + Etot (Ry)")
    ):
        ax.plot(steps, y, linewidth=0.9)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
    axes[-1].set_xlabel("AIMD step")
    fig.suptitle("Quantum ESPRESSO AIMD: Thermodynamic quantities")
    fig.savefig(f"{args.prefix}.png", dpi=200)
    plt.close(fig)

    print(f"Parsed {len(steps)} complete ionic steps; step range {steps[0]}-{steps[-1]}.")
    print(f"Mean temperature: {np.mean(temp):.2f} K")
    print(f"Mean DFT Etot: {np.mean(epot):.6f} Ry")
    print(f"Mean ionic Ekin: {np.mean(ekin):.6f} Ry")
    print(f"Wrote {args.prefix}.png and {csv_name}")
    print("Note: Ekin + Etot is not expected to remain strictly constant for Nose NVT MD.")


if __name__ == "__main__":
    main()
