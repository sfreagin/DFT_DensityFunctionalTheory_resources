
"""
Make BORN.txt-like output from a Quantum ESPRESSO ph.x output.

Output format (like your BORN.txt):
  # epsilon and Z* of atoms 1 3 5 9 11
  <9 numbers for epsilon (flattened 3x3)>
  <9 numbers for Z* atom 1>
  <9 numbers for Z* atom 3>
  <9 numbers for Z* atom 5>
  <9 numbers for Z* atom 9>
  <9 numbers for Z* atom 11>

Usage:
  python qe_ph_to_born.py /path/to/bivo4.ph.out /path/to/BORN.txt

Optional: choose atoms
  python qe_ph_to_born.py in.out out.txt --atoms 1 3 5 9 11
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


EPS_TITLE = "Dielectric constant in cartesian axis"
ZSTAR_TITLE = "Effective charges (d Force / dE) in cartesian axis with asr applied:"

ATOM_RE = re.compile(r"^\s*atom\s+(\d+)\s+([A-Za-z]+)\b")

# Lines like:
#   E*x (        5.78212        0.08501       -0.00002 )
ZROW_RE = re.compile(
    r"""^\s*E\*[xyz]\s*\(\s*
        (?P<a>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?)\s+
        (?P<b>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?)\s+
        (?P<c>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?)\s*
        \)\s*$""",
    re.VERBOSE,
)

# Lines like:
#   (       7.323297815      -0.000351078       0.000003442 )
EPS_ROW_RE = re.compile(
    r"""^\s*\(\s*
        (?P<a>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?)\s+
        (?P<b>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?)\s+
        (?P<c>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?)\s*
        \)\s*$""",
    re.VERBOSE,
)


def find_line(lines: List[str], needle: str) -> int:
    for i, line in enumerate(lines):
        if needle in line:
            return i
    raise ValueError(f"Could not find section header: {needle!r}")


def parse_3x3_parentheses(lines: List[str], start: int, row_re: re.Pattern) -> Tuple[int, List[float]]:
    """Parse 3 consecutive '( a b c )' style rows starting at/after start."""
    rows: List[List[float]] = []
    i = start
    while i < len(lines) and len(rows) < 3:
        m = row_re.match(lines[i])
        if m:
            rows.append([float(m.group("a")), float(m.group("b")), float(m.group("c"))])
        i += 1
    if len(rows) != 3:
        raise ValueError("Failed to parse a complete 3x3 matrix (expected 3 rows).")
    flat = [v for r in rows for v in r]
    return i, flat


def extract_epsilon(lines: List[str]) -> List[float]:
    i = find_line(lines, EPS_TITLE) + 1
    _, eps_flat = parse_3x3_parentheses(lines, i, EPS_ROW_RE)
    return eps_flat


def extract_zstars(lines: List[str]) -> Dict[int, List[float]]:
    sec = find_line(lines, ZSTAR_TITLE) + 1
    i = sec
    z_by_atom: Dict[int, List[float]] = {}

    while i < len(lines):
        # skip blanks
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            break

        m_atom = ATOM_RE.match(lines[i])
        if not m_atom:
            break

        atom_idx = int(m_atom.group(1))
        i += 1

        # parse 3 E* rows
        i, flat = parse_3x3_parentheses(lines, i, ZROW_RE)

        # keep first occurrence for each atom index
        z_by_atom.setdefault(atom_idx, flat)

    if not z_by_atom:
        raise ValueError("Found Z* section header, but extracted zero tensors.")
    return z_by_atom


def fmt9(vals: List[float]) -> str:
    # Match BORN.txt style: fixed 8 decimals
    return "  ".join(f"{v: .8f}".replace(" -", "-").strip() for v in vals)


def parse_atoms(argv: List[str]) -> List[int]:
    if "--atoms" not in argv:
        return [1, 3, 5, 9, 11]
    j = argv.index("--atoms")
    atoms: List[int] = []
    for tok in argv[j + 1 :]:
        if tok.startswith("-"):
            break
        atoms.append(int(tok))
    if not atoms:
        raise ValueError("You provided --atoms but no atom indices.")
    return atoms


def main() -> int:
    if len(sys.argv) < 3:
        print("Usage: python qe_ph_to_born.py INPUT_PH_OUT OUTPUT_TXT [--atoms 1 3 5 9 11]", file=sys.stderr)
        return 2

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    atoms = parse_atoms(sys.argv[3:])

    lines = in_path.read_text(errors="replace").splitlines()

    eps = extract_epsilon(lines)
    z_by_atom = extract_zstars(lines)

    missing = [a for a in atoms if a not in z_by_atom]
    if missing:
        raise ValueError(f"Missing Z* tensors for atoms: {missing}. Available: {sorted(z_by_atom.keys())}")

    out_lines = []
    out_lines.append("# epsilon and Z* of atoms " + " ".join(str(a) for a in atoms))
    out_lines.append(fmt9(eps))
    for a in atoms:
        out_lines.append(fmt9(z_by_atom[a]))

    out_path.write_text("\n".join(out_lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
