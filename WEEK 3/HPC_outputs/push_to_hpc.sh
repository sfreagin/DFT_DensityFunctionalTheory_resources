#!/bin/bash
# Copy locally created Week 3 homework files to the HPC cluster.
# Run this from your Mac (not inside an already-open SSH session).
# You will be prompted for your SJSU HPC password.

set -euo pipefail

HOST="019119984@coe-hpc2.sjsu.edu"
LOCAL_ROOT="/Users/stephenreagin/Desktop/SJSU/MATE 250/WEEK3/HPC_outputs"
REMOTE_ROOT="mate250/week3/silicon/autobash"

cd "$LOCAL_ROOT"

scp \
  "ecut_convergence.png" \
  "convergence_caption.txt" \
  "${HOST}:${REMOTE_ROOT}/"

scp \
  "k_points/k_convergence.png" \
  "${HOST}:${REMOTE_ROOT}/k_points/"

echo "Done. On HPC, check:"
echo "  ls -la ~/${REMOTE_ROOT}/ecut_convergence.png ~/${REMOTE_ROOT}/convergence_caption.txt"
echo "  ls -la ~/${REMOTE_ROOT}/k_points/k_convergence.png"
