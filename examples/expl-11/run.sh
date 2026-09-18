#!/bin/bash
# Optimize all structures in a trajectory/ensemble file with GFN2-xTB.
# Output: crest_ensemble.xyz (optimized structures, not sorted)

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
crest -mdopt xtb.trj

# --- TOML run (equivalent settings) ---
# crest input.toml
