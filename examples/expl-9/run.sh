#!/bin/bash
# Standalone CREGEN sorting of an ensemble/trajectory file.
# Removes duplicate structures (by RMSD and energy) and sorts by energy.
# Output: crest_conformers.xyz (unique), crest_rotamers.xyz (all), xtb.trj.sorted

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# --- CLI run ---
crest -cregen xtb.trj -ewin 100.0

# --- TOML run (equivalent settings) ---
# crest input.toml
