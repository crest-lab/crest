#!/bin/bash
# Constrained iMTD-GC conformer search of 1-propanol:
# the C-C-C-O backbone (atoms 1-4) is frozen; only the OH dihedral is sampled.

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# --- CLI run ---
crest struc.xyz -gfnff -freeze 1-4 -ewin 2.0 --imtdgc

# --- TOML run (constraints defined directly in the input file) ---
# crest input.toml
