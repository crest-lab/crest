#!/bin/bash
# Geometry optimization followed by numerical Hessian (vibrational frequencies)
# of 1-propanol with GFN2-xTB.
# Output: crest_best.xyz (optimized structure), vibspectrum (frequencies)

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
crest struc.xyz -ohess

# --- TOML run (equivalent settings) ---
# crest input.toml
