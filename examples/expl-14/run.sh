#!/bin/bash
# Metal/ion adduct generation for alpha-D-glucose with GFN2-xTB.
# Replaces H+ with Cs+ (via -swel) to generate Cs+ adducts.
# Other ions (Na+, Li+, Ca2+, …) can be used the same way.
# Output: protonated.xyz

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
crest struc.xyz -protonate -swel Cs+

# --- TOML run (equivalent settings) ---
# crest input.toml
