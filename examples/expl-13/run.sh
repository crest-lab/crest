#!/bin/bash
# Protonation site sampling of uracil with GFN2-xTB.
# Generates protomers by adding H+ to basic sites on the molecule.
# Expected output: 3 major protomers in protonated.xyz within 30 kcal/mol.

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
crest struc.xyz -protonate

# --- TOML run (equivalent settings) ---
# crest input.toml
