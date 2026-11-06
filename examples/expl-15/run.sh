#!/bin/bash
# Tautomer screening of guanine via protonation/deprotonation sequences.
# Explores prototropic tautomers at GFN2-xTB level.
# Expected output: 5 major tautomers within 10 kcal/mol in tautomers.xyz.

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
crest struc.xyz -tautomerize -ewin 10.0

# --- TOML run (equivalent settings) ---
# crest input.toml
