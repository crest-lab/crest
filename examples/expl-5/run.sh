#!/bin/bash
# Default iMTD-GC conformer search of 1-propanol with GFN-FF.
# Expected output: ~4 unique conformers in crest_conformers.xyz within 2.0 kcal/mol.

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
crest struc.xyz -imtdgc -gfnff -ewin 2.0

# --- TOML run (equivalent settings) ---
# crest input.toml
