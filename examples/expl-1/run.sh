#!/bin/bash
# Single-point GFN2-xTB energy evaluation of 1-propanol.
# Output: energy printed to stdout; no structure files are written.

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
crest struc.xyz -sp

# --- TOML run (equivalent settings) ---
# crest input.toml
