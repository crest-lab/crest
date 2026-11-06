#!/bin/bash
# Single-point energy of 1-propanol with a method-independent implicit-solvation
# add-on: a GFN2-xTB gas-phase parent plus a composite solvation contribution
# (EEQ-BC charges + ddX/CPCM continuum + GFN2/ALPB nonpolar term) in water.
#
# The composite solvation calculator is configured through the TOML input only
# (there is no dedicated CLI flag for it).

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- TOML run ---
crest input.toml
