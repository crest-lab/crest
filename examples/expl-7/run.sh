#!/bin/bash
# iMTD-GC conformer search of 1-propanol with GFN2-xTB and ALPB implicit solvation (water).
# Conformers in solution differ in relative energy from the gas phase.
# Uses 4 CPU threads.  Output: crest_conformers.xyz

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# --- CLI run ---
crest struc.xyz -gfn2 -alpb h2o -T 4 -ewin 2.0 -imtdgc

# --- TOML run (equivalent settings) ---
# crest input.toml
