#!/bin/bash
# Standalone molecular dynamics (MD) simulation of 1-propanol with GFN-FF.
# Runs a 20 ps NVT trajectory at 400 K.
# Output: crest_dynamics.trj (trajectory), crest_property.out (energies)

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# --- CLI run ---
crest struc.xyz -dyn -gfnff -mdtemp 400 -mdlen 20 -tstep 1.0

# --- TOML run (equivalent settings) ---
# crest input.toml
