#!/bin/bash
# Non-covalent interaction (NCI / iMTD-NCI) conformer sampling of a water trimer.
# A wall potential is generated automatically to prevent cluster dissociation.
# Output: crest_conformers.xyz, crest_rotamers.xyz

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# --- CLI run ---
crest struc.xyz --gfnff --imtdgc -nci

# --- TOML run (equivalent settings) ---
# crest input.toml
