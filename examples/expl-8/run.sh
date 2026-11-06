#!/bin/bash
# Quick iMTD-GC conformer search of 1-propanol with reduced MTD simulation length.
# Useful for a fast first survey of conformational space.
# Output: crest_conformers.xyz
# Also available: -squick (super quick) and -mquick (mega quick) modes.
#
# -finalhess adds a post-search Hessian on each conformer in the final ensemble
# and re-ranks by Gibbs free energy.

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# --- CLI run (quick search only) ---
crest struc.xyz -quick -gfnff -ewin 2.0 -imtdgc -finalhess

# --- TOML run (equivalent settings) ---
# crest input.toml
