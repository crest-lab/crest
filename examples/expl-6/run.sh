#!/bin/bash
# Two-level iMTD-GC conformer search of 1-propanol using the A//B scheme:
# GFN-FF handles the fast MD/MTD sampling phase; GFN2 single-points
# re-rank the final ensemble.  This gives good accuracy at reduced cost.
# Output: crest_conformers.xyz

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# --- CLI run (A//B: B for sampling, A for single-point re-ranking) ---
crest struc.xyz -imtdgc --gfn2//gfnff -ewin 6.0

# Alternative: GFN-FF sampling + GFN2 geometry refinement of each conformer:
# crest struc.xyz -imtdgc --gfnff/opt/gfn2 -ewin 6.0

# --- TOML run (equivalent settings) ---
# crest input.toml
