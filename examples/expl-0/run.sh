#!/bin/bash
# Dry run of CREST: prints settings and thresholds without running any calculation.
# Use this to preview the iMTD-GC setup before committing to a full run.
# Note: -dry is a CLI-only flag; see input.toml for the equivalent full-run settings.

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- CLI run ---
crest struc.xyz -dry

# --- TOML run (equivalent full run without -dry) ---
# crest input.toml
