#!/bin/bash
# Geometry optimisation of the caffeine molecule using the MACE-OFF23 organic force
# field, served through the fmlip-relay persistent backend.
# fmlip-relay spawns a persistent server process; CREST communicates with it
# over a local TCP socket, avoiding repeated Python startup overhead.
# Output: crest_best.xyz (optimised geometry)

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# Check for the fmlip-relay server; suggest pip install if absent
if ! command -v fmlip-relay-server >/dev/null 2>&1; then
  echo >&2 ""
  echo >&2 "ERROR: 'fmlip-relay-server' not found."
  echo >&2 "Install fmlip-relay (with MACE extras) from the CREST subproject directory:"
  echo >&2 ""
  echo >&2 "    pip install \"../../subprojects/fmlip_relay[mace]\""
  echo >&2 ""
  echo >&2 "MACE-OFF23 checkpoints download automatically on first use (cache: ~/.cache/mace)."
  exit 1
fi

# --- TOML run ---
# (No CLI equivalent; mlip settings are TOML-only.)
crest input.toml
