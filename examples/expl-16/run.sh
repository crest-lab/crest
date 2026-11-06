#!/bin/bash
# Geometry optimisation of an Ar4 cluster using the Lennard-Jones potential
# served through the fmlip-relay persistent Python backend.
# fmlip-relay spawns a persistent server process; CREST communicates with it
# over a local TCP socket, avoiding repeated Python startup overhead.
# Output: crest_best.xyz (optimised Ar4 geometry near the LJ minimum ~3.82 Å)

command -v crest >/dev/null 2>&1 || {
  echo >&2 "Cannot find crest binary."
  exit 1
}

# Check for the fmlip-relay server; suggest pip install if absent
if ! command -v fmlip-relay-server >/dev/null 2>&1; then
  echo >&2 ""
  echo >&2 "ERROR: 'fmlip-relay-server' not found."
  echo >&2 "Install fmlip-relay from the CREST subproject directory:"
  echo >&2 ""
  echo >&2 "    pip install ../../subprojects/fmlip_relay"
  echo >&2 ""
  echo >&2 "For a user-local install add '--user', or activate a virtual environment first."
  exit 1
fi

# --- TOML run ---
# (No CLI equivalent; mlip settings are TOML-only.)
crest input.toml
