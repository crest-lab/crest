#!/bin/bash
# ONIOM multi-level embedding: the alcohol head group of 1-propanol at
# GFN2-xTB ("high"), embedded in the GFN-FF description of the full
# molecule ("low"), with a link H atom saturating the cut C-C bond.
# Geometry optimization + vibrational frequencies on the ONIOM surface.
# Output: crestopt.xyz (optimized structure), vibspectrum (frequencies),
#         fragment.N.xyz (the generated ONIOM model systems)

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- TOML run ---
# The ONIOM setup has no CLI flags; it is configured through the
# [lwoniom] block of the TOML input.
crest input.toml
