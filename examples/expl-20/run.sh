#!/bin/bash
# Geometry optimization of n-pentane using an external ORCA subprocess as the
# energy+gradient backend (r2scan-3c, assembled from the TOML input).
#
# CREST runs its ANCOPT optimizer and calls ORCA for every engrad evaluation.
# The ORCA calculator is configured entirely through the TOML input (method,
# cores and memory); no ORCA template file is required.
#
# NOTE: the ORCA executable path in input.toml is a placeholder
# (/PATH/TO/ORCA/orca) and must be set to your local ORCA installation for
# this example to actually run.

command -v crest >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary."; exit 1; }

# --- TOML run ---
crest input.toml
