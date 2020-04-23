#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

$crst struc.xyz -cregen xtb.trj -ewin 100.0


# The sorting routine from the CREST conformational search can be
# used as a standalone to sort any .xyz or .trj ensemble file.
# The above command will sort the file xtb.trj according to
# its energy and determine duplicate structures.
# Two files are written analogous to 'crest_conformers.xyz'
# and 'crest_rotamers.xyz'.
# The new file 'crest_ensemble.xyz' will contain only unique
# structures from xtb.trj, while the new file 'xtb.trj.sorted'
# is just a sorted version of the original file (without the
# -ewin flag the default 6.0 kcal/mol window will be used)
# The routine requires a reference structure which is given
# by 'struc.xyz'.

