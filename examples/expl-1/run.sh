#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -ewin 2.0
 else
    $crst struc.xyz -ewin 2.0 -xnam $xtbin
 fi


# This will execute a conformational search with default settings
# for the 1-propanol molecule.
# The energy window is set to 2.0 kcal/mol with the '-ewin' flag
# (instead of the default 6.0 kcal/mol window)
# Within this window there should be 4 conformers for 1-propanol
# in the gas phase.
# The 4 unique conformers can be found in the file 'crest_conformers.xyz'.
# All degenerate conformers (rotamers, pseudo-enantiomers) of the 4 structures
# can be found in the file 'crest_rotamers.xyz'


