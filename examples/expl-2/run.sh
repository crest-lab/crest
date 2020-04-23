#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -ewin 2.0 -g h2o -gfn2 -T 4
 else
    $crst struc.xyz -ewin 2.0 -g h2o -gfn2 -T 4 -xnam $xtbin
 fi


# This will execute a conformational search with some manually changed
# settings for the 1-propanol molecule.
# The GBSA implicit solvation model for H2O is employed with the
# '-g' flag.
# Furthermore, the use of GFN2-xTB is requested explicitly ('-gfn2')
# and the program is ordered to use 4 CPU threads ('-T').
# For a 1-propanol the conformers in implicit solvation are the
# same as in the gas phase, but the relative energies should
# differ significantly.
# Unique conformers can be found in the file 'crest_conformers.xyz'.
# All degenerate conformers (rotamers, pseudo-enantiomers)
# can be found in the file 'crest_rotamers.xyz'


