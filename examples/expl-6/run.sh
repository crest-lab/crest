#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -nci
 else
    $crst struc.xyz -nci -xnam $xtbin
 fi


# This will execute the NCI sampling mode of CREST of the
# water trimer with default settings.
# A wall-potential is automatically generated and added to
# the calculation to prevent dissociation.
# The NCI mode is a special case of the constrained sampling.
# Just like the regular conformational search unique conformers
# can be found in the file 'crest_conformers.xyz'.
# All degenerate conformers (rotamers, pseudo-enantiomers)
# can be found in the file 'crest_rotamers.xyz'


