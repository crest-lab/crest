#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -protonate
 else
    $crst struc.xyz -protonate -xnam $xtbin
 fi


# This command will create protomers of the uracil molecule.
# The default energy window for this application is 30 kcal/mol
# Only 3 structures should remain in the gas phase at the
# default GFN2-xTB level.
# The structures can be found in the file 'protonated.xyz'
