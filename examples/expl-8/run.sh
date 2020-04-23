#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -protonate -swel Cs+
 else
    $crst struc.xyz -protonate -swel Cs+ -xnam $xtbin
 fi

# In a modified version of the protonation tool
# other ionization adducts can be created
# (only mono nuclear ions)
# To do this, the flag '-swel' (short for switch element)
# is used to indicate the new ion and its charge,
# e.g., Na+, Ca2+, Li+, etc.
#
# As a example the alpha-D-glucose-Cs+ adducts
# will be created at the GFN2-xTB level with the above command.
# The adducts can be found in the file 'protonated.xyz'

