#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -tautomerize -ewin 10.0
 else
    $crst struc.xyz -tautomerize -ewin 10.0 -xnam $xtbin
 fi

# The -tautomerize flag can be used to request
# a screening of prototropic tautomers.
# The structures are build from a sequence of
# protonating and deprotonating steps of the
# (neutral) input structure.
# In the above example this procedure is
# performed on the guanine molecule to get
# the gas phase tautomers at GFN2-xTB level.
# Within the 10 kcal/mol window 5 tautomers
# should remain at this level.
# The structures can be found in 'tautomers.xyz'.

