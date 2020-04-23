#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -constrain 1-4
    $crst struc.xyz -cinp .xcontrol.sample
 else
    $crst struc.xyz -constrain 1-4
    $crst struc.xyz -cinp .xcontorl.sample -xnam $xtbin
 fi


# Constraint conformational sampling is possible by
# providing the constrainment info as a file
# via the '-cinp' flag.
# For detailed information about the constraining
# options see the online documentation of
# CREST and xTB.
# However, a dummy constraining file '.xcontrol.sample'
# can be written by CREST with a seperate call using
# The '-constrain <atoms>' flag.
# In the above example the carbon atoms and the oxygen
# atom of 1-propanol (atoms 1-4) will be constrained.
# In the resulting "ensemble" only conformers resulting
# from different OH dihedral angles will be present
# (2 conformers total)
