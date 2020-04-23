#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst -mdopt xtb.trj
 else
    $crst -mdopt xtb.trj -xnam $xtbin
 fi


# A ensemble file (or MD trajectory) can also
# be optimized in a standalone application
# of CREST using the '-mdopt' flag.
# The optimized structures are written to a
# file called 'crest_ensemble.xyz', but will
# not be sorted in any way.
