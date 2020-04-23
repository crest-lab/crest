#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -dry
 else
    $crst struc.xyz -dry -xnam $xtbin
 fi

# Before starting any calculation, settings
# can be checked with the '-dry' flag.
# This will only print a summary about the
# selected settings and thresholds to the
# consol and check for the xtb binary.
#
# Every time the input file (struc.xyz) is
# something else than 'coord', a file called
# 'coord' will be (over-)written, containing
# the atomic coordinates in Bohr. CREST will
# then continue to use and overwrite this
# coord file for all further calculations.
