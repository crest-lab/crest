#!/bin/bash

xtbin='xtb'
crst='crest'

command -v $xtbin >/dev/null 2>&1 || { echo >&2 "Cannot find xtb binary. Exit."; exit 1; }
command -v $crst >/dev/null 2>&1 || { echo >&2 "Cannot find crest binary. Exit."; exit 1; }

if [ $xtbin == 'xtb' ]
 then
    $crst struc.xyz -ewin 2.0 -quick -prop ohess
 else
    $crst struc.xyz -ewin 2.0 -quick -prop ohess -xnam $xtbin
 fi

# Some further calculations can be added automatically
# after the conformer search with the '-prop' command.
# In the above example, after searching for the
# conformers of 1-propanol, each conformer is optimized
# again and frequencies are calculated (ohess).
# The conformer ensemble is then re-ranked with free
# energies from RRHO contributions.
#
# There are also some different 'quick'-modes to run
# the conformational search with reduced settings.
# With these modes the conformational space will be
# explored less extensively, but it will speed up
# the calculation. ('-quick','-squick','-mquick')

