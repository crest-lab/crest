# Example applications of the CREST program

This directory contains several examples for 
standard applications of the `crest` program.

Each example directory contains a input structure
(typically called `struc.xyz`) and a bash script
called `run.sh` that includes some information about
the example and will execute the calculation upon
execution.

To run the example scripts simply go to the respective
directory and execute it from the command line:
```bash
./run.sh
```

It is assumed that the `xtb` and `crest` binaries
are present in the *PATH* variable as such.
The `run.sh` scripts will check for this, however.


## Examples

0. *dry run* of the `crest` program
1. default conformational search (iMTD-GC)
2. example for different CMD settings
3. sorting an ensemble file (CREGEN)
4. constrained conformational sampling
5. standalone optimization along a trajectory
6. NCI sampling mode (iMTD-NCI)
7. protonation site sampling
8. modified protonation site sampling
9. tautomer sampling
