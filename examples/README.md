# Example applications of the CREST program

This directory contains examples covering the most common workflows
of the `crest` program.

Each example directory contains an input structure (`struc.xyz` or
similar), a shell script `run.sh`, and a TOML input file `input.toml`.

## Running an example

Go to the example directory and execute the script:
```bash
cd expl-6
./run.sh
```

The `run.sh` scripts show CLI usage.  Alternatively, every example can
be run through its TOML input file:
```bash
crest input.toml
```

TOML files are detected automatically by their `.toml` extension.  They
offer the same settings as the CLI flags but in a structured, documented
format that is easier to modify and reuse.

It is assumed that the `crest` binary is available in `$PATH`.


## Examples

| # | Topic | Molecule |
|---|-------|---------|
| **0** | *Dry run* — print settings without computing | 1-propanol |
| **1** | Single-point energy | 1-propanol |
| **2** | Geometry optimization | 1-propanol |
| **3** | Optimization + Hessian (vibrational frequencies) | 1-propanol |
| **4** | Standalone MD simulation | 1-propanol |
| **5** | Default iMTD-GC conformer search | 1-propanol |
| **6** | Two-level conformer search (GFN2//GFN-FF, A//B) | 1-propanol |
| **7** | iMTD-GC with ALPB implicit solvation (GFN2) | 1-propanol |
| **8** | Quick iMTD-GC conformer search (with -finalhess) | 1-propanol |
| **9** | Standalone CREGEN ensemble sorting | 1-propanol |
| **10** | Constrained conformer search | 1-propanol |
| **11** | Ensemble optimization (mdopt) | 1-propanol |
| **12** | NCI sampling mode (iMTD-NCI) | water trimer |
| **13** | Protonation site sampling | uracil |
| **14** | Metal/ion adducts (Cs+) | alpha-D-glucose |
| **15** | Tautomer screening | guanine |
| **16** | fmlip-relay: geometry optimisation with LJ potential | Ar4 cluster |
| **17** | fmlip-relay: geometry optimisation with FairChem UMA model | caffeine |
| **18** | fmlip-relay: geometry optimisation with MACE-OFF23 model | caffeine |
| **19** | Implicit-solvation add-on (ddX/EEQ-BC composite, GFN2 parent) | 1-propanol |
| **20** | External ORCA subprocess as engrad backend (geometry opt.) | n-pentane |
| **21** | ONIOM embedding (GFN2 head group in GFN-FF, opt. + freq.) | 1-propanol |
