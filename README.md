# Conformer-Rotamer Ensemble Sampling Tool

[![Latest Version](https://img.shields.io/github/v/release/grimme-lab/crest)](https://github.com/grimme-lab/crest/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.1039%2Fc9cp06869d%20-blue)](http://dx.doi.org/10.1039/c9cp06869d)

This is the offical repository of the `crest` program developed by the Grimme group in Bonn.

`crest` is an extension to the [`xtb`](https://github.com/grimme-lab/xtb) program.
It functions as an IO based OMP scheduler (*i.e.*, calculations are
performed by the `xtb` program) and tool for the creation and analysation of
structure ensembles.

<div align="center">
<img src="./assets/crest.png" alt="CREST" width="200">
</div>

## Installation

We are currently preparing the source code and will be providing it in this repository
with one of the upcoming releases.

Until then please use the statically linked binaries (Intel Compiler 17.0.7)
that can be found at the [release page](https://github.com/grimme-lab/crest/releases),
of this repository.
Simply unpack the binary and add it to your *PATH* variable.
```bash
tar -xvzf crest.tgz
```

Also make sure that you have correctly installed and sourced the [`xtb`](https://github.com/grimme-lab/xtb) program before attempting any calculations with `crest`.

## Examples

This repository contains several examples for default applications of `crest`.

See [`examples`](examples). To test the installation please try to run [example 0](examples/expl-0).

## Documentation

The `crest` documentation is hosted at [read-the-docs](https://xtb-docs.readthedocs.io/en/latest/crest.html).

## Citations

1. P. Pracht, F. Bohle, S. Grimme, *Phys. Chem. Chem. Phys.*, **2020**, 22, 7169-7192.
  DOI: [10.1039/C9CP06869D](https://dx.doi.org/10.1039/C9CP06869D)

2. S. Grimme, *J. Chem. Theory Comput.*, **2019**, 155, 2847-2862.
  DOI: [10.1021/acs.jctc.9b00143](https://dx.doi.org/10.1021/acs.jctc.9b00143)


## License

`crest` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.

