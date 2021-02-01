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

For any installation make sure that you have correctly installed and sourced the [`xtb`](https://github.com/grimme-lab/xtb) program before attempting any calculations with `crest`.

There are multiple possible ways of installing `crest`. 
For building the program from source we recommend the Intel `ifort` and `icc` compilers.


### The easy way

Use the statically linked binaries (Intel Compiler 2019.6.324)
that can be found at the [release page](https://github.com/grimme-lab/crest/releases),
of this repository.
Simply unpack the binary and add it to your *PATH* variable.
```bash
tar -xvzf crest.tgz
```
The program should be directly executable in most cases.

### Building via `make`

In the `src` directory a `Makefile` can be found to build a statically linked binary. Modify the makefile to your requirements and build the program via
```bash
make
```

### Cmake

For the setup of Cmake see also the [Cmake setup](https://github.com/grimme-lab/xtb/blob/master/cmake/README.adoc) page hosted at the `xtb` repository.
Building `crest` with CMake works with the following chain of commands:

```bash
export FC=ifort CC=icc
cmake -B _build_intel -DCMAKE_BUILD_TYPE=Release
make -C _build_intel
```

To install the `crest` binaries to `/usr/local` use (might require `sudo`)

```bash
make -C _build_intel install
```

### Meson

For the setup an configuration of meson see also the [meson setup](https://github.com/grimme-lab/xtb/blob/master/meson/README.adoc) page hosted at the `xtb` repository.
The chain of commands to build `crest` with meson is:

```bash
export FC=ifort CC=icc
meson setup _build_intel --prefix=$PWD/_dist
meson install -C _build_intel
```


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

`crest` is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

`crest` is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in `crest` by you, as defined in the GNU Lesser General Public license, shall be licensed as above, without any additional terms or conditions
