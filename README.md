# Conformer-Rotamer Ensemble Sampling Tool

[![Latest Version](https://img.shields.io/github/v/release/crest-lab/crest)](https://github.com/crest-lab/crest/releases/latest)
[![DOI](https://img.shields.io/badge/DOI-10.1039%2Fc9cp06869d%20-blue)](http://dx.doi.org/10.1039/c9cp06869d)

This is the offical repository of the CREST program developed by the Grimme group in Bonn.

CREST is an extension to the [`xtb`](https://github.com/grimme-lab/xtb) program.
It functions as an IO based OMP scheduler (*i.e.*, calculations are
performed by the `xtb` program) and tool for the creation and analysation of
structure ensembles.

<div align="center">
<img src="./assets/crest.png" alt="CREST" width="200">
</div>


## Documentation

The CREST documentation with installation instructions and application examples is hosted at <https://crest-lab.github.io/crest-docs/>.


---

## Installation quick guide

For any installation make sure that you have correctly installed and sourced the [`xtb`](https://github.com/grimme-lab/xtb) program before attempting any calculations with CREST.

There are multiple possible ways of installing CREST. 
For building the program from source we recommend the Intel `ifort` and `icc` compilers (tested with the 2021 version).

Detailed build instructions can be found at <https://crest-lab.github.io/crest-docs/page/installation>.

For builds with subprojects see [here](./subprojects/README.md).


### Precompiled binaries

To use the statically linked binaries (Intel compilers)
that can be found at the [release page](https://github.com/crest-lab/crest/releases),
of this repository.
The most recent program version is automatically build (`meson`/`ifort`) from the main branch and can be found at the [continous release page](https://github.com/crest-lab/crest/releases/tag/latest).
Simply unpack the binary and add it to your *PATH* variable.
```bash
unzip crest.zip
```
or
```bash
tar -xf crest-latest.tar.xz
```
The program should be directly executable.

### Tested builds
Working and tested builds of CREST (mostly on Ubuntu 20.04 LTS):

| Build System | Compiler | Linear Algebra Backend | Build type     | Status     |
|--------------|----------|------------------------|----------------|------------|
| CMake | GNU (gcc 10.3.0)  | [OpenBLAS](https://github.com/xianyi/OpenBLAS) (with OpenMP) | dynamic | :white_check_mark: |
| CMake | GNU (gcc 10.3.0)  |  MKL shared (oneAPI 2023.1)| dynamic | :white_check_mark:
| Meson | [Intel (ifort 2021.9.0)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)   | MKL static (oneAPI 2023.1) | static  | :white_check_mark: |


### `meson`

For the setup an configuration of meson see also the [meson setup](https://github.com/grimme-lab/xtb/blob/master/meson/README.adoc) page hosted at the `xtb` repository.
The chain of commands to build CREST with meson is:

```bash
export FC=ifort CC=icc
meson setup _build --prefix=$PWD/_dist
meson install -C _build
```

The `meson` build of CREST is mainly focused on and tested with the Intel `ifort`/`icc` compilers.
When attempting to build with `gfortran` and `gcc`, add `-Dla_backend=mkl` to the meson setup command. Compatibility with the GNU compilers might be limited, however.

By default the `meson` build will create a **statically** linked binary.

### `cmake`

For the setup of CMake see also the [CMake setup](https://github.com/grimme-lab/xtb/blob/master/cmake/README.adoc) page hosted at the `xtb` repository.
Building CREST with CMake works with the following chain of commands:
```bash
export FC=gfortran CC=gcc
cmake -B _build -DCMAKE_BUILD_TYPE=Release
```
and then to build the CREST binary
```bash
make -C _build
```

The CMake build of CREST is focused on and tested with the GNU `gfortran`/`gcc` compilers. The Intel compilers could technically be used as well, but in our experience 

By default the `cmake` build will create a **dynamically** linked binary.

### Conda

A [conda-forge](https://github.com/conda-forge) feedstock is maintained at <https://github.com/conda-forge/crest-feedstock>.

Installing CREST from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

```
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Once the `conda-forge` channel has been enabled, CREST can be installed with `conda`:

```
conda install crest
```

The confa-forge distribution is based on a CMake/`gfortran` build. 

---

## Citations

1. P. Pracht, F. Bohle, S. Grimme, *Phys. Chem. Chem. Phys.*, **2020**, 22, 7169-7192.
  DOI: [10.1039/C9CP06869D](https://dx.doi.org/10.1039/C9CP06869D)

2. S. Grimme, *J. Chem. Theory Comput.*, **2019**, 155, 2847-2862.
  DOI: [10.1021/acs.jctc.9b00143](https://dx.doi.org/10.1021/acs.jctc.9b00143)

3. P. Pracht, S. Grimme, *Chem. Sci.*, **2021**, 12, 6551-6568.
  DOI: [10.1039/d1sc00621e](https://dx.doi.org/10.1039/d1sc00621e)

4. P. Pracht, C.A. Bauer, S. Grimme, *J. Comput. Chem.*, **2017**, *38*, 2618-2631. 
  DOI: [10.1002/jcc.24922](https://dx.doi.org/10.1002/jcc.24922)

5. S. Spicher, C. Plett, P. Pracht, A. Hansen, S. Grimme,  *J. Chem. Theory Comput.*, **2022**,
  *18*, 3174-3189. DOI: [10.1021/acs.jctc.2c00239](https://dx.doi.org/10.1021/acs.jctc.2c00239)

6. P. Pracht, C. Bannwarth, *J. Chem. Theory Comput.*, **2022**, *18 (10)*, 6370-6385. DOI: [10.1021/acs.jctc.2c00578](https://dx.doi.org/10.1021/acs.jctc.2c00578)

## License

CREST is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

CREST is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in CREST by you, as defined in the GNU Lesser General Public license, shall be licensed as above, without any additional terms or conditions
