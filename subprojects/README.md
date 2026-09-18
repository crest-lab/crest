## Using CREST subprojects

Newer versions of CREST use external projects. Pinned versions/commits below
reflect the current `subprojects/*.wrap` revisions (and matching git submodule
checkouts).

**Direct subprojects** (toggled via a CREST build option, all enabled by default):

| Library | Description | Build Option | Pinned version | git submodule | CMake build | `meson` build |
| ------- | ----------- | ------------ | ------------- | :-----------: | :---------: | :-----------: |
| [`tblite`](https://github.com/tblite/tblite) | A lightweight implementation of the GFN1 and GFN2-xTB Hamiltonians | `-DWITH_TBLITE=true` (default) | `v0.7.0` (`b244a0f`) | ✅ | ✅ | ✅ |
| [`toml-f`](https://github.com/toml-f/toml-f) | A TOML parser for Fortran | `-DWITH_TOMLF=true` (default) | `v0.5.2` (`28f4601`) | ✅ | ✅ | ✅ |
| [`gfn0`](https://github.com/pprcht/gfn0) | A GFN0-xTB standalone library | `-DWITH_GFN0=true` (default) | `v0.1.1` (`d77fea8`) | ✅ | ✅ | ✅ |
| [`gfnff`](https://github.com/pprcht/gfnff) | A GFN-FF standalone library | `-DWITH_GFNFF=true` (default) | `v0.2.0` (`272ed91`) | ✅ | ✅ | ✅ |
| [`libpvol`](https://github.com/pprcht/libpvol) | Molecular volume/surface library (PV, XHCFF) | `-DWITH_LIBPVOL=true` (default) | branch `build-update` (`010bdda`) | ✅ | ✅ | ✅ |
| [`lwoniom`](https://github.com/crest-lab/lwoniom) | A lightweight ONIOM implementation | `-DWITH_LWONIOM=true` (default) | `v0.0.1` (`ab66c7e`) | ✅ | ✅ | ✅ |
| [`fmlip_relay`](https://github.com/pprcht/fmlip-relay) | Relay interface to the fmlip ML/IP potential | `-DWITH_FMLIP_RELAY=true` (default) | branch `main` (`1f29707`) | ✅ | ✅ | ✅ |

**Transitive dependencies** (pulled in automatically by the above — mostly
through `tblite` — with no separate CREST build option):

| Library | Description | Pulled in via | Pinned version | git submodule | CMake build | `meson` build |
| ------- | ----------- | ------------- | ------------- | :-----------: | :---------: | :-----------: |
| [`mctc-lib`](https://github.com/grimme-lab/mctc-lib) | Modular computation tool chain library (I/O, data types) | `tblite` | `v0.5.2` (`e9de066`) | ✅ | ✅ | ✅ |
| [`multicharge`](https://github.com/grimme-lab/multicharge) | Electronegativity-equilibration partial charges | `tblite` | `v0.5.0` (`6a5d63f`) | ✅ | ✅ | ✅ |
| [`dftd4`](https://github.com/dftd4/dftd4) | DFT-D4 dispersion correction | `tblite` | `v4.2.0` (`6e1f59c`) | ✅ | ✅ | ✅ |
| [`s-dftd3`](https://github.com/dftd3/simple-dftd3) | Simple DFT-D3 dispersion correction | `tblite` | `v1.4.0` (`6f0b06f`) | ✅ | ✅ | ✅ |
| [`ddx`](https://github.com/ddsolvation/ddX) | Domain-decomposition implicit solvation | `tblite` | `v0.8.0` (`4d79e3d`) | ✅ | ✅ | ✅ |
| [`jonquil`](https://github.com/toml-f/jonquil) | JSON support for TOML Fortran | `mctc-lib` | `v0.3.0` (`4d43ffe`) | ❌ (wrap-redirect) | ✅ | ✅ |
| [`mstore`](https://github.com/grimme-lab/mstore) | Molecular structure store (unit-test data) | `tblite` (tests) | `v0.3.0` (`663245d`) | ✅ | ✅ | ✅ |
| [`test-drive`](https://github.com/fortran-lang/test-drive) | Lightweight unit-testing framework | tests | `v0.6.1` (`c506771`) | ✅ | ✅ | ✅ |

> Note: `libpvol` and `fmlip_relay` have no tagged releases yet, so they track a
> branch tip (commit shown). `jonquil` is not a CREST submodule — its version is
> whatever `mctc-lib` pins via a wrap-redirect.


Both `cmake` and `meson` should be **able to handle the download automatically** (with meson being a little bit better at this). The build option can be specified in the respective setup step.

However, some projects are also set up as `git` submodules (see table) if you want to download the most current commits by hand.
To do so, in the CREST main directory use
```bash
git submodule init
git submodule update
```
which should check out all the subprojects.

To update the submodule sources from the respective remote branches
```bash
git submodule update --remote
```
can be used.

Alternatively, a source directory of the respective project can be placed in the subprojects directory, or a symbolic link can be set.
