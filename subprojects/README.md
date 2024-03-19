## Using CREST subprojects

Newer versions of CREST use external projects:

| Library | Description | Build Option | git submodule | CMake build | `meson` build |
| ------- | ----------- | ------------ | :-----------: | :---------: | :-----------: |
| [`toml-f`](https://github.com/toml-f/toml-f) | A TOML parser for Fortran | `-DWITH_TOMLF=true` (default) | ✅ | ✅ | ✅ |
| [`gfn0`](https://github.com/pprcht/gfn0) | A GFN0-xTB standalone library | `-DWITH_GFN0=true` (default) | ✅ | ✅ | ✅ |
| [`gfnff`](https://github.com/pprcht/gfnff) | A GFN-FF standalone library | `-DWITH_GFNFF=true` (default) | ✅ | ✅ | ✅ |
| [`tblite`](https://github.com/tblite/tblite) | A lightweight implementation of the GFN1 and GFN2-xTB Hamiltonians | `-DWITH_TBLITE=true` (default) | ✅ | ✅ | ✅ |
| [`lwoniom`](https://github.com/crest-lab/lwoniom) | A lightweight ONIOM implementation  | `-DWITH_LWONIOM=true` (default) | ✅ | ✅ | ✅ |
<!--
| [`xhcff`](https://github.com/zellerf/xhcff) | Implementation of the XHCFF force field | `-DWITH_XHCFF=true` | ✅ | ✅ | ✅ |
-->


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
