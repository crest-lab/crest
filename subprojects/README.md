## Using CREST subprojects

Newer versions of CREST use external projects:

| Library | Description | Build Option | git submodule |
| ------- | ----------- | ------------ | :-----------: |
| [`toml-f`](https://github.com/toml-f/toml-f) | A TOML parser for Fortran | `-DWITH_TOMLF=true` (default) | &#9745; |
| [`gfn0`](https://github.com/pprcht/gfn0) | A GFN0-xTB standalone library | `-DWITH_GFN0=true` (default) | &#9745; |
| [`gfnff`](https://github.com/pprcht/gfnff) | A GFN-FF standalone library | `-DWITH_GFNFF=true` | &#9745; |
| [`tblite`](https://github.com/tblite/tblite) | A lightweight implementation of the GFN1 and GFN2-xTB Hamiltonians | `-DWITH_TBLITE=true` | &#9746; |

Both `cmake` and `meson` should be **able to handle the download automatically** (with meson being a little bit better at this). The build option can be specified in the respective setup step.

However, some projects are also set up as `git` submodules (see table) if you want to download the most current commits by hand.
To do so, in the CREST main directory use
```bash
git submodule init
git submodule update
```
which should download all the subprojects.

To update the submodule sources from the respective remote branches
```bash
git submodule update --remote
```
can be used.

---

### `tblite` additional information
The [`tblite`](https://github.com/tblite/tblite) subproject is an important exception.
It is **not** set up as a `git` submodule, but it could still be downloaded manually.
To do so, while in the CREST main directory, use the usual
```bash
git clone https://github.com/tblite/tblite.git subprojects/tblite
```
to clone `tblite` to the correct place.
However, to make the build work after downloading, some of the `meson` build instructions of `tblite` must be updated.
We have prepared a patch file for this located at [packagefiles/tblite/](./packagefile/tblite/)
Change to the directory and apply the patches via
```bash
cd subprojects/tblite
git apply ../packagefile/tblite/tblite_patch.patch
```

Note that `meson` will download and apply the patch automatically *if you don't download `tblite` yourself*!
