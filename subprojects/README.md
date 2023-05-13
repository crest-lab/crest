## Using subprojects

Newer versions of CREST use external projects:
- [toml-f](https://github.com/toml-f/toml-f), a TOML parser for Fortran
- [gfn0](https://github.com/pprcht/gfn0), a GFN0-xTB standalone library
- [gfnff](https://github.com/pprcht/gfnff), a GFN-FF standalone library
- [tblite](https://github.com/tblite/tblite), a lightweight implementation of the GFN1 and GFN2-xTB Hamiltonians


Both `cmake` and `meson` should be **able to handle the download automatically** (with meson being a little bit better at this).

However, some of the projects are also set up as `git` submodules, in case you want to download the most current commits by hand.
To do so, in the CREST main directory use
```bash
git submodule init
git submodule update
```
Which should download all the subprojects.

To update the submodule sources from the respective remote branches
```bash
git submodule update --remote
```
can be used.

---

An important exception that requires special attention is the `tblite` subproject when building with `meson`. 
To make the build work, some of the build meson instructions of `tblite` must be updated.
We have prepared patch files for this located at [packagefiles/tblite/](./packagefile/tblite/)
After downloading the `tblite` subproject with the above `git` commands, change to the directory, and apply the patches via
```bash
cd subprojects/tblite
git apply ../packagefile/tblite/tblite_patch.patch
```

