# Building crest with Meson

## Prerequisites

| Tool | Minimum version | Notes |
|------|----------------|-------|
| [Meson](https://mesonbuild.com) | 0.57.0 | `pip install meson` or distro package |
| [Ninja](https://ninja-build.org) | 1.10 | usually installed alongside meson |
| Fortran compiler | — | gfortran ≥ 10, ifort ≥ 2021, or ifx ≥ 2023 |
| C compiler | — | gcc or the matching Intel C compiler |
| LAPACK + BLAS | — | OpenBLAS, Intel MKL, or netlib |
| OpenMP | — | libgomp (GNU) or libomp/libiomp5 (Intel) |

Optional (auto-detected or pulled from `subprojects/`):

- tblite, toml-f, GFN-FF, GFN0-xTB, libpvol, lwONIOM, fmlip-relay
- test-drive (only for `--tests`)

---

## Quick start

```sh
# Configure (defaults: release build, OpenBLAS/auto LAPACK, OpenMP on, all
# optional libs auto-detected, unit tests enabled)
meson setup build

# Compile
ninja -C build

# Run tests
ninja -C build test

# Install to /usr/local
ninja -C build install
```

---

## Compiler selection

The build system auto-detects whatever Fortran/C compilers are first on
`$PATH`. Use the native-file templates in `config/` to pin a specific
toolchain:

```sh
# Pure GNU (gfortran + gcc)
meson setup build --native-file config/gnu.ini

# Intel LLVM (ifx + icx) — oneAPI 2023+
source /opt/intel/oneapi/setvars.sh
meson setup build --native-file config/intel-llvm.ini

# Intel classic (ifort + icc) — oneAPI 2022 or earlier
meson setup build --native-file config/intel-classic.ini

# Mixed: Intel Fortran (ifx) + GNU C (gcc)
meson setup build --native-file config/intel-fortran-gnu-c.ini
```

Or set compilers directly via environment variables (older style):

```sh
FC=gfortran CC=gcc meson setup build
FC=ifort    CC=icc meson setup build
FC=ifx      CC=icx meson setup build
```

---

## Build options

Pass options with `-D` at configure time, or modify with `meson configure`:

| Option | Default | Description |
|--------|---------|-------------|
| `openmp` | `true` | Enable OpenMP parallelisation |
| `lapack` | `auto` | LAPACK/BLAS provider: `auto`, `openblas`, `mkl`, `netlib`, `custom` |
| `lapack_libs` | `[]` | Library names for `lapack=custom` |
| `blas_libs` | `[]` | Library names for `lapack=custom` |
| `static` | `false` | Link a fully static binary |
| `tblite` | `auto` | tblite semiempirical library |
| `toml-f` | `auto` | TOML-Fortran (file-based input) |
| `gfn0` | `auto` | GFN0-xTB library |
| `gfnff` | `auto` | GFN-FF library |
| `libpvol` | `auto` | libpvol (volume computation) |
| `lwoniom` | `auto` | lwONIOM |
| `fmlip-relay` | `auto` | fmlip-relay ML/IP interface |
| `tests` | `true` | Build unit tests |

Feature options (`auto` / `enabled` / `disabled`):

- `auto` — use it if found, silently skip if not
- `enabled` — require it; fail the build if not found
- `disabled` — never use it even if installed

Examples:

```sh
# Disable all optional libraries (minimal build)
meson setup build -Dtblite=disabled -Dtoml-f=disabled \
                  -Dgfn0=disabled -Dgfnff=disabled \
                  -Dlibpvol=disabled -Dlwoniom=disabled \
                  -Dfmlip-relay=disabled

# Require tblite (fail if not found)
meson setup build -Dtblite=enabled

# Debug build with bounds checking
meson setup build --buildtype=debug

# Change an option after configuration
meson configure build -Dopenmp=false
```

---

## LAPACK / BLAS selection

### Auto (default)

- For Intel compilers: tries MKL first, then OpenBLAS, then netlib
- For GNU: tries OpenBLAS first, then MKL (with `mkl_gnu_thread`), then netlib

```sh
meson setup build -Dlapack=auto   # (this is the default)
```

### OpenBLAS

```sh
meson setup build -Dlapack=openblas
```

OpenBLAS bundles both BLAS and LAPACK. Make sure `libopenblas-dev` (Debian/Ubuntu)
or `openblas-devel` (RHEL/Fedora) is installed, or that `pkg-config --exists openblas`
succeeds.

### Intel MKL

```sh
source /opt/intel/oneapi/setvars.sh
meson setup build -Dlapack=mkl
```

The build system selects the correct threading layer automatically:

| Fortran compiler | MKL threading layer | OpenMP runtime |
|-----------------|---------------------|----------------|
| gfortran | `mkl_gnu_thread` | libgomp |
| ifort / ifx | `mkl_intel_thread` | libiomp5 / libomp |

**Do not mix** `mkl_gnu_thread` with Intel OpenMP or vice versa — this causes
silent wrong results or crashes.

If `pkg-config` can see `mkl-sdl`, that single-dynamic-library interface is
used instead and no threading-layer selection is needed.

### Custom libraries

For non-standard LAPACK installations (e.g. a vendor-tuned LAPACK on a
cluster module):

```sh
meson setup build -Dlapack=custom \
  -Dlapack_libs=lapack,blas \
  -Dblas_libs=blas
# or with full paths via pkg-config / LIBRARY_PATH
```

---

## Fully static binary

A static binary embeds all libraries including the OpenMP runtime and LAPACK.
This is the most portable output for distribution on HPC clusters.

### GNU static

Requires: `libgfortran.a`, `libgomp.a`, `libopenblas.a` (or `liblapack.a` +
`libblas.a`) to be available as static `.a` archives.  On Debian/Ubuntu
install `gfortran-static`, `libgomp1` (usually comes with `libgomp-staticdev`),
and `libopenblas-dev`.

```sh
meson setup build_static      \
  --buildtype=release         \
  --native-file config/gnu.ini \
  -Dstatic=true               \
  -Dlapack=openblas
ninja -C build_static
# Result: build_static/crest  — fully self-contained ELF
ldd build_static/crest  # should print "not a dynamic executable"
```

### Intel classic static

```sh
source /opt/intel/oneapi/setvars.sh
meson setup build_static               \
  --native-file config/intel-classic.ini \
  -Dstatic=true                         \
  -Dlapack=mkl
ninja -C build_static
```

Intel's `-static-intel -qopenmp-link=static` flags are applied automatically.
The Intel static libraries (`libifcore.a`, `libimf.a`, `libsvml.a`,
`libiomp5.a`) must be present — they are typically in
`$ONEAPI_ROOT/compiler/latest/linux/compiler/lib/intel64_lin/`.

### Intel LLVM (ifx) static

```sh
source /opt/intel/oneapi/setvars.sh
meson setup build_static              \
  --native-file config/intel-llvm.ini  \
  -Dstatic=true                       \
  -Dlapack=mkl
ninja -C build_static
```

\`-static-intel -qopenmp-link=static\` are applied automatically for ifx.

---

## Subprojects

Optional chemistry libraries are resolved in this order:

1. **System-installed** — found via `pkg-config` or in standard library paths
2. **Git submodule** — if `subprojects/<name>/` exists and contains a
   `meson.build` (populate with `git submodule update --init --recursive`)
3. **Wrap file** — `subprojects/<name>.wrap` tells Meson to clone the repo on
   first use with `meson subprojects download`

To pre-fetch all wrap-defined subprojects:

```sh
meson subprojects download
```

To update existing subproject clones:

```sh
meson subprojects update
```

---

## Metadata generation

At configure time Meson fills `assets/template/metadata.f90` and writes the
result to `<builddir>/crest_metadata.fh`.  The placeholders populated are:

| Placeholder | Value |
|-------------|-------|
| `@version@` | project version from `meson.build` |
| `@commit@`  | short git hash (or `unknown-commit`) |
| `@date@`    | configure timestamp |
| `@author@`  | `$USER` / `$USERNAME` |
| `@origin@`  | hostname |
| `@fcid@`    | Fortran compiler name |
| `@fcver@`   | Fortran compiler version |
| `@ccid@`    | C compiler name |
| `@ccver@`   | C compiler version |
| `@bsystem@` | `meson <version>` |
| `@tomlfvar@` | `true` / `false` |
| `@gfn0var@` | `true` / `false` |
| `@gfnffvar@` | `true` / `false` |
| `@tblitevar@` | `true` / `false` |
| `@libpvolvar@` | `true` / `false` |
| `@lwoniomvar@` | `true` / `false` |

---

## Compiler / LAPACK cross-compatibility reference

```
Fortran compiler │ C compiler │ Recommended LAPACK │ OpenMP runtime
─────────────────┼────────────┼────────────────────┼───────────────
gfortran         │ gcc        │ OpenBLAS (default) │ libgomp
gfortran         │ gcc        │ MKL                │ libgomp + mkl_gnu_thread
ifort / ifx      │ icc / icx  │ MKL  ← best match │ libiomp5 / libomp
ifort / ifx      │ gcc        │ MKL                │ libiomp5 (Intel wins link)
ifort / ifx      │ icc / icx  │ OpenBLAS (seq.)    │ libiomp5
```

**Rule of thumb:** the compiler that drives the *final link step* owns the
OpenMP runtime.  With mixed toolchains, the Fortran compiler always drives the
link step here (crest is a Fortran-primary project), so use the LAPACK
threading layer that matches the *Fortran* compiler.
