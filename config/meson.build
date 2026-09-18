# This file is part of crest.
# SPDX-Identifier: LGPL-3.0-or-later
#
# Included from the root meson.build via subdir('config').
# Runs in the same variable scope as the root — all variables set here
# are visible in the root after this file returns.
#
# Responsibilities:
#   - Compiler identification and friendly display names
#   - Per-compiler build flags (release and debug)
#   - Static-binary link arguments
#   - Metadata header generation (crest_metadata.fh)
#   - inc_dirs (build include path carrying the generated header)

# ═══════════════════════════════════════════════════════════════════════════════
#  Compiler identification
# ═══════════════════════════════════════════════════════════════════════════════
fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

fc_id = fc.get_id()   # 'gcc' | 'intel' | 'intel-llvm' | ...
cc_id = cc.get_id()

_fc_name_map = {
  'gcc'        : 'GNU Fortran',
  'intel'      : 'Intel Fortran (Classic)',
  'intel-llvm' : 'Intel Fortran (LLVM/ifx)',
  'flang'      : 'Flang',
  'flang-new'  : 'Flang-new',
}
_cc_name_map = {
  'gcc'        : 'GNU C',
  'intel'      : 'Intel C (Classic)',
  'intel-llvm' : 'Intel C (LLVM/icx)',
  'clang'      : 'Clang',
}

fcid_str = _fc_name_map.get(fc_id, fc_id)
ccid_str = _cc_name_map.get(cc_id, cc_id)

if fc_id not in ['gcc', 'intel', 'intel-llvm']
  warning(
    'Unsupported Fortran compiler "@0@". Proceeding, but results may be wrong.'
      .format(fc_id),
  )
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Compiler flags
# ═══════════════════════════════════════════════════════════════════════════════
_is_debug = get_option('buildtype') in ['debug', 'debugoptimized']

# --- GNU (gfortran / gcc) -----------------------------------------------------
if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language : 'fortran',
  )
  if _is_debug
    add_project_arguments(
      '-fcheck=all',
      '-ffpe-trap=invalid,zero,overflow',
      '-fbounds-check',
      '-finit-real=snan',
      '-finit-integer=-999',
      '-Wall',
      '-Wextra',
      language : 'fortran',
    )
  endif
endif

# --- Intel classic (ifort / icc) -----------------------------------------------
if fc_id == 'intel'
  add_project_arguments(
    '-r8',                    # treat all REAL as REAL(8)
    '-align', 'array64byte',
    '-traceback',
    language : 'fortran',
  )
  if _is_debug
    add_project_arguments(
      '-check', 'all',
      '-fpe0',
      language : 'fortran',
    )
  endif
endif

# --- Intel LLVM (ifx / icx) ----------------------------------------------------
if fc_id == 'intel-llvm'
  add_project_arguments(
    '-r8',
    '-align', 'array64byte',
    language : 'fortran',
  )
  if _is_debug
    add_project_arguments(
      '-traceback',
      '-check', 'all',
      '-fpe0',
      language : 'fortran',
    )
  endif
endif

static_build = get_option('static')

# Static subproject linking: selects static archives for subproject deps
# (tblite, gfn0, toml-f, …) without the fully-static `-static` syslib linking.
# A fully static build (-Dstatic) implies it.
link_deps_static = static_build or get_option('static-deps')

# Keyword set for dependency() lookups of subprojects/external libraries.
# Passed via kwargs: so that a disabled static-deps adds no 'static' keyword at
# all.  An explicit dependency(static : false) would force
# default_library=shared onto every fallback subproject, overriding a global
# -Ddefault_library=static (needed e.g. by package builds that link the
# bundled subprojects statically but take system libraries as shared objects).
static_dep_kw = link_deps_static ? {'static' : true} : {}

if fc_id == 'intel-llvm'
  ## -lifport  : POSIX intrinsics (getcwd, chdir, …) not auto-linked by ifx.
  ##             In static builds -static-intel already bundles libifport.a, so
  ##             adding a bare -lifport here would cause the linker to resolve it
  ##             to the shared library (before -static-intel takes effect).
  ## libgcc_s  : was previously added to supply _Unwind_* symbols for libifcore's
  ##             signal handler. No longer needed since signal registration was
  ##             moved to crest_install_signal() via ISO_C_BINDING, bypassing
  ##             libifcore's handler entirely.
  _libgcc_s = run_command('gcc', '-print-file-name=libgcc_s.so',
                           check : false).stdout().strip()
  _gcc_s_link = (_libgcc_s != 'libgcc_s.so') ? [_libgcc_s] : []
  if _gcc_s_link.length() == 0
    warning('intel-llvm: libgcc_s.so not found — _Unwind_* symbols will be NULL, ' +
            'Intel Fortran runtime will SIGSEGV at startup.')
  endif
  # For static builds libgcc_s.so is not needed (signal handler uses ISO_C_BINDING
  # and no longer goes through libifcore's _Unwind_* path), and adding an explicit
  # .so path would defeat -static.  Only add it for dynamic builds.
  _lifport_arg = static_build ? [] : ['-lifport']
  _gcc_s_arg   = static_build ? [] : _gcc_s_link
  add_project_link_arguments(_lifport_arg + _gcc_s_arg, language : 'fortran')
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Static-binary link arguments
# ═══════════════════════════════════════════════════════════════════════════════
# Injected globally so they propagate to every link step (library + executable).
# See BUILD.md for notes on required system .a archives.

if static_build
  if fc_id == 'gcc'
    # Full static: bundles glibc, libgfortran, libgomp, libopenblas/lapack.
    # Requires libgfortran-static and glibc-static system packages.
    add_project_link_arguments('-static', language : 'fortran')
    add_project_link_arguments('-static', language : 'c')

  elif fc_id == 'intel'
    # -static-intel  : statically links Intel runtime (ifcore, imf, svml)
    # -qopenmp-link=static : statically links Intel OpenMP (libiomp5)
    # -static        : forces all -l flags (including MKL) to resolve to .a
    add_project_link_arguments(
      '-static-intel',
      '-qopenmp-link=static',
      '-static',
      language : 'fortran',
    )

  elif fc_id == 'intel-llvm'
    # -static only: with LIBRARY_PATH including the Intel compiler lib dir,
    # the linker finds libifcore.a/libimf.a/libiomp5.a directly.
    # -static-intel + -qopenmp-link=static are NOT used because combined with
    # -fopenmp link args from non-Intel subprojects they cause libifcore.a vs
    # libifcoremt.a duplicate-symbol errors (ifx 2025.3).  Allow-multiple-
    # definition suppresses the remaining conflict (same pattern as ifort above).
    add_project_link_arguments(
      '-static',
      '-Wl,--allow-multiple-definition',
      language : 'fortran',
    )
  endif
endif

# ═══════════════════════════════════════════════════════════════════════════════
#  Metadata: fill assets/template/metadata.f90 → <builddir>/crest_metadata.fh
# ═══════════════════════════════════════════════════════════════════════════════
_git = find_program('git', required : false)
_commit = 'unknown-commit'
if _git.found()
  _r = run_command(_git, 'show', '-s', '--format=%h', check : false)
  if _r.returncode() == 0
    _commit = _r.stdout().strip()
  endif
endif

_py = find_program('python3', 'python', required : false)
_date = 'unknown-date'
if _py.found()
  _r = run_command(
    _py, '-c',
    'import datetime; print(datetime.datetime.now().strftime("%a, %d %B %H:%M:%S, %m/%d/%Y"))',
    check : false,
  )
  if _r.returncode() == 0
    _date = _r.stdout().strip()
  endif
endif

_hostname = 'unknown-host'
_hostname_prog = find_program('hostname', required : false)
if _hostname_prog.found()
  _r = run_command(_hostname_prog, check : false)
  if _r.returncode() == 0
    _hostname = _r.stdout().strip()
  endif
endif

_sh = find_program('sh', required : false)
_user = 'unknown'
if _sh.found()
  _r = run_command(
    _sh, '-c', 'printf "%s" "${USER:-${USERNAME:-unknown}}"',
    check : false,
  )
  if _r.returncode() == 0 and _r.stdout().strip() != ''
    _user = _r.stdout().strip()
  endif
endif

_conf = configuration_data()
_conf.set('version',    meson.project_version())
_conf.set('commit',     _commit)
_conf.set('date',       _date)
_conf.set('author',     _user)
_conf.set('origin',     _hostname)
_conf.set('fcid',       fcid_str)
_conf.set('fcver',      fc.version())
_conf.set('ccid',       ccid_str)
_conf.set('ccver',      cc.version())
_conf.set('bsystem',    'meson ' + meson.version())
# The with_* booleans are set in the root after optional deps are resolved;
# use string placeholders here and let the root patch them in if needed.
# In practice configure_file runs at configure-time so these are all known
# by the time the root calls subdir('config') — see note below.
#
# NOTE: the with_* variables (with_tomlf, with_tblite, …) are set in the root
# meson.build AFTER this file runs, so they are not yet available here.
# The root therefore calls configure_file() itself using this _conf object
# after appending the feature flags.  We expose _conf so the root can do that.
_conf.set('tomlfvar',   'false')
_conf.set('gfn0var',    'false')
_conf.set('gfnffvar',   'false')
_conf.set('tblitevar',  'false')
_conf.set('libpvolvar', 'false')
_conf.set('lwoniomvar', 'false')
_conf.set('fmliprelayvar', 'false')

# Expose _conf to the root so it can set the with_* keys and then call
# configure_file().  The root owns the final configure_file() call.
metadata_conf = _conf


