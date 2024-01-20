# This file is part of crest.
# SPDX-Identifier: LGPL-3.0-or-later
#
# crest is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# crest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with crest.  If not, see <https://www.gnu.org/licenses/>.

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')
fc_id = fc.get_id()
cc_id = cc.get_id()

if fc.get_id() != cc.get_id()
  warning('FC and CC are not from the same vendor')
endif

## ========================================= ##
## Compiler specific default arguments
## ========================================= ##  
if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    '-Wno-maybe-uninitialized',
    '-Wno-uninitialized',
    '-Wno-unused-variable',
    '-Wno-unused-dummy-argument',
    '-Wno-unused-function',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_link_arguments(
    '-Wl,--allow-multiple-definition',
    language: 'fortran',
  )
  add_global_link_arguments(
    '-Wl,--allow-multiple-definition',
    language: 'fortran',
  )
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
  add_project_arguments(
    '-DLINUX',
    language: 'c',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'flang'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    language: 'fortran',
  )
endif

add_project_arguments('-D_Float128=__float128', language: 'c')


## ========================================= ##
## build type option arguments
## ========================================= ## 
if ( get_option('default_library') == 'static')
  message('Static linking selected')
  add_project_link_arguments('-static', language: 'fortran')
  add_project_link_arguments('-static', language: 'c') # icc will do linking
endif


## ========================================= ##
## OpenMP
## ========================================= ##
if get_option('openmp')
  omp_dep = dependency('openmp', required: fc.get_id() != 'intel')
  if not omp_dep.found()
    omp_dep = declare_dependency(
      compile_args: '-qopenmp',
      link_args: '-fopenmp',
    )
  endif
  exe_deps += omp_dep
endif

## ======================================== ##
## Linear Algebra Libraries
## ======================================== ##
la_backend = get_option('la_backend')
message('Linear algebra backend: '+get_option('la_backend'))
if la_backend == 'mkl' or la_backend == 'mkl-static'
  add_project_arguments('-DWITH_MKL', language: 'fortran')
  if la_backend == 'mkl-static'
    add_project_link_arguments('-static', language: 'fortran')
    add_project_link_arguments('-static', language: 'c') # icc will do linking
  endif

  if get_option('default_library') == 'shared'
    mkl_rt_dep = cc.find_library('mkl_rt', required: true)
    exe_deps += mkl_rt_dep
  else
    if fc.get_id() == 'gcc'
      libmkl_exe = [cc.find_library('mkl_gf_lp64')]
      if get_option('openmp')
        libmkl_exe += cc.find_library('mkl_gnu_thread')
      endif
    elif fc.get_id() == 'intel' or fc.get_id() == 'intel-cl'
      libmkl_exe = [cc.find_library('mkl_intel_lp64')]
      if get_option('openmp')
        libmkl_exe += cc.find_library('mkl_intel_thread')
      endif
    elif fc.get_id() == 'pgi' or fc.get_id() == 'nvidia_hpc'
      libmkl_exe = [cc.find_library('mkl_intel_lp64')]
      if get_option('openmp')
        libmkl_exe += cc.find_library('mkl_pgi_thread')
      endif
    endif
    if not get_option('openmp')
      libmkl_exe += cc.find_library('mkl_sequential')
    endif
    libmkl_exe += cc.find_library('mkl_core')
    exe_deps += libmkl_exe
  endif

elif la_backend == 'mkl-rt'
  add_project_arguments('-DWITH_MKL', language: 'fortran')

  mkl_rt_dep = fc.find_library('mkl_rt', required: true)
  if fc.get_id() == 'intel'
    exe_deps += fc.find_library('ifcore')
  endif

  exe_deps += mkl_rt_dep

elif la_backend == 'openblas'
  # search for OpenBLAS
  blas_dep = dependency('openblas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('openblas', required: true)
  endif
  exe_deps += blas_dep
  # some OpenBLAS versions can provide lapack, check if we can find dsygvd
  openblas_provides_lapack = fc.links(
    'external dsygvd; call dsygvd(); end',
    dependencies: blas_dep,
  )
  # otherwise we fall back to LAPACK
  if not openblas_provides_lapack
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack', required: true)
    endif
    exe_deps += lapack_dep
  endif

elif la_backend == 'custom'
  foreach lib: get_option('custom_libraries')
    exe_deps += fc.find_library(lib)
  endforeach

else
  # Find BLAS (usually netlib, but in conda also OpenBLAS/MKL)
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas', required: true)
  endif
  exe_deps += blas_dep
  # Find LAPACK (usually netlib, but in conda also MKL)
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack', required: true)
  endif
  exe_deps += lapack_dep
endif


## ========================================= ##
## External subprojects
## ========================================= ##  

# GFN1-xTB and GFN2-xTB via TBLITE
if get_option('WITH_TBLITE')
 add_project_arguments('-DWITH_TBLITE', language: 'fortran')
 tblite_dep = dependency(
   'tblite',
   version: '>=0.2',
   fallback: ['tblite', 'tblite_dep'],
   default_options: ['default_library=static', 'api=false'],
 )
 exe_deps += tblite_dep
endif

# TOML-F
if get_option('WITH_TOMLF')
 add_project_arguments('-DWITH_TOMLF', language: 'fortran')
 tomlf_dep = dependency(
   'toml-f',
   version: '>=0.2.0',
   fallback: ['toml-f', 'tomlf_dep'],
   default_options: ['default_library=static'],
 )
 exe_deps += tomlf_dep
endif

# GFN0-xTB
if get_option('WITH_GFN0')
 add_project_arguments('-DWITH_GFN0', language: 'fortran')
 gfn0_dep = dependency(
   'gfn0',
#   version: '>=0.2',
   fallback: ['gfn0', 'gfn0_dep'],
   default_options: ['default_library=static','with_gbsa=true'],
 )
 exe_deps += gfn0_dep
endif

# GFN-FF
if get_option('WITH_GFNFF')
 add_project_arguments('-DWITH_GFNFF', language: 'fortran')
 gfnff_dep = dependency(
   'gfnff',
   fallback: ['gfnff', 'gfnff_dep'],
   default_options: ['default_library=static','with_gbsa=true'],
 )
 exe_deps += gfnff_dep
endif


# XHCFF
if get_option('WITH_XHCFF')
 add_project_arguments('-DWITH_XHCFF', language: 'fortran')
 xhcff_dep = dependency(
   'xhcff',
   fallback: ['xhcff', 'xhcfflib_dep'],
   default_options: ['default_library=static'],
 )
 exe_deps += xhcff_dep
endif


# lwONIOM
if get_option('WITH_LWONIOM')
 add_project_arguments('-DWITH_LWONIOM', language: 'fortran')
 lwoniom_dep = dependency(
   'lwoniom',
   fallback: ['lwoniom', 'lwoniom_dep'],
   default_options: ['default_library=static'],
 )
 exe_deps += lwoniom_dep
endif


## ========================================= ##
## populate the data for crest_metadata.fh
## ========================================= ##  
commit = 'unknown commit'
git = find_program('git', required: false)
if git.found()
  git_commit = run_command(git, 'show', '-s', '--format=%h',check:true)
  if git_commit.returncode() == 0
    commit = git_commit.stdout().strip()
  endif
endif
# create configuration data
config = configuration_data({
  'name': meson.project_name(),
  'version': meson.project_version(),
  'description': 'Conformer Rotamer Ensemble Sampling Tool',
  'commit': commit,
  'date': run_command('date',check:true).stdout().     strip(),
  'author': run_command('id','-u','-n', check:true).stdout().strip(),
  'origin': run_command('hostname', check:true).stdout().strip(),
  'fcid': fc.get_id(),
  'fcver': fc.version(),
  'ccid': cc.get_id(),
  'ccver': cc.version(),
  'bsystem': 'meson '+meson.version(),
  'tomlfvar': get_option('WITH_TOMLF'),
  'gfn0var': get_option('WITH_GFN0'),
  'gfnffvar': get_option('WITH_GFNFF'), 
  'tblitevar': get_option('WITH_TBLITE'),
  'xhcffvar': get_option('WITH_XHCFF'),
  'lwoniomvar': get_option('WITH_LWONIOM'),
})


