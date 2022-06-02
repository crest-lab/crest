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

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    '-Wno-maybe-uninitialized',
    '-Wno-uninitialized',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_link_arguments(
    '-Wl,--allow-multiple-definition',
    language: 'fortran',
  )
  add_project_arguments(
    '-traceback',
    language: 'fortran',
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

omp_dep = dependency('openmp', required: fc.get_id() != 'intel')
if not omp_dep.found()
  omp_dep = declare_dependency(
    compile_args: '-qopenmp',
    link_args: '-qopenmp',
  )
endif
exe_deps += omp_dep

la_backend = get_option('la_backend')
if la_backend == 'mkl' or la_backend == 'mkl-static'
  add_project_arguments('-DWITH_MKL', language: 'fortran')
  if la_backend == 'mkl-static'
    add_project_link_arguments('-static', language: 'fortran')
    add_project_link_arguments('-static', language: 'c') # icc will do linking
  endif

  if fc.get_id() == 'gcc'
    libmkl_exe = [fc.find_library('mkl_gf_lp64')]
    libmkl_exe += fc.find_library('mkl_gnu_thread')
  elif fc.get_id() == 'intel'
    libmkl_exe = [fc.find_library('mkl_intel_lp64')]
    libmkl_exe += fc.find_library('mkl_intel_thread')
  elif fc.get_id() == 'pgi'
    libmkl_exe = [fc.find_library('mkl_intel_lp64')]
    libmkl_exe += fc.find_library('mkl_pgi_thread')
  endif
  libmkl_exe += fc.find_library('mkl_core')
  exe_deps += libmkl_exe

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
