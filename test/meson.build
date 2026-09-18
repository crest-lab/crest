# This file is part of crest.
# SPDX-Identifier: LGPL-3.0-or-later

tests = [
  'tblite',
  'gfnff',
  'gfn0',
  'gfn0occ',
  'CN',
  'optimization',
  'molecular_dynamics',
  'metadynamics',
  'pbc_cregen',
  'getsym',
  'irmsd',
  'graphs',
  'ddx',
  'external',
  'oniom',
]

test_srcs = files(
  'testmol.f90',
  'helpers.f90',
  'main.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.F90'.format(t))
endforeach

crest_tester = executable(
  'crest-tester',
  sources             : test_srcs,
  include_directories : inc_dirs,
  link_with           : lib_crest,
  dependencies        : [testdrive_dep, lapack_dep, blas_dep]
                        + _omp_link_dep + _optional_deps,
  link_language       : 'fortran',
)

foreach t : tests
  test(
    'crest/@0@'.format(t),
    crest_tester,
    args    : [t],
    timeout : 120,
  )
endforeach
