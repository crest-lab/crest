
tests = [
  'molreader',
]

test_srcs = files(
  'main.f90',
  'mctc_accuracy.f90',
  'mctc_error.f90',
  'mctc_testing.f90', 
  'testsuite_data.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'cresttest',
  sources: test_srcs + srcs,
  dependencies: [crest_dep],
  link_language: 'fortran',
)

foreach t : tests
  test(t, tester, args: t)
endforeach
