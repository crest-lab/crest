# Unit testing
set(
  tests
  "tblite"
  "gfnff"
  "gfn0"
  "gfn0occ"
  "optimization"
)
set(
  test-srcs
  "testmol.f90"
  "main.f90"
)
foreach(t IN LISTS tests)
  string(MAKE_C_IDENTIFIER ${t} t) 
  list(APPEND test-srcs "test_${t}.F90")
endforeach()

add_executable(
  "${PROJECT_NAME}-tester"
  "${test-srcs}"
)
target_link_libraries(
  "${PROJECT_NAME}-tester"
  PRIVATE
  "lib-${PROJECT_NAME}-static"
  "test-drive::test-drive"
)

foreach(t IN LISTS tests)
  add_test("${PROJECT_NAME}/${t}" "${PROJECT_NAME}-tester" "${t}")
endforeach()
