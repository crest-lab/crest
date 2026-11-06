set(_lib "mstore")
set(_pkg "MSTORE")
set(_url "https://github.com/grimme-lab/mstore")
set(_branch "663245d739be0123da61c917e55116b0c3db4c74")

# Discovery method order can be overridden by the parent project, e.g.:
#   set(mstore_FIND_METHOD "subproject" "cmake")
if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

# Reuse whichever utils macro your main project already provides.
# Replace "crest-utils" with the actual name if yours differs.
include("${CMAKE_CURRENT_LIST_DIR}/crest-utils.cmake")

crest_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_branch}")

set(found FALSE)
if(TARGET "mstore::mstore")
  set(found TRUE)
endif()
message(STATUS "Found mstore: ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
