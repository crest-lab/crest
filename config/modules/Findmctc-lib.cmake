set(_lib "mctc-lib")
set(_pkg "MCTCLIB")
set(_url "https://github.com/grimme-lab/mctc-lib")
set(_branch "e9de066d89f250d1cfb6de3a33f0c27c0e2f855d")

# Discovery method order can be overridden by the parent project, e.g.:
#   set(mctc-lib_FIND_METHOD "subproject" "cmake")
if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

# Reuse whichever utils macro your main project already provides.
# Replace "crest-utils" with the actual name if yours differs.
include("${CMAKE_CURRENT_LIST_DIR}/crest-utils.cmake")

crest_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_branch}")

set(found FALSE)
if(TARGET "mctc-lib::mctc-lib")
  set(found TRUE)
endif()
message(STATUS "Found mctc-lib: ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
