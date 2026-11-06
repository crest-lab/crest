set(_lib "multicharge")
set(_pkg "MULTICHARGE")
set(_url "https://github.com/grimme-lab/multicharge")
set(_branch "6a5d63f9e9e29dcf13cc47cc27f33bf9015681bf")

# Discovery method order can be overridden by the parent project, e.g.:
#   set(multicharge_FIND_METHOD "subproject" "cmake")
if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

# Reuse whichever utils macro your main project already provides.
# Replace "crest-utils" with the actual name if yours differs.
include("${CMAKE_CURRENT_LIST_DIR}/crest-utils.cmake")

crest_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_branch}")

set(found FALSE)
if(TARGET "multicharge::multicharge")
  set(found TRUE)
endif()
message(STATUS "Found multicharge: ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
