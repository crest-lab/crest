set(_lib "fmlip_relay")
set(_pkg "FMLIP_RELAY")
set(_url "https://github.com/pprcht/fmlip-relay")
set(_branch "06608dea8b7deb53f4a25630f6dcb8be0cdc6e76")

# Discovery method order can be overridden by the parent project, e.g.:
#   set(FMLIP_RELAY_FIND_METHOD "subproject" "cmake")
if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

# Reuse whichever utils macro your main project already provides.
# Replace "crest-utils" with the actual name if yours differs.
include("${CMAKE_CURRENT_LIST_DIR}/crest-utils.cmake")

crest_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_branch}")

set(found FALSE)
if(TARGET "fmlip_relay::fmlip_relay")
  set(found TRUE)
endif()
message(STATUS "Found fmlip_relay: ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
