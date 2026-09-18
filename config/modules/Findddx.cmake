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

set(_lib "ddx")
set(_pkg "DDX")
set(_url "https://github.com/ddsolvation/ddX")
set(_branch "4d79e3d9caeae5e602683572a71cb550414f9b09")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/crest-utils.cmake")

# ddX builds its examples, tests and driver by default; only the library is
# needed here.  The options are forced off for the subproject build.
set(EXAMPLES OFF CACHE BOOL "Disable ddx examples" FORCE)
set(TESTS OFF CACHE BOOL "Disable ddx tests" FORCE)
set(WARNING_FLAGS OFF CACHE BOOL "Disable ddx warning flags" FORCE)
crest_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_branch}")

set(found FALSE)
if(TARGET "ddx::ddx")
  set (found TRUE)
endif()
message(STATUS "Found ddx: ${found}")

# ddX's own CMake writes its Fortran .mod files into <binary>/subprojects/ddx/src
# but only exposes an install-interface include dir ($<INSTALL_INTERFACE:include>).
# In-tree consumers (tblite, crest) therefore cannot locate ddx.mod during the
# build.  Expose the module output directory through the build interface so the
# subproject build resolves like the sibling subprojects (mctc-lib, s-dftd3).
if(TARGET "ddx")
  target_include_directories("ddx" INTERFACE
    "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/subprojects/ddx/src>")
endif()

unset(_lib)
unset(_pkg)
unset(_url)
unset(_branch)
