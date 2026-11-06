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

set(_lib "dftd4")
set(_pkg "DFTD4")
set(_url "https://github.com/dftd4/dftd4")
set(_branch "6e1f59c3f39d919a2dbef0601d2576727c8b30e8")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/crest-utils.cmake")

set(temp_with_tests ${WITH_TESTS}) # Save the current value of WITH_TESTS
set(WITH_TESTS FALSE CACHE BOOL "Temporarily disable tests for the dftd4 subproject" FORCE)
set(WITH_API FALSE)
crest_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_branch}")

set(found FALSE)
if(TARGET "dftd4::dftd4")
  set (found TRUE)
endif()
message(STATUS "Found dftd4: ${found}")

set(WITH_TESTS ${temp_with_tests} CACHE BOOL "Enable tests for the main project" FORCE)

unset(_lib)
unset(_pkg)
unset(_url)
