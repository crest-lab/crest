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

set(_lib "toml-f")
set(_pkg "TOML-F")
set(_url "https://github.com/toml-f/toml-f")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  set("${_pkg}_FIND_METHOD" "subproject" "cmake" "fetch" "pkgconf")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/crest-utils.cmake")

crest_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}")

if(TARGET "toml-f::toml-f")
  set (found TRUE)
else()
  set (found FALSE)
endif()
message(STATUS "Found toml-f: ${found}")

unset(_lib)
unset(_pkg)
unset(_url)
