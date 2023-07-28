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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

set(
  module-dir
  "${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
)
set(module-dir "${module-dir}" PARENT_SCOPE)


# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified."
  )
endif()


### Options and defaults

include("${CMAKE_CURRENT_SOURCE_DIR}/modules/crest-utils.cmake")
set(fortran_minimal_versions "GNU;7.5" "Intel;19.0")
check_minimal_compiler_version("Fortran" "${fortran_minimal_versions}")

option(WITH_OpenMP "Enable support for shared memory parallelisation with OpenMP" TRUE)

option(WITH_TBLITE "Enable build with the lightweight tight-binding library" FALSE)

option(WITH_TOMLF "Enable build with toml-f support" TRUE)

option(WITH_GFN0 "Enable build with GFN0-xTB support" TRUE)

option(WITH_GFNFF "Enable build with GFN-FF support" TRUE)

if(NOT DEFINED "${PROJECT_NAME}-dependency-method")
  set(
    "${PROJECT_NAME}-dependency-method"
    "subproject" "cmake" "pkgconf" "fetch"
  )
endif()

#
# Compiler settings
#
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#  set(dialect "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fbacktrace")
  set(dialect "-g -O0 -fbacktrace -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#  set(dialect "-axAVX2 -r8 -traceback")
  set(dialect "-g -O2 -r8 -align array64byte -traceback") 
  set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(dialect "-Mbackslash -Mallocatable=03 -r8 -traceback")
endif()
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}" PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}" PARENT_SCOPE)


#
# Populate crest_metadata.fh
#
set(version ${PROJECT_VERSION}${SOVERSION})
execute_process(COMMAND git show -s --format=%h
  RESULT_VARIABLE git_return
  OUTPUT_VARIABLE commit
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(git_return)
  set(commit "unknown-commit")
endif()
string(TIMESTAMP date "%a, %d %B %H:%M:%S, %m/%d/%Y")
set(author $ENV{USERNAME})
execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
set(origin ${HOSTNAME})
set(fcid ${CMAKE_Fortran_COMPILER_ID})
set(fcver ${CMAKE_Fortran_COMPILER_VERSION})
set(ccid ${CMAKE_C_COMPILER_ID})
set(ccver ${CMAKE_C_COMPILER_VERSION})
set(bsystem "cmake ${CMAKE_VERSION}")
set(tomlfvar "${WITH_TOMLF}")#string(TOLOWER "${MY_BOOL}" ${WITH_TOMLF}))
set(gfn0var "${WITH_GFN0}")
set(gfnffvar "${WITH_GFNFF}")
set(tblitevar "${WITH_TBLITE}")
configure_file(
  "${PROJECT_SOURCE_DIR}/assets/template/metadata.f90"
  "${PROJECT_BINARY_DIR}/crest_metadata.fh"
  @ONLY
)
# gfortran needs the file at another place...
configure_file(
  "${PROJECT_SOURCE_DIR}/assets/template/metadata.f90"
  "${PROJECT_BINARY_DIR}/include/crest_metadata.fh"
  @ONLY
)



