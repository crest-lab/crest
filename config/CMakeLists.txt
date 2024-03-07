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
# Set the module path for CMake includes
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules;${CMAKE_MODULE_PATH}")

# Options for enabling or disabling features
option(STATICLINK "Build a static binary" ON)
option(WITH_OPENMP "Enable OpenMP support" ON)
option(WITH_TBLITE "Enable support for tblite" OFF)
option(WITH_TOMLF "Enable support for toml-f" OFF)
option(WITH_GFN0 "Enable support for GFN0-xTB" OFF)
option(WITH_GFNFF "Enable support for GFN-FF" OFF)
option(WITH_XHCFF "Enable support for XHCFF" OFF)
option(WITH_LWONIOM "Enable support for lwONIOM" OFF)

#########################################################################################
#########################################################################################

# Compiler settings for GNU and Intel Fortran compilers
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-g -O0 -fbacktrace -ffree-line-length-none -fbacktrace")
  set(bounds "-fbounds-check")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(dialect "-g -O2 -r8 -align array64byte -traceback")
  set(bounds "-check bounds")
endif()

# Apply the compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#########################################################################################
#########################################################################################

# Handling dependencies and features

# Find LAPACK & BLAS
if(STATICLINK)
  set(BLA_STATIC ON)
endif()
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  list(APPEND libs ${LAPACK_LIBRARIES})
endif()
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  list(APPEND libs ${BLAS_LIBRARIES})
endif()

# OpenMP
if(WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_Fortran_FOUND)
    #list(APPEND libs OpenMP::OpenMP_Fortran)
    list(APPEND omp ${OpenMP_Fortran_LIBRARIES})
  endif()
endif()

# tblite
if(WITH_TBLITE)
  find_package(tblite REQUIRED)
  if(tblite_FOUND)
    add_definitions(-DWITH_TBLITE)
    include_directories(${tblite_INCLUDE_DIRS})
    list(APPEND libs ${tblite_LIBRARIES})
  endif()
endif()

# toml-f
if(WITH_TOMLF)
  find_package(toml-f REQUIRED)
  if(toml-f_FOUND)
    add_definitions(-DWITH_TOMLF)
    include_directories(${toml-f_INCLUDE_DIRS})
    list(APPEND libs ${tomlf_LIBRARIES})
  endif()
endif()

# GFN0-xTB
if(WITH_GFN0)
  find_package(GFN0 REQUIRED)
  if(GFN0_FOUND)
    add_definitions(-DWITH_GFN0)
    include_directories(${GFN0_INCLUDE_DIRS})
    list(APPEND libs ${GFN0_LIBRARIES})
  endif()
endif()

# GFN-FF
if(WITH_GFNFF)
  find_package("gfnff" REQUIRED)
  if(gfnff_FOUND)
    add_definitions(-DWITH_GFNFF)
    list(APPEND libs gfnff::gfnff)
  endif()
endif()

# lwONIOM
if(WITH_LWONIOM)
  find_package(LWONIOM REQUIRED)
  if(LWONIOM_FOUND)
    add_definitions(-DWITH_LWONIOM)
    include_directories(${LWONIOM_INCLUDE_DIRS})
    list(APPEND libs ${LWONIOM_LIBRARIES})
  endif()
endif()

# XHCFF
if(WITH_XHCFF)
  find_package(XHCFF REQUIRED)
  if(XHCFF_FOUND)
    add_definitions(-DWITH_XHCFF)
    include_directories(${XHCFF_INCLUDE_DIRS})
    list(APPEND libs ${XHCFF_LIBRARIES})
  endif()
endif()

# IMPORTANT
set(omp ${omp} PARENT_SCOPE)
set(libs ${libs} PARENT_SCOPE)

##########################################################################################
##########################################################################################

# Populate crest_metadata.fh
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
set(tomlfvar "${WITH_TOMLF}")
set(gfn0var "${WITH_GFN0}")
set(gfnffvar "${WITH_GFNFF}")
set(tblitevar "${WITH_TBLITE}")
set(xhcffvar "${WITH_XHCFF}")
set(lwoniomvar "${WITH_LWONIOM}")

configure_file(
  "${PROJECT_SOURCE_DIR}/assets/template/metadata.f90"
  "${PROJECT_BINARY_DIR}/crest_metadata.fh"
  @ONLY
)

#########################################################################################
#########################################################################################
