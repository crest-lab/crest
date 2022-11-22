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

cmake_minimum_required(VERSION 3.14)

project(
  "crest"
  LANGUAGES "Fortran" "C"
  VERSION "3.0"
)


# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
set(libs)
add_subdirectory("config")


#
# Libraries
#

# OpenMP dependency
if(NOT TARGET OpenMP::OpenMP_Fortran AND WITH_OpenMP)
  find_package(OpenMP REQUIRED)
  list(
    APPEND libs
    OpenMP::OpenMP_Fortran
  )
endif()

# BLAS and LAPACK
if(NOT TARGET BLAS::BLAS)
  find_package(BLAS REQUIRED)
  if(NOT TARGET BLAS::BLAS AND BLAS_FOUND)
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    target_link_libraries(BLAS::BLAS INTERFACE "${BLAS_LIBRARIES}")
    target_link_options(BLAS::BLAS INTERFACE "${BLAS_LINKER_FLAGS}")
  endif()
endif()
if(NOT TARGET LAPACK::LAPACK)
  find_package(LAPACK REQUIRED)
  if(NOT TARGET LAPACK::LAPACK AND LAPACK_FOUND)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    target_link_libraries(LAPACK::LAPACK INTERFACE "${LAPACK_LIBRARIES}")
    target_link_options(LAPACK::LAPACK INTERFACE "${LAPACK_LINKER_FLAGS}")
  endif()
endif()
list(
  APPEND libs
  LAPACK::LAPACK
  BLAS::BLAS
)

## mctc-lib
#if(NOT TARGET "mctc-lib::mctc-lib")
#  find_package("mctc-lib" REQUIRED)
#  list(
#    APPEND libs
#    mctc-lib::mctc-lib
#  )
#endif()

# tblite
if(WITH_TBLITE)
  find_package("tblite" REQUIRED)
  add_definitions(-DWITH_TBLITE)
  list(
    APPEND libs
    tblite::tblite
  )
endif()


# GFN0-xTB
if(WITH_GFN0)
  find_package("gfn0" REQUIRED)
  add_definitions(-DWITH_GFN0)
  list(
    APPEND libs
    gfn0::gfn0
  )
endif()



#
# CREST sources
#
set(prog)
set(srcs)
add_subdirectory("src")


#
# Executables
#
add_executable(
  ${PROJECT_NAME}-exe
  "${srcs}" "${prog}"
)
set_target_properties(
  ${PROJECT_NAME}-exe
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_link_libraries(
  ${PROJECT_NAME}-exe
  PRIVATE
  "${libs}"
)


#
# Binary installing option
#
install(
  TARGETS
  "${PROJECT_NAME}-exe"
  EXPORT
  "${PROJECT_NAME}-targets"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
