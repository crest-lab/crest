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
  VERSION "2.11.0"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
set(exe-deps)
add_subdirectory("config")

# OpenMP dependency
if(NOT TARGET OpenMP::OpenMP_Fortran)
  find_package(OpenMP REQUIRED)
  list(
    APPEND exe-deps
    OpenMP::OpenMP_Fortran
  )
endif()

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
  APPEND exe-deps
  LAPACK::LAPACK
  BLAS::BLAS
)

# General configuration information
set(prog)
set(srcs)
add_subdirectory("src")

add_executable(
  "${PROJECT_NAME}-exe"
  "${srcs}" "${prog}"
)
set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${exe-deps}"
)

install(
  TARGETS
  "${PROJECT_NAME}-exe"
  EXPORT
  "${PROJECT_NAME}-targets"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
