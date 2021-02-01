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

find_package(LAPACK REQUIRED)
list(
  APPEND exe-deps
  LAPACK::LAPACK
)

set(exe-deps "${exe-deps}" PARENT_SCOPE)

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

install(
  PROGRAMS
  "${CMAKE_CURRENT_SOURCE_DIR}/combi-scripts/crest_combi"
  DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
