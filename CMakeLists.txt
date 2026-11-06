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

cmake_minimum_required(VERSION 3.17)
option(WITH_OBJECT "To build using object library" TRUE)
option(INSTALL_MODULES "Install Fortran module files to include directory." FALSE)

# Buggy CMake versions
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.27.0 AND CMAKE_VERSION VERSION_LESS 3.28.0)
   set(WITH_OBJECT FALSE)
endif()

# Setup the crest Project
project(
  crest
  LANGUAGES "C" "CXX" "Fortran"
  VERSION  3.1.0
  DESCRIPTION "A tool for the exploration of low-energy chemical space"
)

# Apply debug flags when building in Debug mode
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} \
    -g -O0 -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999")

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# Include further configurations
add_subdirectory("config")


###############################################################################
####################### SUBPROJECTS & DEPENDENCIES ############################
###############################################################################
# Check a specific CMake targets and execute the
# corresponding Find scripts (located in config/modules/)

# LAPACK, BLAS, OpenMP (usually via shared dependencies)
if(STATICBUILD)
 set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
 set(CMAKE_EXE_LINKER_FLAGS "-static")
 set(CMAKE_LINK_SEARCH_START_STATIC TRUE)
 set(CMAKE_LINK_SEARCH_END_STATIC TRUE)
 set(BLA_STATIC ON)
endif()
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
  get_target_property(OpenMP_Fortran_LIBRARIES OpenMP::OpenMP_Fortran INTERFACE_LINK_LIBRARIES)
  message(STATUS "OpenMP::OpenMP_Fortran is linking the following libraries:")
  foreach(lib ${OpenMP_Fortran_LIBRARIES})
    message(STATUS "${lib}")
  endforeach()
endif()
# Propagate OpenMP link flags globally so subproject executables (e.g. pvol tests)
# that bypass CMake target dependencies still link the OpenMP runtime correctly.
if(WITH_OpenMP AND OpenMP_Fortran_FOUND)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# Check if we are using OpenBLAS (need a precompiler definition if yes)
# Match "libopenblas" specifically to avoid false positives from conda env
# path names that contain "openblas" (e.g. crest-gnu-openblas/lib/libgomp.so).
if(LAPACK_LIBRARIES OR BLAS_LIBRARIES)
  set(_all_blas_libs "${LAPACK_LIBRARIES};${BLAS_LIBRARIES}")
  string(FIND "${_all_blas_libs}" "libopenblas" _openblas_in_lapack)

  if(NOT _openblas_in_lapack EQUAL -1)
    message(STATUS "libopenblas was found as part of LAPACK/BLAS")
    add_compile_definitions(WITH_OPENBLAS)
  endif()
endif()

# Fortran unit test interface.  Imported unconditionally (not gated on
# WITH_TESTS): several subprojects build their test suites without a guard and
# pull in test-drive themselves.  By resolving the bundled subprojects/test-drive
# here first, the `test-drive::test-drive` target already exists when those
# subprojects look for it, so they reuse the local copy instead of falling back
# to fetching/downloading it.
if(NOT TARGET "test-drive::test-drive")
   find_package("test-drive" REQUIRED)
endif()

# TOML-F (needs to be imported before tblite)
if(NOT TARGET "toml-f::toml-f" AND WITH_TOMLF)
   find_package("toml-f" REQUIRED)
   add_compile_definitions(WITH_TOMLF)
endif()

# tblite
if(NOT TARGET "tblite::tblite" AND WITH_TBLITE)
   find_package("mctc-lib"    REQUIRED)
   find_package("mstore"      REQUIRED)
   find_package("multicharge" REQUIRED)
   find_package("dftd4"       REQUIRED)
   find_package("s-dftd3"     REQUIRED)
   find_package("ddx"         REQUIRED)
   find_package("tblite"      REQUIRED)
   add_compile_definitions(WITH_TBLITE)
   # ddX is a hard requirement of this block, so its modules are available to
   # crest sources here; WITH_DDX guards the standalone solvation components.
   add_compile_definitions(WITH_DDX)
endif()

# g-xTB via tblite
if(WITH_GXTB)
   if(NOT WITH_TBLITE)
      message(FATAL_ERROR "WITH_GXTB requires WITH_TBLITE to be enabled")
   endif()
   add_compile_definitions(WITH_GXTB)
endif()

# GFN-FF
if(NOT TARGET "gfnff::gfnff" AND WITH_GFNFF)
   find_package("gfnff" REQUIRED)
   add_compile_definitions(WITH_GFNFF)
endif()

# GFN0-xTB
if(NOT TARGET "gfn0::gfn0" AND WITH_GFN0)
   find_package("gfn0" REQUIRED)
   add_compile_definitions(WITH_GFN0)
endif()

# libpvol
if(NOT TARGET "pvol::pvol" AND WITH_LIBPVOL)
   find_package("libpvol" REQUIRED)
   add_compile_definitions(WITH_LIBPVOL)
endif()

# lwONIOM
if(NOT TARGET "lwoniom::lwoniom" AND WITH_LWONIOM)
   find_package("lwoniom" REQUIRED)
   add_compile_definitions(WITH_LWONIOM)
endif()

if(NOT TARGET "fmlip_relay:fmlip_relay" AND WITH_FMLIP_RELAY)
  find_package("fmlip_relay" REQUIRED)
  add_compile_definitions(WITH_FMLIP_RELAY)
endif()

# Sources: initialize program sources (prog) and library sources (srcs) empty
set(prog)
set(srcs)
add_subdirectory("src")


if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
   file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
endif()

###############################################################################
########################## OBJECT LIBRARY #####################################
###############################################################################
if(WITH_OBJECT AND NOT STATICBUILD)
   
   add_library(
      "${PROJECT_NAME}-object"
      OBJECT
      ${srcs}
   )

   # customize object library
   set_target_properties(
      "${PROJECT_NAME}-object"
      PROPERTIES
      Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
      POSITION_INDEPENDENT_CODE ON
   )

   # link object library conditionally against required
   target_link_libraries(
      "${PROJECT_NAME}-object"
      PUBLIC
      $<$<BOOL:${WITH_TBLITE}>:tblite::tblite>
      $<$<BOOL:${WITH_GFN0}>:gfn0::gfn0>
      $<$<BOOL:${WITH_GFNFF}>:gfnff::gfnff>
      $<$<BOOL:${WITH_LIBPVOL}>:pvol::pvol>
      $<$<BOOL:${WITH_TOMLF}>:toml-f::toml-f>
      $<$<BOOL:${WITH_LWONIOM}>:lwoniom::lwoniom>
      $<$<BOOL:${WITH_FMLIP_RELAY}>:fmlip_relay::fmlip_relay> 
      $<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>
   )

   # include directories
   target_include_directories(
      "${PROJECT_NAME}-object"
      PUBLIC
      ${crest-config-dir}
      ${PROJECT_BINARY_DIR}
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
   )

endif()

###############################################################################
############################### Static Library ################################
###############################################################################
if(WITH_OBJECT AND NOT STATICBUILD)
   add_library(
      "lib-${PROJECT_NAME}-static"
      STATIC
      $<TARGET_OBJECTS:${PROJECT_NAME}-object>
   )
else()
   add_library(
      "lib-${PROJECT_NAME}-static"
      STATIC
      ${srcs}
   )
endif()

# workaround
set(LINK_OpenMP FALSE)
if(WITH_OpenMP AND NOT STATICBUILD)
  set(LINK_OpenMP TRUE)
endif()

target_link_libraries(
   "lib-${PROJECT_NAME}-static"
   PUBLIC
   ${BLAS_LIBRARIES}
   ${LAPACK_LIBRARIES}
   $<$<BOOL:${LINK_OpenMP}>:OpenMP::OpenMP_Fortran>
   $<$<BOOL:${WITH_TBLITE}>:tblite::tblite>
   $<$<BOOL:${WITH_GFN0}>:gfn0::gfn0>
   $<$<BOOL:${WITH_GFNFF}>:gfnff::gfnff>
   $<$<BOOL:${WITH_LIBPVOL}>:pvol::pvol>
   $<$<BOOL:${WITH_TOMLF}>:toml-f::toml-f>
   $<$<BOOL:${WITH_LWONIOM}>:lwoniom::lwoniom>
   $<$<BOOL:${WITH_FMLIP_RELAY}>:fmlip_relay::fmlip_relay> 
   $<$<BOOL:${STATICBUILD}>:-static>
)


set_target_properties(
   "lib-${PROJECT_NAME}-static"
   PROPERTIES
   Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
   ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
   POSITION_INDEPENDENT_CODE ON
   OUTPUT_NAME "${PROJECT_NAME}"
 )


target_include_directories(
  "lib-${PROJECT_NAME}-static"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
)


####################################################################################
############################# Shared Library #######################################
####################################################################################
if (WITH_OBJECT AND NOT STATICBUILD) 
   add_library(
      "lib-${PROJECT_NAME}-shared"
      SHARED
      $<TARGET_OBJECTS:${PROJECT_NAME}-object>
   )

   target_link_libraries(
      "lib-${PROJECT_NAME}-shared"
      PUBLIC
      ${BLAS_LIBRARIES}
      ${LAPACK_LIBRARIES}
      $<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>
      $<$<BOOL:${WITH_TBLITE}>:tblite::tblite>
      $<$<BOOL:${WITH_GFN0}>:gfn0::gfn0>
      $<$<BOOL:${WITH_GFNFF}>:gfnff::gfnff>
      $<$<BOOL:${WITH_LIBPVOL}>:pvol::pvol>
      $<$<BOOL:${WITH_TOMLF}>:toml-f::toml-f>
      $<$<BOOL:${WITH_LWONIOM}>:lwoniom::lwoniom>
      $<$<BOOL:${WITH_FMLIP_RELAY}>:fmlip_relay::fmlip_relay> 
   )

   set_target_properties(
      "lib-${PROJECT_NAME}-shared"
      PROPERTIES
      Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
      LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
      OUTPUT_NAME "${PROJECT_NAME}"
      VERSION "${PROJECT_VERSION}"
      SOVERSION "${PROJECT_VERSION_MAJOR}"
   )

   target_include_directories(
      "lib-${PROJECT_NAME}-shared"
      PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
   )
endif()

###############################################################################
############################### Executables ###################################
###############################################################################
add_executable(
  ${PROJECT_NAME}-exe
  ${prog}
)

target_link_libraries(
   ${PROJECT_NAME}-exe
   PRIVATE
   "lib-${PROJECT_NAME}-static"
   $<$<BOOL:${STATICBUILD}>:-static>
)

set_target_properties(
  ${PROJECT_NAME}-exe
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  OUTPUT_NAME "${PROJECT_NAME}"
)

target_include_directories(
   ${PROJECT_NAME}-exe 
   PRIVATE 
   ${PROJECT_SOURCE_DIR}/include
   $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)


###############################################################################
########################### Install Instructions ##############################
###############################################################################

# Build output artifacts
if (WITH_OBJECT AND NOT STATICBUILD)
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
#      "lib-${PROJECT_NAME}-shared"
      "${PROJECT_NAME}-exe"
      EXPORT "${PROJECT_NAME}-targets"
      LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )
else()
   install(
      TARGETS
      "lib-${PROJECT_NAME}-static"
      "${PROJECT_NAME}-exe"
      EXPORT "${PROJECT_NAME}-targets"
      ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
      RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
      INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
   )
endif()
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

##############################################################################
############################# Unit Tests #####################################
##############################################################################

if (WITH_TESTS)
  enable_testing()
  add_subdirectory("test")
endif()
