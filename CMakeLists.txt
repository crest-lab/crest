# SPDX-License-Identifier: LGPL-3.0-or-later
# Project configuration
cmake_minimum_required(VERSION 3.14)
project(crest VERSION 3.0 LANGUAGES Fortran C)
set(SOVERSION "pre")

# Include GNU installation directories
include(GNUInstallDirs)

# General configuration and addition of subdirectories for further configurations and source
set(libs)  # external dependencies will be added
set(omp)   # omp dependency
add_subdirectory(config) # For detailed compiler options, feature toggles, and dependencies
set(srcs)  # sources will be added
set(prog)  # src/crest_main.f90 will be added
add_subdirectory(src)    # Source files



# Define executable with the main program source included
add_executable(${PROJECT_NAME} ${prog} ${srcs})
set_target_properties(${PROJECT_NAME} PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  OUTPUT_NAME ${PROJECT_NAME}
)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE ${libs})

# Installation rules for the executable
install(TARGETS ${PROJECT_NAME}
  EXPORT "${PROJECT_NAME}-targets"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

