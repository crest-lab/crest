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

project(
  'crest',
  'fortran', 'c',
  version: '3.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.63',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=static',
    'c_link_args=-static',
    'fortran_link_args=-static',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')


# =================================== #
## General configuration information ##
# =================================== #
exe_deps = []
subdir('config')

# create the metadata file with the configured data
configure_file(
  input: files('assets/template/metadata.f90'),
  output: 'crest_metadata.fh',
  configuration : config,
)


# Documentation
#subdir('docs')

# Collect source of the project
prog = []
srcs = []
subdir('src')


# Create library target
crest_lib = library(
  meson.project_name(),
  sources: srcs,
  dependencies: exe_deps,
#  include_directories: include_directories('include'),
)

# Export as dependency
crest_inc = crest_lib.private_dir_include()
crest_dep = declare_dependency(
  link_with: crest_lib,
  include_directories: crest_inc,
  dependencies: exe_deps,
)


# Create executable target
crest_exe = executable(
  meson.project_name(),
  sources: prog,
  dependencies: crest_dep,
  install: install,
  link_language: 'fortran',
)



# add the testsuite separate meson.build 
#subdir('testsuite') #has to be filled with tests, not availabel yet
