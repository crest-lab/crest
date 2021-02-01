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
  version: '2.11.0',
  meson_version: '>=0.53',
  default_options: [
    'buildtype=debugoptimized',
  ],
)

# General configuration information
exe_deps = []
subdir('config')

# Collect source of the project
prog = []
srcs = []
subdir('src')

crest_exe = executable(
  meson.project_name(),
  sources: srcs + prog,
  dependencies: exe_deps,
  install: true,
  link_language: 'fortran',
)

install_data(
  files('combi-scripts/crest_combi'),
  install_dir: get_option('bindir'),
  install_mode: 'rwxr--r--',
)
