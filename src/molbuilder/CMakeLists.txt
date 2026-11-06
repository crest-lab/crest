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

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/rigidconf.f90"
  "${dir}/analyze.f90"
  "${dir}/reconstruct.f90"
  "${dir}/construct.f90"
  "${dir}/construct_list.f90"
  "${dir}/classify.f90"
  "${dir}/classify_type.f90"
  "${dir}/classify_func.f90"
  "${dir}/ttconf_settings.f90"
  "${dir}/ttconf_ringsample.f90"
  "${dir}/ttconf_rings.f90"
  "${dir}/ttconf_light_mod.f90"
)

set(srcs ${srcs} PARENT_SCOPE)






