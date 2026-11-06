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
  "${dir}/aiss_param.F90"
  "${dir}/iff_lmo.f90"
  "${dir}/iff_ini.F90"
  "${dir}/iff_energy.f90"
  "${dir}/iff_prepare.F90"
  "${dir}/aiss_search.F90"
  "${dir}/aiss_dock.F90"
  "${dir}/aiss_standalone.F90"
)

set(srcs ${srcs} PARENT_SCOPE)
