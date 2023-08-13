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
  "${dir}/singlepoint.f90"
  "${dir}/dynamics.f90"
  "${dir}/numhess.f90"
  "${dir}/optimization.f90"
  "${dir}/playground.f90"
  "${dir}/crossing.f90"
  "${dir}/refine.f90"
  "${dir}/scan.f90" 
  "${dir}/search_1.f90"
  "${dir}/search_mecp.f90"
  "${dir}/setuptest.f90"
  "${dir}/hessian_tools.f90"
  "${dir}/search_conformers.f90"
  "${dir}/search_entropy.f90"
  "${dir}/parallel.f90"
)

set(srcs ${srcs} PARENT_SCOPE)






