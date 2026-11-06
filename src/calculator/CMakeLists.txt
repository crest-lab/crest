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
  "${dir}/calculator.F90"
  "${dir}/external_engrad.f90"
  "${dir}/calc_type.f90"
  "${dir}/printouts.F90"
  "${dir}/constraints.f90"
  "${dir}/gfn0_api.F90"
  "${dir}/gfnff_api.F90"
  "${dir}/lj.f90"
  "${dir}/lwoniom.F90"
  "${dir}/oniom_hessian.F90"
  "${dir}/nonadiabatic.f90"
  "${dir}/tblite_api.F90"
  "${dir}/components/rmsdpot.f90"
  "${dir}/api_helpers.F90"
  "${dir}/api_engrad.f90"
  "${dir}/gradreader.f90"
  "${dir}/components/libpvol.F90"
  "${dir}/subprocess_types.f90"
  "${dir}/xtb_sc.f90"
  "${dir}/orca_sc.f90"
  "${dir}/generic_sc.f90"
  "${dir}/turbom_sc.f90"
  "${dir}/subprocess_engrad.f90"
  "${dir}/hessian_reconstruct.f90"
  "${dir}/hr_utils.f90"
  "${dir}/modelhessians.f90"
  "${dir}/approxg.f90"
  "${dir}/penalty.f90"
  "${dir}/mlip_sc.F90"
  "${dir}/components/ddx_pc.F90"
  "${dir}/components/electrostatic.F90"
  "${dir}/components/surface.F90"
  "${dir}/components/solvation.F90"
)

set(srcs ${srcs} PARENT_SCOPE)








