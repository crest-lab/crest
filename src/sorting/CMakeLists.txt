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
  "${dir}/canonical.f90"
  "${dir}/ccegen_interface.f90"
  "${dir}/ccegen_utils.f90"
  "${dir}/ccegen.f90"
  "${dir}/cregen_interfaces.f90"
  "${dir}/cregen.f90"
  "${dir}/cregen_utils.f90"
  "${dir}/ensemblecomp.f90"
  "${dir}/hungarian.f90"
  "${dir}/irmsd_module.f90"
  "${dir}/ls_rmsd.f90"
  "${dir}/pbc_fingerprint.f90"
  "${dir}/quicksort.f90"
  "${dir}/rotcompare.f90"
"${dir}/unionize.f90"
  "${dir}/zdata.f90"
  "${dir}/ztopology.f90"
)

set(srcs ${srcs} PARENT_SCOPE)






