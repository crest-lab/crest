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
  "${dir}/acidbase.f90"
  "${dir}/confscript2_main.f90"
  "${dir}/confscript2_misc.f90"
  "${dir}/confscript3.f90"
  "${dir}/confopt.f90"
  "${dir}/cregen_old.f90"
  "${dir}/deprotonate.f90"
  "${dir}/flexi.f90"
  "${dir}/ligands.f90"
  "${dir}/normmd.f90"
  "${dir}/pka.f90"
  "${dir}/protonate.f90"
  "${dir}/reactor.f90"
  "${dir}/relaxensemble.f90"
  "${dir}/stereoisomers.f90"
  "${dir}/tautomerize.f90"
  "${dir}/trialmd_legacy.f90"
  "${dir}/zsort.f90"
)

set(srcs ${srcs} PARENT_SCOPE)






