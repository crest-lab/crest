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


add_subdirectory("parsing")
add_subdirectory("optimize")
add_subdirectory("calculator")
add_subdirectory("algos")
add_subdirectory("dynamics")
add_subdirectory("qcg")
add_subdirectory("qmhelpers")
add_subdirectory("graphs")
add_subdirectory("rigidconf")
add_subdirectory("discretize")
add_subdirectory("entropy")
add_subdirectory("legacy_algos")
add_subdirectory("msreact")


set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/atmasses.f90"
  "${dir}/axis_module.f90"
  "${dir}/biasmerge.f90"
  "${dir}/bondconstraint.f90"
  "${dir}/ccegen.f90"
  "${dir}/choose_settings.f90"
  "${dir}/classes.f90"
  "${dir}/cleanup.f90"
  "${dir}/compress.f90"
  "${dir}/confparse.f90"
  "${dir}/cregen.f90"
  "${dir}/crest_pars.f90"
  "${dir}/ensemblecomp.f90"
  "${dir}/eval_timer.f90"
  "${dir}/filemod.f90"
  "${dir}/flexi.F90"
  "${dir}/freqmasses.f90"
  "${dir}/geo.f90"
  "${dir}/identifiers.f90"
  "${dir}/internals.f90"
  "${dir}/internals2.f90"
  "${dir}/iomod.F90"
  "${dir}/legacy_wrappers.f90"
  "${dir}/ls_rmsd.f90"
  "${dir}/marqfit.f90"
  "${dir}/minitools.f90"
  "${dir}/miscdata.f90"
  "${dir}/ncigeo.f90"
  "${dir}/ompmklset.F90"
  "${dir}/printouts.f90"
  "${dir}/prmat.f90"
  "${dir}/propcalc.f90"
  "${dir}/quicksort.f90"
  "${dir}/readl.f90"
  "${dir}/restartlog.f90"
  "${dir}/rotcompare.f90"
  "${dir}/scratch.f90"
  "${dir}/sdfio.f90"
  "${dir}/select.f90"
  "${dir}/signal.c"
  "${dir}/sigterm.f90"
  "${dir}/sortens.f90"
  "${dir}/strucreader.f90"
  "${dir}/symmetry2.f90"
  "${dir}/symmetry_i.c"
  "${dir}/testmol.f90"
  "${dir}/timer.f90"
  "${dir}/trackorigin.f90"
  "${dir}/utilmod.f90"
  "${dir}/wallsetup.f90"
  "${dir}/zdata.f90"
  "${dir}/ztopology.f90"
)

list(APPEND prog
  "${dir}/crest_main.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
set(prog ${prog} PARENT_SCOPE)
