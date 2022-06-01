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

add_subdirectory("qcg")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND srcs
  "${dir}/acidbase.f90"
  "${dir}/atmasses.f90"
  "${dir}/axis_module.f90"
  "${dir}/biasmerge.f90"
  "${dir}/bondconstraint.f90"
  "${dir}/ccegen.f90"
  "${dir}/classes.f90"
  "${dir}/cleanup.f90"
  "${dir}/compress.f90"
  "${dir}/confcross.f90"
  "${dir}/confopt.f90"
  "${dir}/confparse.f90"
  "${dir}/confscript1.f90"
  "${dir}/confscript2_main.f90"
  "${dir}/confscript2_misc.f90"
  "${dir}/confscript3.f90"
  "${dir}/constraining.f90"
  "${dir}/cregen.f90"
  "${dir}/data.f90"
  "${dir}/deprotonate.f90"
  "${dir}/dft_propcalc.f90"
  "${dir}/ensemblecomp.f90"
  "${dir}/entropic.f90"
  "${dir}/entropy.f90"
  "${dir}/eval_timer.f90"
  "${dir}/filemod.f90"
  "${dir}/flexi.f90"
  "${dir}/floydwarshall.f90"
  "${dir}/freqmasses.f90"
  "${dir}/geo.f90"
  "${dir}/identifiers.f90"
  "${dir}/intmodes.f"
  "${dir}/io.f"
  "${dir}/iomod.F90"
  "${dir}/ligands.f90"
  "${dir}/ls_rmsd.f90"
  "${dir}/marqfit.f"
  "${dir}/minitools.f90"
  "${dir}/msmod.f90"
  "${dir}/msreact.f90"
  "${dir}/ncigeo.f90"
  "${dir}/neighbor.f"
  "${dir}/newcregen.f90"
  "${dir}/normmd.f90"
  "${dir}/ompmklset.F90"
  "${dir}/optinplace.f90"
  "${dir}/pka.f90"
  "${dir}/printouts.f90"
  "${dir}/prmat.f"
  "${dir}/progressbar.f90"
  "${dir}/propcalc.f90"
  "${dir}/protonate.f90"
  "${dir}/quicksort.f90"
  "${dir}/reactor.f90"
  "${dir}/readl2.f"
  "${dir}/readl.f"
  "${dir}/relaxensemble.f90"
  "${dir}/rotcompare.f90"
  "${dir}/scratch.f90"
  "${dir}/sdfio.f90"
  "${dir}/select.f90"
  "${dir}/setdg.f90"
  "${dir}/shortpath.f90"
  "${dir}/signal.c"
  "${dir}/sigterm.f90"
  "${dir}/sortens.f90"
  "${dir}/spline.f90"
  "${dir}/stereoisomers.f90"
  "${dir}/strucreader.f90"
  "${dir}/symmetry2.f90"
  "${dir}/symmetry_i.c"
  "${dir}/syscheck.f90"
  "${dir}/tautomerize.f90"
  "${dir}/thermocalc.f90"
  "${dir}/thermo.f90"
  "${dir}/trackorigin.f90"
  "${dir}/trialmd.f90"
  "${dir}/utilities.f"
  "${dir}/xtbmodef.f90"
  "${dir}/zdata.f90"
  "${dir}/zsort.f90"
  "${dir}/ztopology.f90"
)

list(APPEND prog
  "${dir}/crest_main.f90"
)

set(srcs ${srcs} PARENT_SCOPE)
set(prog ${prog} PARENT_SCOPE)
