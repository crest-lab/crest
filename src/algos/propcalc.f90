!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

!> Skeleton stub for property calculations on an ensemble.
!> Each mode needs to be implemented via modern algo routines.
!> The legacy implementation (system-call/I/O based) has been removed.
!> See git history for reference.

subroutine propcalc(iname,imode,env,tim)
  use crest_parameters
  use crest_data
  use cregen_interface
  implicit none
  character(len=*),intent(in) :: iname
  integer,intent(in) :: imode
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout) :: tim

  select case (imode)
  case (p_prop_hess)
    !> TODO: Hessian calculations for all conformers (was: xtb --hess)
  case (p_prop_autoir)
    !> TODO: IR spectrum averaging over populated conformers (was: autoir + xtb --ohess)
  case (p_prop_ohess)
    !> TODO: Optimization + Hessian for all conformers (was: xtb --ohess)
  case (p_prop_gsolv)
    !> TODO: Free energy in solvation, 2-step (was: xtb --sp + xtb --ohess)
  case (p_prop_reopt)
    !> TODO: Vtight reoptimization for all conformers (was: xtb --opt vtight)
  case (p_prop_finalhess)
    call crest_finalhess(iname,env,tim)

  case (p_prop_multilevel:p_prop_multilevel+9)
    !> Post-search processing of the conformer ensemble at a higher level.
    !> Dispatched by job number; input is typically crest_rotamers.xyz.
    block
      integer :: saved_stage
      saved_stage = env%calc%refine_stage
      select case (imode)
      case (p_prop_multilevel+1)  !> A@B post-search geo-opt (existing)
        env%calc%refine_stage = refine%post_opt
        call crest_multilevel_reopt(iname,env,tim)
      case (p_prop_multilevel+2)  !> --rerank post-search SP re-ranking
        call crest_rerank_sp(iname,env,tim)
      case (p_prop_multilevel+3)  !> --reopt post-search geo-opt (standalone)
        env%calc%refine_stage = refine%post_reopt
        call crest_multilevel_reopt(iname,env,tim)
      case default
        write (stdout,'(a,i0,a)') 'propcalc: multilevel mode ',imode,' not implemented'
      end select
      env%calc%refine_stage = saved_stage
    end block
  case (p_prop_dipole)
    !> TODO: Singlepoint + dipole extraction (was: xtb --sp, grep molecular dipole)
  case (p_prop_rerank)
    !> TODO: Singlepoint + reranking (was: xtb --sp + newcregen)
  case default
    write (stdout,'(a,i0,a)') 'propcalc: mode ',imode,' not yet implemented'
  end select

end subroutine propcalc

!========================================================================================!

subroutine crest_multilevel_reopt(iname,env,tim)
!*******************************************************************
!* Read the ensemble iname, optimize all structures using the
!* calculator whose refine_lvl matches env%calc%refine_stage,
!* sort via CREGEN, and write crest_reopt.xyz.
!*
!* The caller is responsible for setting env%calc%refine_stage
!* to the desired level before this routine is called.
!*
!* Input:
!*   iname  - path to input ensemble (e.g. crest_rotamers.xyz)
!* Output:
!*   crest_reopt.xyz (sorted conformer ensemble at the higher level)
!*******************************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use parallel_interface
  use cregen_interface
  use iomod,only:drawbox,catdel
  implicit none
  character(len=*),intent(in) :: iname
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout) :: tim
  integer :: nat,nall,T,Tn
  type(coord),allocatable :: structures(:)
  character(len=*),parameter :: outname = 'crest_reopt.xyz'
  logical :: ex

  inquire(file=iname,exist=ex)
  if (.not.ex) then
    write(stdout,'(a,a,a)') '**WARNING** ',trim(iname),' not found, skipping multilevel reopt'
    return
  end if

  call tim%start(16,'Multilevel reopt')

  call rdensembleparam(iname,nat,nall)
  if (nall < 1) then
    write(stdout,*) '**WARNING** empty ensemble, skipping multilevel reopt'
    call tim%stop(16)
    return
  end if
  call rdensemble(iname,nall,structures)

  call new_ompautoset(env,'auto',nall,T,Tn)

  write(stdout,*)
  call drawbox(stdout,'MULTILEVEL ENSEMBLE REOPT',charset=7,width=51,ltab=10)
  write(stdout,'(1x,a,i0,a,1x,a)') &
    & 'Re-optimizing ',nall,' structures of file ',trim(iname)

! ── refine_stage is set by the caller; run geo-opt ───────────────
  call crest_oloop(env,nall,structures,.true.)

! ── write output ─────────────────────────────────────────────────
  call wrensemble(outname,nall,structures)

  write(stdout,'(/,a,a,a)') 'Re-optimized ensemble written to <',outname,'>'

! ── sort via CREGEN ──────────────────────────────────────────────
  call newcregen(env,0,outname)
  call catdel('cregen.out.tmp')

  deallocate(structures)
  call tim%stop(16)
end subroutine crest_multilevel_reopt

!========================================================================================!

subroutine crest_rerank_sp(iname,env,tim)
!*******************************************************************
!* Read the ensemble iname, run single-point energies using the
!* calculator tagged with refine_lvl == refine%post_sp (= 11,
!* set by the --rerank keyword), re-sort via CREGEN, and write
!* crest_reopt.xyz. Geometries are not changed.
!*
!* Input:
!*   iname  - path to input ensemble (e.g. crest_rotamers.xyz)
!* Output:
!*   crest_reopt.xyz (ensemble re-ranked by higher-level SP energies)
!*******************************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use parallel_interface
  use cregen_interface
  use iomod,only:drawbox,catdel
  implicit none
  character(len=*),intent(in) :: iname
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout) :: tim
  integer :: nat,nall,T,Tn,old_stage
  type(coord),allocatable :: structures(:)
  character(len=*),parameter :: outname = 'crest_reopt.xyz'
  logical :: ex

  inquire(file=iname,exist=ex)
  if (.not.ex) then
    write(stdout,'(a,a,a)') '**WARNING** ',trim(iname),' not found, skipping SP rerank'
    return
  end if

  call tim%start(16,'Post-search SP rerank')

  call rdensembleparam(iname,nat,nall)
  if (nall < 1) then
    write(stdout,*) '**WARNING** empty ensemble, skipping SP rerank'
    call tim%stop(16)
    return
  end if
  call rdensemble(iname,nall,structures)

  call new_ompautoset(env,'auto',nall,T,Tn)

  write(stdout,*)
  call drawbox(stdout,'POST-SEARCH SP RE-RANKING',charset=7,width=51,ltab=10)
  write(stdout,'(1x,a,i0,a,1x,a)') &
    & 'Re-ranking ',nall,' structures of file ',trim(iname)

! ── activate only the SP reranking calculator ────────────────────
  old_stage = env%calc%refine_stage
  env%calc%refine_stage = refine%post_sp

  call crest_sploop(env,nall,structures)

  env%calc%refine_stage = old_stage

! ── write output ─────────────────────────────────────────────────
  call wrensemble(outname,nall,structures)

  write(stdout,'(/,a,a,a)') 'Re-ranked ensemble written to <',outname,'>'

! ── sort via CREGEN ──────────────────────────────────────────────
  call newcregen(env,0,outname)
  call catdel('cregen.out.tmp')

  if (allocated(structures)) deallocate(structures)
  call tim%stop(16)
end subroutine crest_rerank_sp

!========================================================================================!

subroutine crest_finalhess(iname,env,tim)
!*******************************************************************
!* Run Hessians + thermochemistry on the final conformer ensemble,
!* then re-sort via CREGEN using Gibbs free energies.
!* Falls back from crest_conformers.xyz to crest_ensemble.xyz.
!* Uses the main calculator as-is (no separate refine level).
!* Input:
!*   iname  - primary input file (crest_conformers.xyz)
!* Output:
!*   crest_hess.xyz (intermediate), then CREGEN overwrites
!*   crest_conformers.xyz / crest_ensemble.xyz sorted by Gfree
!*******************************************************************
  use crest_parameters,only:wp,stdout,bohr
  use crest_data
  use crest_calculator
  use strucrd
  use parallel_interface
  use cregen_interface
  use iomod,only:drawbox,catdel
  implicit none
  character(len=*),intent(in) :: iname
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout) :: tim
  character(len=:),allocatable :: infile
  integer :: nat,nall,T,Tn
  real(wp),allocatable :: xyz(:,:,:),eread(:),etmp(:)
  integer,allocatable  :: at(:)
  character(len=*),parameter :: outname = 'crest_hess.xyz'
  logical :: ex

! ── select input file with fallback to crest_ensemble.xyz ────────
  inquire(file=iname,exist=ex)
  if (ex) then
    infile = iname
  else
    inquire(file=ensemblefile,exist=ex)
    if (ex) then
      infile = ensemblefile
    else
      write(stdout,'(a)') '**WARNING** no conformer ensemble found, skipping --finalhess'
      return
    end if
  end if

  call tim%start(16,'Final ensemble Hessians')

  call rdensembleparam(infile,nat,nall)
  if (nall < 1) then
    write(stdout,*) '**WARNING** empty ensemble, skipping --finalhess'
    call tim%stop(16)
    return
  end if
  allocate(xyz(3,nat,nall),at(nat),eread(nall),etmp(nall))
  call rdensemble(infile,nat,nall,at,xyz,eread)
! ── crest_hessloop requires coordinates in Bohr ──────────────────
  xyz = xyz/bohr

  call new_ompautoset(env,'auto',nall,T,Tn)

  write(stdout,*)
  call drawbox(stdout,'FINAL ENSEMBLE HESSIANS',charset=7,width=51,ltab=10)
  write(stdout,'(1x,a,i0,a,1x,a)') &
    & 'Computing Hessians for ',nall,' structures of file ',trim(infile)

  call crest_hessloop(env,nat,nall,at,xyz,etmp)
  eread(:) = eread(:) + etmp(:)

! ── back to Angstrom, write intermediate file ─────────────────────
  xyz = xyz*bohr
  call wrensemble(outname,nat,nall,at,xyz,eread)
  write(stdout,'(/,a,a,a)') 'Hessian ensemble written to <',outname,'>'

! ── sort via CREGEN using Gibbs free energies ────────────────────
  call newcregen(env,0,outname)
  call catdel('cregen.out.tmp')

  deallocate(xyz,at,eread)
  call tim%stop(16)
end subroutine crest_finalhess
