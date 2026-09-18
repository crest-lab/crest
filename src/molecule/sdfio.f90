!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2025 Philipp Pracht
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

!---------------------------------------------------------------------------------
! Routines for the handling of molecules in the chemoinformatical *.SDF format
!---------------------------------------------------------------------------------
subroutine inpsdf(env,fname)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=*) :: fname
  type(coord) :: mol
  integer :: i
  env%sdfformat = .true.
  call checkcoordtype(fname,i)
  if (any((/31,32/) == i)) then
    call mol%open(fname)
  else
    error stop 'file not in sdf format'
  end if
  return
end subroutine inpsdf

!================================================================================!

subroutine new_wrsdfens(env,fname,oname,conf)
  !***********************************************************************
  !* Write a conformer ensemble as an SDF file.
  !* Bond orders (WBO) are obtained via a singlepoint calculation using
  !* the calculator configured in env.
  !*
  !* Input:
  !*  env  -  crest's systemdata object (provides calculator and charge)
  !*  fname - input XYZ ensemble file
  !*  oname - output SDF file name
  !*  conf  - if .true., run a separate SP for each structure (loopwbo)
  !***********************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  implicit none
  type(systemdata),intent(inout) :: env
  character(len=*),intent(in) :: fname
  character(len=*),intent(in) :: oname
  logical,intent(in) :: conf
  !> local
  type(coord),allocatable :: structures(:)
  type(calcdata) :: tmpcalc
  real(wp),allocatable :: wbo(:,:)
  real(wp),allocatable :: grad(:,:)
  real(wp),allocatable :: icharges(:)
  integer :: nall,nat,ich,io,i
  real(wp) :: energy,er
  logical :: loopwbo,atmchrg
  character(len=120) :: sdfcomment

  atmchrg = .false.
  loopwbo = conf

  ! ── read ensemble as array of coord objects (xyz in Bohr) ─────────────
  call rdensemble(fname,nall,structures)
  nat = structures(1)%nat

  ! ── set up a minimal GFN0 singlepoint calculator for WBOs ─────────────
  call tmpcalc%create('gfn0',chrg=env%chrg,uhf=env%uhf)
  tmpcalc%calcs(1)%rdwbo = .true.
  allocate (wbo(nat,nat),grad(3,nat),source=0.0_wp)
  energy = 0.0_wp

  ! ── for non-loopwbo: one SP on the first structure, WBO reused for all ─
  if (.not.loopwbo) then
    call engrad(structures(1),tmpcalc,energy,grad,io)
    if (allocated(tmpcalc%calcs(1)%wbo)) wbo = tmpcalc%calcs(1)%wbo
  end if

  ! ── (optional) special per-atom charge handling ───────────────────────
  if (env%properties == p_protonate) atmchrg = .true.
  if (atmchrg) allocate (icharges(nat),source=0.0_wp)

  ! ── write SDF output ──────────────────────────────────────────────────
  open (newunit=ich,file=oname)
  do i = 1,nall
    write (sdfcomment,'(a,i0,a,i0)') 'structure ',i,' of ',nall
    er = structures(i)%energy
    if (loopwbo) then
      ! ── per-structure SP for bond orders (protonation/tautomer modes) ─
      wbo = 0.0_wp
      call engrad(structures(i),tmpcalc,energy,grad,io)
      if (allocated(tmpcalc%calcs(1)%wbo)) wbo = tmpcalc%calcs(1)%wbo
    end if
    if (atmchrg) then
      if (env%properties == p_protonate) then
        call set_prot_icharges(nat,wbo,icharges)
      end if
      !> wrsdf expects Angstrom: multiply Bohr coordinates by bohr (Å/bohr)
      call wrsdf(ich,nat,structures(i)%at,structures(i)%xyz*bohr, &
      &          er,env%chrg,wbo,sdfcomment,icharges)
    else
      call wrsdf(ich,nat,structures(i)%at,structures(i)%xyz*bohr, &
      &          er,env%chrg,wbo,sdfcomment)
    end if
  end do
  close (ich)

  call tmpcalc%reset()
  if (allocated(icharges)) deallocate (icharges)
  deallocate (wbo,grad,structures)

contains
  subroutine set_prot_icharges(nat,wbo,icharges)
    !***********************************************
    !* For protonation mode: locate the heavy atom
    !* bonded to the added proton (last in list)
    !* and assign it a formal charge of +1.
    !***********************************************
    integer,intent(in) :: nat
    real(wp),intent(in) :: wbo(nat,nat)
    real(wp),intent(out) :: icharges(nat)
    integer :: i,k
    icharges = 0.0_wp
    k = nat
    do i = 1,nat
      if (nint(wbo(i,k)) .ne. 0) then
        icharges(i) = 1.0_wp
      end if
    end do
  end subroutine set_prot_icharges
end subroutine new_wrsdfens

!================================================================================!

subroutine crest_ensemble_reformat(env)
  !***************************************************************
  !* Reformat the conformer ensemble into requested alternative
  !* file formats (currently SDF) after the run completes.
  !*
  !* Input:
  !*  env  -  crest's systemdata object
  !***************************************************************
  use crest_parameters
  use crest_data
  implicit none
  type(systemdata),intent(inout) :: env

  ! ── SDF ensemble output ───────────────────────────────────────
  if (env%outputsdf .or. env%sdfformat) then
    if (any((/crest_mfmdgc,crest_imtd,crest_imtd2/) == env%crestver)) then
      call new_wrsdfens(env,conformerfile,conformerfilebase//'.sdf',.false.)
    end if
    if (any((/crest_screen,crest_mdopt/) == env%crestver)) then
      call new_wrsdfens(env,'crest_ensemble.xyz','crest_ensemble.sdf',.false.)
    end if
  end if

end subroutine crest_ensemble_reformat
