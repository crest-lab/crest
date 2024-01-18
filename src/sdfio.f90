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

subroutine new_wrsdfens(env,fname,oname,conf)
  use iso_fortran_env,only:wp => real64
  use iomod
  use crest_data
  use strucrd
  use zdata,only:readwbo
  implicit none
  type(systemdata) :: env
  character(len=*),intent(in) :: fname
  character(len=*),intent(in) :: oname
  logical,intent(in),optional :: conf
  integer  :: nat,nall
  integer,allocatable  :: at(:)
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  real(wp),allocatable :: c0(:,:)
  real(wp),allocatable :: wbo(:,:)
  real(wp),allocatable :: icharges(:)
  integer :: i,j,ich
  real(wp) :: er
  logical :: ex,loopwbo,atmchrg
  character(len=120) :: sdfcomment
  loopwbo = .false.
  atmchrg = .false.
  if (present(conf)) loopwbo = conf
  !>--- read existing ensemble
  call rdensembleparam(fname,nat,nall)
  allocate (at(nat),eread(nat),xyz(3,nat,nall))
  call rdensemble(fname,nat,nall,at,xyz,eread)
  allocate (wbo(nat,nat),c0(3,nat),source=0.0_wp)
  !>--- determine how to obtain wbos
  wbo = 0.0_wp
  inquire (file='wbo',exist=ex)
  if (ex.and..not.loopwbo) then
    call readwbo('wbo',nat,wbo)
  elseif (.not.ex.and..not.loopwbo) then
    call xtbsp(env,0) !> gfn0 singlepoint
    call readwbo('wbo',nat,wbo)
  end if
  !>--- (optional) some special settings
  if (env%properties == p_protonate) atmchrg = .true.

  if (atmchrg) allocate (icharges(nat),source=0.0_wp)

  !>--- open sdf output file
  open (newunit=ich,file=oname)
  do i = 1,nall
    write (sdfcomment,'(a,i0,a,i0)') 'structure ',i,' of ',nall
    c0(1:3,1:nat) = xyz(1:3,1:nat,i)
    er = eread(i)
    if (loopwbo) then
      call wrxyz('tmpstruc.xyz',nat,at,c0)
      call xtbsp2('tmpstruc.xyz',env) !> singlepoint for wbos
      call readwbo('wbo',nat,wbo)
    end if
    if (atmchrg) then
      if (env%properties == p_protonate) then
        call set_prot_icharges(nat,wbo,icharges)
      end if
      call wrsdf(ich,nat,at,c0,er,env%chrg,wbo,sdfcomment,icharges)
    else
      call wrsdf(ich,nat,at,c0,er,env%chrg,wbo,sdfcomment)
    end if
  end do
  close (ich)

  if (allocated(icharges)) deallocate (icharges)
  deallocate (c0,wbo,xyz,eread,at)

contains
  subroutine set_prot_icharges(nat,wbo,icharges)
    !>--- special routine for protonation mode
    !     find the atom on which the proton was set
    !     and set its charge to 1. The added proton is
    !     always the last in the list, k=nat
    integer,intent(in) :: nat
    real(wp),intent(in) :: wbo(nat,nat)
    real(wp),intent(out) :: icharges(nat)
    integer :: i,j,k
    icharges = 0.0_wp
    k = nat
    do i = 1,nat
      if (nint(wbo(i,k)) .ne. 0) then
        icharges(i) = 1.0_wp
      end if
    end do
    return
  end subroutine set_prot_icharges
end subroutine new_wrsdfens
