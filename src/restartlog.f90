!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

!> global variables to keep track of restart

module crest_restartlog
  use crest_parameters
  use crest_data
  implicit none
  private

  logical,parameter :: debug = .false.

!>--- tracking variables
  integer :: restart_tracker = 0
  integer :: restart_goal = 0

  logical,allocatable :: last_processed(:)
  character(len=300) :: last_dumped

  !> a backup of the crest envrionment
  type(systemdata),allocatable :: restart_env

  !> backup of the last processed ensemble
  character(len=300) :: last_file
  integer :: last_nat = 0
  integer :: last_nall = 0
  integer,allocatable  :: last_at(:)
  real(wp),allocatable :: last_xyz(:,:,:)
  character(len=128),allocatable :: last_comments(:)

!>--- routines/functions
  public :: trackrestart
  public :: restart_save_env
  public :: trackensemble

  public :: dump_restart,read_restart

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!
!> Tracking routines to be called within the algos

  function trackrestart(env) result(skip)
!****************************************
!* This function is to be called both
!* to increment the restart tracker, and
!* to check if a step needs to be skiped
!****************************************
    implicit none
    logical :: skip
    type(systemdata),intent(in),optional :: env
    skip = .false.
    restart_tracker = restart_tracker+1
    if (restart_tracker < restart_goal) skip = .true.
    if (.not.skip.and.present(env)) then
      call restart_save_env(env)
    end if
    if (debug) write (stdout,*) 'RESTART_TRACKER =',restart_tracker
  end function trackrestart


  subroutine trackensemble(fname,nat,nall,at,xyz,comments)
!*******************************************************
!* This subroutine decides wether to track the ensemble
!*******************************************************
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in)  :: nat,nall
     integer,intent(in)  :: at(nat)
     real(wp),intent(in) :: xyz(3,nat,nall)
     character(len=*),intent(in) :: comments(nall)

    if (restart_tracker > restart_goal)then
      call restart_save_ensemble(fname,nat,nall,at,xyz,comments)
    if (debug) write (stdout,*) 'RESTART_ENSEMBLE = ',trim(fname)
    endif
  end subroutine trackensemble



!========================================================================================!
!========================================================================================!
!> DUMP to binary routines

  subroutine dump_restart()
    implicit none
    integer :: ich,i,j,k,l
    character(len=250) :: atmp
    if (debug) write (stdout,*) 'RESTART DEBUG dump summary'

    open (newunit=ich,file='crest.restart',status='replace',form='unformatted')

    write (ich) restart_tracker
    if (debug) write (stdout,*) restart_tracker

    if (allocated(restart_env)) then
      atmp = restart_env%cmd
      write (ich) atmp
      if (debug) write (stdout,*) 'cmd: ',trim(atmp)

      atmp = restart_env%inputcoords
      write (ich) atmp
      if (debug) write (stdout,*) 'inputcoords: ',trim(atmp)

      write (ich) restart_env%eprivious
      if (debug) write (stdout,*) 'eprivious: ',restart_env%eprivious

      write (ich) restart_env%elowest
      if (debug) write (stdout,*) 'elowest: ',restart_env%elowest

      j = restart_env%ref%nat
      write(ich) j
      if (debug) write (stdout,*) 'ref natoms: ', j
      do i=1,j
        write(ich) restart_env%ref%at(i)
      enddo
      do i=1,j
        write(ich) restart_env%ref%xyz(1:3,i)
      enddo
   
    end if

    call dump_last_ensemble(ich)

    close (ich)
  end subroutine dump_restart

  subroutine dump_last_ensemble(ich)
!******************************************
!* dump last saved ensemble as binary data
!******************************************
    implicit none
    integer, intent(in) :: ich
    integer :: nat,nall,i,j,k,l
    write(ich) last_file
    nat = last_nat
    write(ich) nat
    nall = last_nall
    write(ich) nall
    if(allocated(last_comments) .and. allocated(last_xyz))then
    do k=1,nat
      write(ich) last_at(k)
    enddo
    do i=1,nall
      write(ich) last_comments(i)
      do j=1,nat
        write(ich) last_xyz(1:3,j,i)
      enddo
    enddo
    endif
  end subroutine dump_last_ensemble

!========================================================================================!
!========================================================================================!
!> read from binary subroutines

  subroutine read_restart(env)
    implicit none
    type(systemdata),intent(inout) :: env
    integer :: ich,i,j,k,l
    character(len=250) :: atmp
    real(wp) :: rdum,xyzdum(3)
    integer :: idum
    open (newunit=ich,file='crest.restart',status='old',form='unformatted')

    read (ich) restart_goal

    read (ich) env%cmd
    read (ich) env%inputcoords

    read (ich) env%eprivious
    read (ich) env%elowest

    read (ich) j
    do i=1,j
      read(ich) idum
    enddo
    do i=1,j
      read(ich) xyzdum(1:3)
    enddo
   
    close (ich)
  end subroutine read_restart

  subroutine read_last_ensemble(ich)
!******************************************
!* dump last saved ensemble as binary data
!******************************************
    implicit none
    integer, intent(in) :: ich
    integer :: nat,nall,i,j,k,l
    read(ich) last_file
    read(ich) nat
    nat = last_nat
    read(ich) nall
    nall = last_nall
    if(nat > 0 .and. nall > 0) then
      allocate(last_at(nat))
      allocate(last_xyz(3,nat,nall))
      allocate(last_comments(nall)) 
      do k=1,nat
        read(ich) last_at(k)
      enddo
      do i=1,nall
        read(ich) last_comments(i)
        do j=1,nat
          read(ich) last_xyz(1:3,j,i)
        enddo
      enddo
    endif
  end subroutine read_last_ensemble



!========================================================================================!
!========================================================================================!
!> some routines to create backup data

  subroutine restart_save_env(env)
!*************
!* backup env
!*************  
    implicit none
    type(systemdata),intent(in) :: env
    if (.not.allocated(restart_env)) then
      allocate (restart_env,source=env)
    end if
    restart_env = env
  end subroutine restart_save_env


  subroutine restart_save_ensemble(fname,nat,nall,at,xyz,comments)
!*********************************
!* backup last processed ensemble
!*********************************
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in)  :: nat,nall
     integer,intent(in)  :: at(nat)
     real(wp),intent(in) :: xyz(3,nat,nall)
     character(len=*),intent(in) :: comments(nall)
     integer :: i
     !> backup of the last processed ensemble
     last_file = trim(fname)
     last_nat = nat
     last_nall = nall
     if(.not.allocated(last_at)) allocate(last_at(nat))
     last_at(:) = at(:)
     if(allocated(last_xyz)) deallocate(last_xyz)
     if(allocated(last_comments)) deallocate(last_comments)
     allocate(last_xyz(3,nat,nall))
     allocate(last_comments(nall))
     last_xyz(:,:,:) = xyz(:,:,:)
     last_comments(:) = comments(:)
  end subroutine restart_save_ensemble

!========================================================================================!
!========================================================================================!
end module crest_restartlog
