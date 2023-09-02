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
  use iso_fortran_env,wp => real64,stdout => output_unit
  use crest_data
  implicit none
  private

  logical,parameter :: debug = .true.

!>--- tracking variables
  integer :: restart_tracker = 0
  integer :: restart_goal = 0

  character(len=300) :: restart_inputcoords
  character(len=300) :: restart_cmd

  character(len=300) :: last_file
  integer :: last_nall
  logical,allocatable :: last_processed(:)
  character(len=300) :: last_dumped

  !> a backup of the crest envrionment
  type(systemdata),allocatable :: restart_env


!>--- routines/functions
  public :: trackrestart
  public :: restart_save_env

  public :: dump_restart,read_restart

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

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
    if(.not.skip .and. present(env))then
      call restart_save_env(env)
    endif
    if(debug) write(stdout,*) 'RESTART_TRACKER =',restart_tracker
  end function trackrestart

!========================================================================================!

  subroutine dump_restart()
    implicit none
    integer :: ich
    character(len=250) :: atmp
    if(debug) write(stdout,*) 'RESTART DEBUG dump summary'

    open (newunit=ich,file='crest.restart',status='replace',form='unformatted')

    write(ich) restart_tracker
    if(debug) write(stdout,*) restart_tracker    


    if(allocated(restart_env))then
    atmp = restart_env%cmd
    write(ich) atmp
    if(debug) write(stdout,*) 'cmd: ',trim(atmp)


    atmp = restart_env%inputcoords
    write(ich) atmp
    if(debug) write(stdout,*) 'inputcoords: ',trim(atmp)

    write(ich) restart_env%eprivious
    if(debug) write(stdout,*) 'eprivious: ',restart_env%eprivious
    write(ich) restart_env%elowest
    if(debug) write(stdout,*) 'elowest: ',restart_env%elowest

    endif

    close (ich)
  end subroutine dump_restart


!========================================================================================!

  subroutine read_restart()
    implicit none
    integer :: ich
    character(len=250) :: atmp
    real(wp) :: rdum
    open (newunit=ich,file='crest.restart',status='old',form='unformatted')

    read(ich) restart_goal

    read(ich) atmp
    read(ich) atmp

    read(ich) rdum
    read(ich) rdum

    close (ich)
  end subroutine read_restart

!========================================================================================!
!> some routines to create backup data

  subroutine restart_save_input(inp,cmd)
    implicit none
    character(len=*),intent(in) :: inp
    character(len=*),intent(in) :: cmd
    restart_inputcoords = trim(inp)
    restart_cmd = trim(cmd)
  end subroutine restart_save_input

  subroutine restart_save_env(env)
     implicit none
     type(systemdata),intent(in) :: env
     if(.not.allocated(restart_env))then 
     allocate(restart_env, source = env)
     endif
     restart_env = env
  end subroutine restart_save_env 

!========================================================================================!
!========================================================================================!
end module crest_restartlog
