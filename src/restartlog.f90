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
  use miscdata, only: PSE
  use iomod, only: command
  implicit none
  private

  !logical,parameter :: debug = .true.
  logical,parameter :: debug = .false.
  logical,parameter :: saveensembles = .true.

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
  interface trackensemble
    module procedure :: trackensemble_comments
    module procedure :: trackensemble_energy
  end interface trackensemble
  public :: restart_write_dummy

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
    if (debug) write (stdout,*) '%%% RESTART_TRACKER =',restart_tracker

    if (restart_tracker < restart_goal) skip = .true.
    if (.not.skip.and.present(env)) then
      call restart_save_env(env)
      call dump_restart()
    end if
    if(restart_tracker == restart_goal-1)then
      if (debug) write (stdout,*) '%%% RESTART_RESTORE %%%'
      call restore_ensemble()
    endif
  end function trackrestart


  subroutine trackensemble_comments(fname,nat,nall,at,xyz,comments)
!*******************************************************
!* This subroutine decides wether to track the ensemble
!* Typically, this routine is called in CREGEN
!*******************************************************
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in)  :: nat,nall
     integer,intent(in)  :: at(nat)
     real(wp),intent(in) :: xyz(3,nat,nall)
     character(len=*),intent(in) :: comments(nall)

    if (restart_tracker > restart_goal)then
      call restart_save_ensemble(fname,nat,nall,at,xyz,comments)
    if (debug) write (stdout,*) '%%% RESTART_ENSEMBLE = ',trim(fname)
    endif
  end subroutine trackensemble_comments

  subroutine trackensemble_energy(fname,nat,nall,at,xyz,eread)
!*******************************************************
!* This subroutine decides wether to track the ensemble
!* Typically, this routine is called in CREGEN
!*******************************************************
     implicit none
     character(len=*),intent(in) :: fname
     integer,intent(in)  :: nat,nall
     integer,intent(in)  :: at(nat)
     real(wp),intent(in) :: xyz(3,nat,nall)
     real(wp),intent(in) :: eread(nall)
     character(len=50),allocatable  :: comments(:)
     integer :: i
    if (restart_tracker > restart_goal)then
      allocate(comments(nall))
      do i=1,nall
        write(comments(i),*) eread(i)
      enddo
      call restart_save_ensemble(fname,nat,nall,at,xyz,comments)
      deallocate(comments)
    if (debug) write (stdout,*) '%%% RESTART_ENSEMBLE = ',trim(last_file)
    endif
  end subroutine trackensemble_energy


  subroutine restart_write_dummy(fname)
!*******************************************************
!* This subroutine produces a placeholder file with
!* only one structure
!*******************************************************
     implicit none
     character(len=*),intent(in) :: fname
     integer :: i,ich
     if (restart_tracker < restart_goal-1)then
     if(.not.debug) write (stdout,'(a,a)') 'CREST_RESTART> writing DUMMY file ',trim(fname)
      open(newunit=ich, file=trim(fname))   
      write(ich,*) last_nat
      write(ich,*) trim(last_comments(1))
      do i=1,last_nat
      write(ich,'(a2,3f20.10)') PSE(last_at(i)),last_xyz(1:3,i,1)
      enddo
      close(ich)
      if (debug) write (stdout,*) '%%% RESTART_DUMMY = ',trim(fname)
     endif
  end subroutine restart_write_dummy


  subroutine restore_ensemble()
!*******************************************************
!* This subroutine produces a placeholder file with
!* only one structure
!*******************************************************
     implicit none
     integer :: i,ich,j,k
     character(len=:),allocatable :: fname
     character(len=:),allocatable :: atmp
     !if(index(last_file,crefile).ne.0)then
     ! call command('rm -f '//crefile//'* 2>/dev/null')
     ! fname = crefile//'_0.xyz'
     !else
      fname=trim(last_file)
     !endif

      if(.not.debug)then
         atmp = 'CREST_RESTART> RESTORING file '//trim(fname)
         k = len_trim(atmp)+2
         write (stdout,'(a,/,a,/,a)') repeat(':',k),trim(atmp),repeat(':',k) 
      endif
      open(newunit=ich, file=trim(fname))
      do j=1,last_nall
      write(ich,*) last_nat
      write(ich,*) trim(last_comments(j))
      do i=1,last_nat
      write(ich,'(a2,3f20.10)') PSE(last_at(i)),last_xyz(1:3,i,j)
      enddo
      enddo
      close(ich)
      if (debug) write (stdout,*) '%%% RESTORE_ENSEMBLE = ',trim(fname)
  end subroutine restore_ensemble





!========================================================================================!
!========================================================================================!
!> DUMP to binary routines

  subroutine dump_restart()
    implicit none
    integer :: ich,i,j,k,l
    character(len=250) :: atmp
    if (debug) write (stdout,*) '%%% RESTART DEBUG dump summary'

    !> DO NOT OVERWRITE IF WE HAVEN'T REACHED THE PREVIOUS RESTART ENTRY POINT
    if( restart_tracker < restart_goal) return

    open (newunit=ich,file='crest.restart',status='replace',form='unformatted')

    write (ich) restart_tracker
    if (debug) write (stdout,*) '%%% RESTART_TRACKER =',restart_tracker

    if (allocated(restart_env)) then
      atmp = restart_env%cmd
      write (ich) atmp
      if (debug) write (stdout,*) '%%% cmd: ',trim(atmp)

      atmp = restart_env%inputcoords
      write (ich) atmp
      if (debug) write (stdout,*) '%%% inputcoords: ',trim(atmp)

      write (ich) restart_env%eprivious
      if (debug) write (stdout,*) '%%% eprivious: ',restart_env%eprivious

      write (ich) restart_env%elowest
      if (debug) write (stdout,*) '%%% elowest: ',restart_env%elowest

      j = restart_env%ref%nat
      write(ich) j
      if (debug) write (stdout,*) '%%% ref natoms: ', j
      do i=1,j
        write(ich) restart_env%ref%at(i)
      enddo
      do i=1,j
        write(ich) restart_env%ref%xyz(1:3,i)
      enddo
   
    end if

    call dump_last_ensemble(ich)
    if (debug) write (stdout,'(1x,a,a)') '%%% ensemble: ',trim(last_file)
    if (debug) write (stdout,*) '%%% nall: ',last_nall

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
    logical :: ex
    integer,allocatable :: at(:) 
    real(wp),allocatable :: xyz(:,:)
    inquire(file='crest.restart', exist=ex)

    if(.not.ex)then
     write(stderr,'(a)') '**ERROR** while attempting to read crest.restart: file does not exist'
     error stop
    endif

    open (newunit=ich,file='crest.restart',status='old',form='unformatted')
    write(stdout,'(/,a)') repeat(":",80)
    write(stdout,'(a)') 'READING crest.restart ...'
    write(stdout,'(/,a)') '**WARNING**'
    write(stdout,'(1x,a)') "It is a user responsibility to re-use an identical job setup,"
    write(stdout,'(1x,a)') 'either via cmd or input file. The restart option only tracks'
    write(stdout,'(1x,a)') 'structure information and a non-unique restart step ID'
    write(stdout,'(a,/)') '**WARNING**'



    read (ich) restart_goal
    write(stdout,'(1x,a,i0)') 'Target restart step: ',restart_goal  

    read (ich) atmp 
    env%cmd = trim(atmp)
    write(stdout,'(1x,a,2a)') 'Previous crest cmd: "',env%cmd,'"'

    read (ich) atmp
    env%inputcoords = trim(atmp)
    write(stdout,'(1x,a,a)') 'Previous coord input file: ',env%inputcoords

    read (ich) env%eprivious
    read (ich) env%elowest
    write(stdout,'(1x,a,f20.10)') 'Previous lowest energy: ',env%elowest

    read (ich) j
    write(stdout,'(1x,a,i0,a)') 'Original input coordinates for ',j,' atoms (Angström, CMA shifted): '
    allocate(at(j))
    do i=1,j
      read(ich) at(i)
    enddo
    allocate(xyz(3,j))
    do i=1,j
      read(ich) xyzdum(1:3)
      xyz(:,i) = xyzdum(:) 
    enddo
    write(stdout,'(a5,3a16)') 'at','X','Y','Z'
    do i=1,j
      write(stdout,'(a5,3f16.8)') trim(PSE(at(i))),xyz(1:3,i)*autoaa
    enddo
    env%ref%nat = j
    call move_alloc(at, env%ref%at)
    call move_alloc(xyz, env%ref%xyz)
     

    call read_last_ensemble(ich)
    if(last_nat > 0 .and. last_nall > 0)then
      write(stdout,'(1x,a,a)') 'Last processed ensemble file: ',trim(last_file)
      write(stdout,'(1x,a,i0)') 'Number of saved structures: ',last_nall
    endif


    close (ich)
    write(stdout,'(a)') 'FINISHED READING crest.restart ...'
    write(stdout,'(a,/)') repeat(":",80)
    !stop

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
    last_nat = nat
    read(ich) nall
    last_nall = nall
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
     if(.not.saveensembles) return
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
