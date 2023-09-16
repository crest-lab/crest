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

!> module orca_type
!> a minimal implementation for storing an ORCA input file

module orca_type
  use crest_parameters,only:wp,stdout,stderr,autoaa
  use iomod
  use strucrd
  implicit none
  public

  type :: orca_input
    character(len=:),allocatable :: cmd
    integer :: nlines = 0
    character(len=:),allocatable :: input(:)
    logical :: mpi = .false.
  contains
    procedure :: read => read_orca_input
    procedure :: write => write_orca_input
  end type orca_input

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine read_orca_input(self,fname)
    implicit none
    class(orca_input) :: self
    character(len=*),intent(in) :: fname
    logical :: ex,trackcoord
    integer :: nlines,width
    integer :: io,ich,i,j,k,l
    character(len=1056) :: atmp
    logical :: gotengrad

    inquire (file=fname,exist=ex)
    if (.not.ex) then
      write (stderr,'(3a)') '**ERROR** ORCA input template ',fname,' could not be found!'
      error stop
    else
      write (stdout,'(3a)') 'Reading ORCA input template ',fname,' ...'
    end if
    nlines = getlines(fname,width)
    open (newunit=ich,file=fname)
    nlines = 0
    trackcoord = .false.
    gotengrad = .false.
    !> count lines and check keywords
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit  !> EOF
      atmp = adjustl(lowercase(atmp))

      !> ignore coord lines, CREST will write those
      if (trackcoord.and.atmp(1:1) == '*') then
        trackcoord = .false.
        cycle
      end if
      if (atmp(1:1).eq.'*' .and. index(atmp,'xyz ') .ne. 0) then
        trackcoord = .true.
        cycle
      end if
      if (atmp(1:1).eq.'*' .and. index(atmp,'xyzfile') .ne. 0) cycle

      if (.not.trackcoord) nlines = nlines+1

      if (index(atmp,'$new_job') .eq. 1) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': please define only single jobs (no $new_job)!'
        error stop
      end if

      !> check runtypes from the simple input line
      if (atmp(1:1) .eq. '!'.and.index(atmp,'engrad') .ne. 0) then
        gotengrad = .true.
      end if
      if (atmp(1:1) .eq. '!'.and.index(atmp,'md') .ne. 0) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': only run EnGrad jobs!'
        error stop
      end if
      if (atmp(1:1) .eq. '!'.and.index(atmp,'opt') .ne. 0) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': only run EnGrad jobs!'
        error stop
      end if
      if (atmp(1:1) .eq. '!'.and.index(atmp,'freq') .ne. 0) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': only run EnGrad jobs!'
        error stop
      end if

      !> coordinate input check, e.g.  remove Bohr keyword, if necessary
      if (atmp(1:1) .eq. '%'.and.index(atmp,'coords') .ne. 0) then
        write (stderr,'(3a)') '**ERROR** In ORCA input template ',fname,': please remove %coords block!'
        error stop
      end if

      !> check if there is parallelization
      if (atmp(1:1) .eq. '!'.and.index(atmp,'pal') .ne. 0) then
        self%mpi = .true.
      end if
      if (atmp(1:1) .eq. '%'.and.index(atmp,'pal') .ne. 0) then
        self%mpi = .true.
      end if
    end do
    close (ich)

    !> allocate memory space
    width = width+10
    allocate (self%input(nlines),source=repeat(" ",width))
    self%nlines = nlines

    !> Open file from the beginning and read into memory
    open (newunit=ich,file=fname)
    k = 0
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit  !> EOF
      atmp = adjustl(atmp)

      !> ignore coord lines, CREST will write those
      if (trackcoord.and.atmp(1:1) == '*') then
        trackcoord = .false.
        cycle
      end if
      if (atmp(1:1).eq.'*' .and. index(atmp,'xyz ') .ne. 0) then
        trackcoord = .true.
        cycle
      end if
      if (atmp(1:1).eq.'*' .and. index(atmp,'xyzfile') .ne. 0) cycle

      if (atmp(1:1) .eq. '!'.and.index(atmp,'bohr') .ne. 0) then
        j = index(atmp,'bohr')
        atmp(j:j+3) = '    '
      end if

      k = k+1
      if (atmp(1:1) .eq. '!'.and..not.gotengrad) then
        atmp = trim(atmp)//' EnGrad'
        gotengrad = .true.
      end if
      self%input(k) = trim(atmp)

    end do
    close (ich)
  end subroutine read_orca_input

!========================================================================================!

  subroutine write_orca_input(self,fname,mol,chrg,mult)
    implicit none
    class(orca_input),intent(in) :: self
    character(len=*),intent(in) :: fname
    type(coord),intent(in) :: mol
    integer,intent(in) :: chrg
    integer,intent(in) :: mult
    integer :: ich,i,j,k,l

    if(.not.allocated(self%input))then
      write (stderr,'(3a)') '**ERROR** Please provide an ORCA input template!'
      error stop
    endif

    open (newunit=ich,file=fname)
    do i=1,self%nlines
      write(ich,'(a)') trim(self%input(i))
    enddo
    write(ich,*)
    write(ich,'(a,1x,i0,1x,i0,a)') '*xyz',chrg,mult,'  # charge and multiplicity (2S+1)'
    do i=1,mol%nat
      write(ich,'(a2,3F25.15)') asym(mol%at(i)),mol%xyz(1:3,i)*autoaa
    enddo
    write(ich,'("*")')
    close (ich)

  end subroutine write_orca_input

!========================================================================================!
!========================================================================================!
end module orca_type
