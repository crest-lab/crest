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

  !> ORCA simple-input runtype keywords that CREST strips when assembling
  !> an input from a user-provided short line (only EnGrad jobs are run)
  character(len=8),parameter :: orca_runtypes(*) = [character(len=8) :: &
    & 'opt     ','optts   ','copt    ','engrad  ','numgrad ','freq    ', &
    & 'numfreq ','anfreq  ','md      ','aimd    ','sp      ','energy  ', &
    & 'goat    ','irc     ']

  type :: orca_input
    character(len=:),allocatable :: cmd
    integer :: nlines = 0
    character(len=:),allocatable :: input(:)
    logical :: mpi = .false.
    integer :: maxcore = 0    !> ORCA %maxcore per core in MB (0 = unset)
    integer :: srckind = 0    !> 0=unset, 1=template file, 2=assembled from TOML
  contains
    procedure :: read => read_orca_input
    procedure :: build => build_orca_input
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
    self%srckind = 1
  end subroutine read_orca_input

!========================================================================================!

  subroutine build_orca_input(self,simple)
!***********************************************************************
!* Assemble an ORCA simple-input '!' line from a short user string
!* (given via TOML) instead of reading a template file.
!*
!* Any runtype keyword (opt/freq/md/... see orca_runtypes) and any PALn
!* keyword are stripped; ' EnGrad' is appended so CREST always drives an
!* energy+gradient job. The %pal / %maxcore blocks are added later in
!* write_orca_input from the level thread count and self%maxcore.
!***********************************************************************
    implicit none
    class(orca_input) :: self
    character(len=*),intent(in) :: simple
    character(len=:),allocatable :: work,token,low,out
    integer :: i,n,i0,k,width
    logical :: keep

    !> strip a leading '!' if the user provided one
    work = adjustl(simple)
    if (len_trim(work) > 0) then
      if (work(1:1) == '!') work = adjustl(work(2:))
    end if

    out = '!'
    n = len_trim(work)
    i = 1
    do while (i <= n)
      !> skip whitespace between tokens
      do while (i <= n .and. work(i:i) == ' ')
        i = i+1
      end do
      if (i > n) exit
      i0 = i
      do while (i <= n .and. work(i:i) /= ' ')
        i = i+1
      end do
      token = work(i0:i-1)
      low = lowercase(token)
      keep = .true.
      !> drop explicit runtype keywords
      do k = 1,size(orca_runtypes)
        if (low == trim(orca_runtypes(k))) then
          keep = .false.
          exit
        end if
      end do
      !> drop parallelization keywords (PAL2..PAL8) -> handled via threads
      if (len(low) >= 3) then
        if (low(1:3) == 'pal') keep = .false.
      end if
      if (keep) out = trim(out)//' '//trim(token)
    end do

    out = trim(out)//' EnGrad'

    width = len_trim(out)+1
    if (allocated(self%input)) deallocate (self%input)
    allocate (self%input(1),source=repeat(' ',width))
    self%input(1) = trim(out)
    self%nlines = 1
    self%mpi = .false.
    self%srckind = 2
  end subroutine build_orca_input

!========================================================================================!

  subroutine write_orca_input(self,fname,mol,chrg,mult,nthreads)
!***********************************************************************
!* Write the ORCA input file: cached template lines followed by a
!* freshly written coordinate block.
!*
!* If nthreads > 0 is given, CREST takes ownership of the parallel
!* setup: any %pal block (single- or multi-line) and simple-input
!* PALn keyword in the template is stripped, and a single
!*   %pal nprocs <nthreads> end
!* line is appended instead. Likewise, if self%maxcore > 0 any %maxcore
!* line is replaced by '%maxcore <maxcore>'. With nthreads unset (<1)
!* and maxcore unset the template is written verbatim (backward compatible).
!***********************************************************************
    implicit none
    class(orca_input),intent(in) :: self
    character(len=*),intent(in) :: fname
    type(coord),intent(in) :: mol
    integer,intent(in) :: chrg
    integer,intent(in) :: mult
    integer,intent(in),optional :: nthreads
    integer :: ich,i,j,k,l
    integer :: nt
    logical :: override,writemem,inpalblock
    character(len=:),allocatable :: line,low

    if(.not.allocated(self%input))then
      write (stderr,'(3a)') '**ERROR** Please provide an ORCA input template!'
      error stop
    endif

    nt = 0
    if (present(nthreads)) nt = nthreads
    override = (nt > 0)
    writemem = (self%maxcore > 0)

    open (newunit=ich,file=fname)
    inpalblock = .false.
    do i=1,self%nlines
      line = trim(self%input(i))
      low = adjustl(lowercase(line))
      if (override) then
        !> skip the interior/end of a multi-line %pal ... end block
        if (inpalblock) then
          if (index(low,'end') .ne. 0) inpalblock = .false.
          cycle
        end if
        !> %pal block: single-line (contains 'end') or start of a block
        if (low(1:1) .eq. '%' .and. index(low,'pal') .ne. 0) then
          if (index(low,'end') .eq. 0) inpalblock = .true.
          cycle
        end if
        !> simple-input PALn keyword on a '!' line: drop just that token
        if (low(1:1) .eq. '!') then
          j = index(low,'pal')
          if (j .ne. 0) then
            !> blank out 'pal' + up to two trailing chars (e.g. PAL8)
            k = min(j+4,len(line))
            line(j:k) = repeat(' ',k-j+1)
            line = trim(line)
          end if
        end if
      end if
      !> drop any existing %maxcore line, CREST writes its own
      if (writemem .and. low(1:1) .eq. '%' .and. index(low,'maxcore') .ne. 0) cycle
      write(ich,'(a)') trim(line)
    enddo
    if (override) then
      write(ich,'(a,1x,i0,1x,a)') '%pal nprocs',nt,'end   # set by CREST (level threads)'
    end if
    if (writemem) then
      write(ich,'(a,1x,i0,a)') '%maxcore',self%maxcore,'   # set by CREST (per core, MB)'
    end if
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
