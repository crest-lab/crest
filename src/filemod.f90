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


!> fortran module for simple plain-text file handling

module filemod
  use iso_fortran_env,wp => real64

  implicit none

  public :: filetype
  public :: nlines
  public :: lwidth
  public :: getlarg
  public :: clearcomment

  private

!========================================================================================!
!> the File class & procedures
  type :: filetype

    integer :: nlines   !> number of lines in the file
    integer :: lwidth   !> maximum line length (used to allocate the file)
    character(:),allocatable :: f(:)  !> The entire file content

    character(:),allocatable :: filename   !> name of the file
    integer :: lcursor                     !> "cursor" postion within the file lines
    integer :: current_line = 1            !> current line within the file

  contains
    procedure :: allocate => allocate_file      !> allocate memory for file
    procedure :: deallocate => deallocate_file  !> deallocate memory
    procedure :: close => deallocate_file       !> deallocate memory
    procedure :: read => read_file              !> read the file into memory
    procedure :: init => initialize_file        !> combination of allocate + read file
    procedure :: open => initialize_file        !> combination of allocate + read file
    procedure :: print => print_file            !> print the entire file into some io channel
    procedure :: line => getline_file           !> get the n-th line of the file
    procedure :: write => write_to_file         !> append a line to the file
    procedure :: insert => insert_to_file       !> insert line into file at a given position
    procedure :: flush => flush_file            !> instantly write the file from memory to disk under 'filename'
    procedure :: flushclose => flush_file_close !> instantly write the file from memory and deallocate afterwards
    procedure :: findall => grep_all_lines      !> grep all occurences of a substring and return postions
    procedure :: grephead => grep_first_appearance
    procedure :: overwrite => replace_line      !> replace a line in the file by a new one
    procedure :: replace => replace_line
    procedure :: rename => rename_file          !> change the filename
    procedure :: clearblanks => clear_blanklines_from_file !> "clear" blanklines from file
  end type filetype

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine allocate_file(self,fname)
    implicit none
    class(filetype) :: self
    character(len=*) :: fname
    character(len=:),allocatable :: dummy
    integer :: a,b
    integer :: i
    logical :: ex
    inquire (file=fname,exist=ex)
    if (ex) then
      b = lwidth(fname)
      a = nlines(fname)
      self%lcursor = a+1
    else !-for new file
      b = 0
      a = 1
      self%lcursor = 1
    end if
    self%lwidth = b
    dummy = repeat(' ',b+5)
    self%nlines = a
    self%current_line = 1
    allocate (self%f(a),source=dummy)
    self%filename = fname
  end subroutine allocate_file

!========================================================================================!
  subroutine deallocate_file(self)
    implicit none
    class(filetype) :: self
    if (allocated(self%f)) deallocate (self%f)
  end subroutine deallocate_file

!========================================================================================!
  subroutine print_file(self,ich)
    implicit none
    class(filetype) :: self
    integer,intent(in) :: ich
    integer :: i
    logical :: ex
    do i = 1,self%nlines
      write (ich,'(a)') trim(self%f(i))
    end do
  end subroutine print_file

!========================================================================================!
  subroutine read_file(self,fname)
    implicit none
    class(filetype)  :: self
    character(len=*) :: fname
    integer :: ich
    integer :: i
    logical :: ex
    open (newunit=ich,file=fname)
    do i = 1,self%nlines
      read (ich,'(a)') self%f(i)
    end do
    close (ich)
  end subroutine read_file

!========================================================================================!
  subroutine initialize_file(self,fname)
    implicit none
    class(filetype)  :: self
    character(len=*) :: fname
    integer :: i
    logical :: ex
    call self%allocate(fname)
    inquire (file=fname,exist=ex)
    if (ex) then
      call self%read(fname)
    end if
    return
  end subroutine initialize_file

!========================================================================================!
  function getline_file(self,n)
    implicit none
    class(filetype) :: self
    character(len=:),allocatable :: getline_file
    integer :: n
    if (n .le. self%nlines.and.n .gt. 0) then
      getline_file = trim(self%f(n))
    else
      getline_file = ''
    end if
    return
  end function getline_file

!========================================================================================!
  subroutine write_to_file(self,str)
    implicit none
    class(filetype)  :: self
    character(len=*) :: str
    character(len=:),allocatable :: tmp(:)  !dummey file
    character(len=:),allocatable :: dummy
    integer :: maxwidth
    integer :: newlines
    integer :: cursor
    integer :: i
    logical :: ex
    newlines = self%nlines+1
    if (len(str) .lt. 1) then
      maxwidth = self%lwidth+5
    else
      maxwidth = max((len(str)+5),self%lwidth)
    end if
    cursor = self%lcursor
    dummy = repeat(' ',maxwidth)
    allocate (tmp(newlines),source=dummy)
    tmp(1:self%nlines) = self%f(1:self%nlines)
    tmp(cursor) = str
    self%nlines = newlines
    self%lwidth = maxwidth
    self%lcursor = cursor+1
    call move_alloc(tmp,self%f)
    return
  end subroutine write_to_file

!========================================================================================!
  subroutine insert_to_file(self,str,pos)
    implicit none
    class(filetype)  :: self
    character(len=*) :: str
    integer :: pos
    character(len=:),allocatable :: tmp(:)  !dummey file
    character(len=:),allocatable :: dummy
    integer :: maxwidth
    integer :: newlines
    integer :: cursor
    integer :: i
    logical :: ex
    newlines = self%nlines+1
    if (len(str) .lt. 1) then
      maxwidth = self%lwidth+5
    else
      maxwidth = max((len(str)+5),self%lwidth)
    end if
    if (pos .gt. self%nlines) then
      cursor = newlines
    elseif (pos .le. 1) then
      cursor = 1
    else
      cursor = pos
    end if
    dummy = repeat(' ',maxwidth)
    allocate (tmp(newlines),source=dummy)
    if (cursor .gt. 1) tmp(1:cursor-1) = self%f(1:cursor-1)
    tmp(cursor) = str
    if (cursor .lt. newlines) tmp(cursor+1:newlines) = self%f(cursor:self%nlines)
    self%nlines = newlines
    self%lwidth = maxwidth
    self%lcursor = self%lcursor+1
    call move_alloc(tmp,self%f)
    return
  end subroutine insert_to_file

!========================================================================================!
  subroutine flush_file(self)
    implicit none
    class(filetype) :: self
    integer :: ich
    integer :: i
    logical :: ex
    open (newunit=ich,file=self%filename)
    do i = 1,self%nlines
      write (ich,'(a)') trim(self%f(i))
    end do
    close (ich)
    return
  end subroutine flush_file

!========================================================================================!
  subroutine flush_file_close(self)
    implicit none
    class(filetype) :: self
    integer :: ich
    integer :: i
    logical :: ex
    open (newunit=ich,file=self%filename)
    do i = 1,self%nlines
      write (ich,'(a)') trim(self%f(i))
    end do
    close (ich)
    if (allocated(self%f)) deallocate (self%f)
    return
  end subroutine flush_file_close

!========================================================================================!
  function grep_first_appearance(self,str)
    implicit none
    class(filetype) :: self
    character(len=*) :: str
    integer :: grep_first_appearance
    integer :: j
    integer :: i
    logical :: ex
    do j = 1,self%nlines
      if (index(self%line(j),str) .ne. 0) then
        grep_first_appearance = j
        exit
      end if
    end do
    return
  end function grep_first_appearance

!========================================================================================!
! the function return value is an array of rank 1 with the dimension k. the entries are the line numbers
! in which the string str is present.
  function grep_all_lines(self,str,k)
    implicit none
    class(filetype) :: self
    character(len=*) :: str
    integer,allocatable :: grep_all_lines(:)
    integer,allocatable :: pos(:),tmp(:)
    integer,intent(out) :: k
    integer :: j
    integer :: i
    logical :: ex
    k = 0
    do j = 1,self%nlines
      if (index(self%line(j),str) .ne. 0) then
        if (k .gt. 0) then
          allocate (tmp(k+1))
          tmp(1:k) = pos(1:k)
          tmp(k+1) = j
          call move_alloc(tmp,pos)
        else
          allocate (pos(1))
          pos(1) = j
        end if
        k = k+1
      end if
    end do
    if (k == 0) then
      grep_all_lines = (/0/)
      return
    end if
    grep_all_lines = pos
    return
  end function grep_all_lines

!========================================================================================!
!replaces line K in file SELF by STR.
  subroutine replace_line(self,k,str)
    implicit none
    class(filetype) :: self
    character(len=*) :: str
    integer :: k
    character(len=:),allocatable :: tmp(:)  !dummey "file"
    character(len=:),allocatable :: dum
    integer :: maxlen
    integer :: i
    logical :: ex
    if (k .gt. self%nlines.or.k .lt. 1) return !k is not in a valid linenumber
    if (len(str) .gt. self%lwidth) then
      maxlen = len(str)
      dum = repeat(' ',maxlen+5)
      allocate (tmp(self%nlines),source=dum)
      tmp(1:self%nlines) = self%f(1:self%nlines)
      self%lwidth = maxlen
      call move_alloc(tmp,self%f)
    end if
    self%f(k) = str
    return
  end subroutine replace_line

!========================================================================================!
!replaces the internal file name by STR
  subroutine rename_file(self,str)
    implicit none
    class(filetype) :: self
    character(len=*) :: str
    integer :: i
    logical :: ex
    if (len_trim(str) .lt. 0) return
    self%filename = trim(adjustl(str))
    return
  end subroutine rename_file

!========================================================================================!
! removes all blank lines from file
  subroutine clear_blanklines_from_file(self)
    implicit none
    class(filetype)  :: self
    character(len=:),allocatable :: tmp(:)  !dummey file
    character(len=:),allocatable :: dummy
    integer :: maxwidth
    integer :: newlines
    integer :: i,k
    logical :: ex
    !--- set the stage
    newlines = self%nlines
    maxwidth = self%lwidth
    dummy = repeat(' ',maxwidth)
    allocate (tmp(newlines),source=dummy)
    k = 1
    do i = 1,self%nlines
      dummy = self%f(i)
      if (len_trim(dummy) .lt. 1) then
        newlines = newlines-1
        cycle
      end if
      tmp(k) = dummy
      k = k+1
    end do
    self%f(1:self%nlines) = tmp(1:self%nlines)
    self%nlines = newlines
    self%lcursor = newlines+1
    deallocate (tmp)
    return
  end subroutine clear_blanklines_from_file

!========================================================================================!
!========================================================================================!
!> PUBLIC subroutines from the module
!========================================================================================!
!========================================================================================!

! get the number of lines of a file
  function nlines(fname)
    implicit none
    integer :: nlines
    character(len=*) :: fname
    integer :: ich,io
    character(len=1024) :: str
    open (newunit=ich,file=fname)
    nlines = 0
    do
      read (ich,'(a)',iostat=io) str
      if (io < 0) exit !--- EOF
      nlines = nlines+1
    end do
    close (ich)
  end function nlines

!========================================================================================!
! get the maximum file width in a file
  function lwidth(fname)
    implicit none
    integer :: lwidth
    character(len=*) :: fname
    integer :: ich,io
    character(len=5000) :: str
    open (newunit=ich,file=fname)
    lwidth = 0
    do
      read (ich,'(a)',iostat=io) str
      if (io < 0) exit !--- EOF
      lwidth = max(lwidth,len_trim(str))
    end do
    close (ich)
  end function lwidth


!========================================================================================!
!get n-th element of a line (seperated by blanks)
  function getlarg(line,n)
    implicit none
    character(len=*) :: line
    character(len=:),allocatable :: getlarg
    character(len=:),allocatable :: dum
    character(len=1) :: c
    integer :: n
    integer :: lw,nt
    integer :: i
    logical :: rd
    getlarg = ''
    lw = len_trim(line)
    nt = 0
    rd = .false.
    dum = ''
    if (lw .lt. 1.or.n .lt. 1) return
    do i = 1,lw
      c = line(i:i)
      if (c .eq. ' ') then
        if (rd) then
          if (nt .eq. n) then
            getlarg = dum
            exit
          end if
          dum = ''
        end if
        rd = .false.
        cycle
      else
        dum = dum//c
        if (.not.rd) nt = nt+1
        rd = .true.
      end if
      if (rd.and.i .eq. lw) then
        getlarg = dum
      end if
    end do
    if (nt .lt. n) getlarg = ''
    return
  end function

!========================================================================================!
!> remove any comment from a given line.
!> comments are identified by the "id" character
  subroutine clearcomment(str,id)
    implicit none
    character(len=*),intent(inout) :: str
    character(len=*),intent(in) :: id
    character(len=:),allocatable :: atmp
    integer :: k
    atmp = str
    k = index(str,id)
    if (k == 1) then
      atmp = ''
    else if (k > 0) then
      atmp = str(1:k-1)
    end if
    str = trim(atmp)
    deallocate (atmp)
    return
  end subroutine clearcomment

!========================================================================================!
!========================================================================================!
end module filemod
