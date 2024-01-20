!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2023 Philipp Pracht
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

!> This is a collection of various IO subroutines

module iomod
  use iso_fortran_env,wp => real64,stdout => output_unit
  use iso_c_binding
  implicit none
  public

  interface
    function mkdir(path,mode) bind(c,name="mkdir")
      use iso_c_binding
      integer(c_int) :: mkdir
      character(kind=c_char,len=1) :: path(*)
      integer(c_int16_t),value :: mode
    end function mkdir
  end interface

  interface
    integer(kind=c_int) function c_setenv(c_name,c_VALUE) bind(c,name="setenv")
      use iso_c_binding
      character(kind=c_char)   :: c_name(*)
      character(kind=c_char)   :: c_VALUE(*)
    end function
  end interface
  interface setenv
    module procedure setenv_string
    module procedure setenv_int
    module procedure setenv_float
  end interface setenv

  interface
    integer(kind=c_int) function c_symlink(c_from,c_to) bind(c,name="symlink")
      use iso_c_binding
      character(kind=c_char)   :: c_from(*)
      character(kind=c_char)   :: c_to(*)
    end function
  end interface

  interface wrshort
    module procedure wrshort_real
    module procedure wrshort_int
    module procedure wrshort_string
  end interface wrshort

  interface rdshort
    module procedure rdshort_real
    module procedure rdshort_int
    module procedure rdshort_string
  end interface rdshort

  interface uppercase
!    module procedure to_upper
    module procedure to_upper_func
  end interface uppercase

  interface lowercase
!    module procedure to_lower
    module procedure to_lower_func
  end interface lowercase

  interface convert_to_string
    module procedure i8_to_string
    module procedure i16_to_string
    module procedure i32_to_string
    module procedure i64_to_string
    module procedure r32_to_string
    module procedure r64_to_string
    module procedure bool_to_string
  end interface convert_to_string
  interface to_str
    module procedure to_str_i8
    module procedure to_str_i16
    module procedure to_str_i32
    module procedure to_str_i64
    module procedure to_str_r32
    module procedure to_str_r64
    module procedure to_str_bool
  end interface to_str

  public :: checkprog
  public :: checkprog_silent
  private :: getpath

!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!
!-------------------------------------------------------------
! remove file 'fname'
! if 'fname' contains the wildcard '*', the remove is ignored
!--------------------------------------------------------------
  subroutine remove(fname)
    implicit none
    character(len=*) :: fname
    integer :: ich,och
    logical :: ex
    if (index(fname,'*') .eq. 0) then
      open (newunit=ich,file=fname)
      close (ich,status='delete')
    end if
  end subroutine remove

!-------------------------------------------------------------
! creates an empty file "fname", similar to the shell command "touch"
!--------------------------------------------------------------
  subroutine touch(fname)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    open (newunit=ich,file=fname)
    close (ich)
  end subroutine touch

!-------------------------------------------------------------
! prints the file "fname", similar to the shell command "cat"
!--------------------------------------------------------------
  subroutine cat(fname)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    character(len=512) :: str
    integer :: io
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) str
      if (io < 0) exit
      write (stdout,'(a)') trim(str)
    end do
    close (ich)
  end subroutine cat

!--------------------------------------------------------------------------------------
! prints the file "fname", similar to the shell command "cat" but deletes it at the end
!--------------------------------------------------------------------------------------
  subroutine catdel(fname)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    character(len=512) :: str
    integer :: io
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) str
      if (io < 0) exit
      write (stdout,'(a)') trim(str)
    end do
    close (unit=ich,status='delete')
  end subroutine catdel

!-------------------------------------------------------------------------
! copy a file from path "from" to path "to"
!-------------------------------------------------------------------------
  subroutine copy(from,to)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: from
    character(len=*) :: to
    character(len=512) :: line
    integer :: io
    inquire (file=from,exist=ex)
    if (ex) then
      open (newunit=ich,file=from)
      open (newunit=och,file=to)
      do
        read (ich,'(a)',iostat=io) line
        if (io < 0) exit
        write (och,'(a)') trim(line)
      end do
      close (och)
      close (ich)
    end if
  end subroutine copy

!-------------------------------------------------------------------------
! copy a file from path "from" to a specified sub-directory "to"
!-------------------------------------------------------------------------
  subroutine copysub(from,to)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: from
    character(len=*) :: to
    character(len=512) :: line,dir1
    integer :: io
    inquire (file=from,exist=ex)
    if (ex) then
      call getcwd(dir1)
      inquire (file=from,exist=ex)
      if (ex) then
        open (newunit=ich,file=from)
        call chdir(to)
        open (newunit=och,file=from)
        do
          read (ich,'(a)',iostat=io) line
          if (io < 0) exit
          write (och,'(a)') trim(line)
        end do
        close (och)
        close (ich)
      end if
      call chdir(dir1)
    end if
  end subroutine copysub

!-------------------------------------------------------------------------
! move a file from path "from" to path "to". "from" file is deleted.
!-------------------------------------------------------------------------
  subroutine move(from,to)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: from
    character(len=*) :: to
    inquire (file=from,exist=ex)
    if (ex) then
      call remove(to)
      call rename(from,to)
    else
      write (stdout,'(''file '',a,'' does not exist!'')') trim(from)
    end if
  end subroutine move

!-------------------------------------------------------------------------
! converts a path relative to the current directory (i.e. containing "../"
! and such) into the absolute Path
!-------------------------------------------------------------------------
  subroutine getAbsPath(Path)
    implicit none
    character(len=*),intent(inout) :: Path
    character(len=512) :: thispath
    call getcwd(thispath)  !get current working directory
    call chdir(Path)
    call getcwd(Path)  !get the absolute path which was given as a relative path
    call chdir(thispath) !switch back to original directory
  end subroutine getAbsPath

!-------------------------------------------------------------------------
! the "inquire" command, but in a sub-directory
!-------------------------------------------------------------------------
  subroutine inquiresub(f,sub,bool)
    implicit none
    character(len=*),intent(in) :: f
    character(len=*),intent(in) :: sub
    character(len=512) :: thispath
    logical,intent(out) :: bool
    call getcwd(thispath)  !get current working directory
    call chdir(sub) !go to the subdirectory
    inquire (file=f,exist=bool) !do the "inquire"
    call chdir(thispath) !switch back to original directory
  end subroutine inquiresub

!-------------------------------------------------------------------------
! append content of text file "from" into text file "to", similar to "cat A >> B"
!-------------------------------------------------------------------------
  subroutine appendto(from,to)
    implicit none
    integer :: ich,och
    logical :: ex
    integer :: io
    character(len=*) :: from
    character(len=*) :: to
    character(len=1024) :: str
    open (newunit=ich,file=to)
    open (newunit=och,file=from)
    do
      read (ich,*,iostat=io)
      if (io < 0) then
        backspace (ich)
        exit
      end if
    end do
    do
      read (och,'(a)',iostat=io) str
      if (io < 0) then
        exit
      else
        write (ich,'(a)') trim(str)
      end if
    end do
    close (och)
    close (ich)
  end subroutine appendto

!-------------------------------------------------------------------------
! make a directory via iso_c_binding
!-------------------------------------------------------------------------
  function makedir(str)
    implicit none
    integer :: makedir
    character(len=*) :: str
    makedir = mkdir(str//char(0),int(o'770',c_int16_t)) !create new directory
    return
  end function makedir

!-------------------------------------------------------------------------
! put a environment variable via iso_c_binding
!-------------------------------------------------------------------------
  function setenv_string(env,str)
    implicit none
    integer :: setenv_string
    character(len=*) :: env
    character(len=*) :: str
    setenv_string = c_setenv(env//c_null_char,str//c_null_char) !create new directory
    return
  end function setenv_string
  function setenv_int(env,intval)
    implicit none
    integer :: setenv_int
    character(len=*) :: env
    integer :: intval
    character(len=20) :: str
    write (str,'(i0)') intval
    setenv_int = c_setenv(env//c_null_char,trim(str)//c_null_char) !create new directory
    return
  end function setenv_int
  function setenv_float(env,floatval)
    implicit none
    integer :: setenv_float
    character(len=*) :: env
    real(wp) :: floatval
    character(len=30) :: str
    write (str,'(f14.6)') floatval
    setenv_float = c_setenv(env//c_null_char,trim(str)//c_null_char) !create new directory
    return
  end function setenv_float

!-------------------------------------------------------------------------
! set a symlink from path1 to path2  via iso_c_binding
!-------------------------------------------------------------------------
  function sylnk(path1,path2)
    implicit none
    integer :: sylnk
    character(len=*) :: path1
    character(len=*) :: path2
    sylnk = c_symlink(trim(path1)//c_null_char,trim(path2)//c_null_char) !create new directory
    return
  end function sylnk

!-------------------------------------------------------------------------
! minigrep: a grep subroutine that returns true or false, depending on
! if the substring "str" is present in the file "fil"
!-------------------------------------------------------------------------
  subroutine minigrep(fil,str,bool)
    implicit none
    integer :: ich,och
    logical :: ex
    logical :: bool
    character(len=*) :: fil
    character(len=*) :: str
    character(len=512) :: tmp
    integer :: io
    bool = .false.
    open (newunit=ich,file=fil)
    do
      read (ich,'(a)',iostat=io) tmp
      if (io < 0) exit
      if (index(tmp,str) .ne. 0) then
        bool = .true.
        close (ich)
        exit
      end if
    end do
    close (ich)
    !minigrep = bool
    return
  end subroutine minigrep

!-------------------------------------------------------------------------
! grepval: a grep subroutine that returns the float that follows immeadeatly
! after the substring "str". Works only for the first occurence of "str".
! If str is not present, "bool" returns .false.
!-------------------------------------------------------------------------
  subroutine grepval(fil,str,bool,val)
    implicit none
    integer :: ich,och
    logical :: ex
    logical :: bool
    character(len=*) :: fil
    character(len=*) :: str
    character(len=512) :: tmp
    real(wp) :: val
    integer :: io
    bool = .false.
    val = 0.0d0
    open (newunit=ich,file=fil)
    do
      read (ich,'(a)',iostat=io) tmp
      if (io < 0) exit
      io = index(tmp,str,.true.)
      if (io .ne. 0) then
        io = io+len(str)
        tmp = adjustl(tmp(io+1:))
        read (tmp,*) val
        bool = .true.
        close (ich)
        exit
      end if
    end do
    close (ich)
    !minigrep = bool
    return
  end subroutine grepval

!---------------------------------------------------------------------------
! grepcntx: a grep subroutine that returns the entire line containing "str"
! or the specified line after the occurence of "str", given
! If str is not present, "bool" returns .false.
!---------------------------------------------------------------------------
  subroutine grepcntxt(fil,str,bool,line,context)
    implicit none
    integer :: ich,och
    logical :: ex
    logical,intent(out) :: bool
    logical :: track
    character(len=*),intent(in) :: fil
    character(len=*),intent(in) :: str
    character(len=*),intent(out) :: line
    integer,intent(in),optional :: context
    integer :: cntxt
    character(len=512) :: tmp
    real(wp) :: val
    integer :: io,itrack,reftrack
    bool = .false.
    track = .false.
    val = 0.0d0
    line = ''
    itrack = 0
    if (present(context)) then
      reftrack = max(0,context)
    else
      reftrack = 0
    end if
    open (newunit=ich,file=fil)
    do
      read (ich,'(a)',iostat=io) tmp
      if (io < 0) exit
      io = index(tmp,str,.true.)
      if (io .ne. 0) then
        bool = .true.
        track = .true.
      end if
      if (track) then
        if (itrack == reftrack) then
          line = trim(tmp)
          exit
        else
          itrack = itrack+1
        end if
      end if
    end do
    close (ich)
    return
  end subroutine grepcntxt

!=====================================================================================!
  function getlines(fname,maxwidth) result(n)
    character(len=*) :: fname
    logical :: ex
    integer,intent(out),optional :: maxwidth
    integer :: n,io,w,ich
    character(len=1056) :: atmp
    n = 0
    w = 0
    inquire (file=fname,exist=ex)
    if (.not.ex) return
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) atmp
      if (io /= 0) exit
      n = n+1
      if (len_trim(atmp) > w) w = len_trim(atmp)
    end do
    if (present(maxwidth)) maxwidth = w
    close (ich)
  end function getlines

!--------------------------------------------------------------------------------------
! Write a file with a single line (INT)
!--------------------------------------------------------------------------------------
  subroutine wrshort_int(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    integer :: var
    open (newunit=ich,file=fname)
    write (ich,*) var
    close (ich)
  end subroutine wrshort_int

!--------------------------------------------------------------------------------------
! Write a file with a single line (REAL)
!--------------------------------------------------------------------------------------
  subroutine wrshort_real(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    real(wp) :: var
    open (newunit=ich,file=fname)
    write (ich,*) var
    close (ich)
  end subroutine wrshort_real

!--------------------------------------------------------------------------------------
! Write a file with a single line (STRING)
!--------------------------------------------------------------------------------------
  subroutine wrshort_string(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    character(len=*) :: var
    open (newunit=ich,file=fname)
    write (ich,'(a)') var
    close (ich)
  end subroutine wrshort_string

!--------------------------------------------------------------------------------------
! reads the first line of a file (INT)
!--------------------------------------------------------------------------------------
  subroutine rdshort_int(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    integer :: var
    integer :: io
    open (newunit=ich,file=fname)
    read (ich,*,iostat=io) var
    if (io .ne. 0) var = 0
    close (ich)
  end subroutine rdshort_int

!--------------------------------------------------------------------------------------
! reads the first line of a file (REAL)
!--------------------------------------------------------------------------------------
  subroutine rdshort_real(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    real(wp) :: var
    integer :: io
    open (newunit=ich,file=fname)
    read (ich,*,iostat=io) var
    if (io .ne. 0) var = 0.0d0
    close (ich)
  end subroutine rdshort_real

!--------------------------------------------------------------------------------------
! reads the first line of a file (STRING)
!--------------------------------------------------------------------------------------
  subroutine rdshort_string(fname,var)
    implicit none
    integer :: ich,och
    logical :: ex
    character(len=*) :: fname
    character(len=*) :: var
    integer :: io
    open (newunit=ich,file=fname)
    read (ich,'(a)',iostat=io) var
    if (io .ne. 0) var = ''
    close (ich)
  end subroutine rdshort_string

!----------------------------------------------------------------------------
! returns entry x of line as string, or as empty if end-of-line is reached

  subroutine glinex(line,x,string)
    implicit none
    integer :: ich,och
    logical :: ex
    integer,intent(in) :: x
    character(len=*),intent(in)  :: line
    character(len=*),intent(out) :: string
    character(len=80) :: stmp
    character(len=1)  :: digit
    integer :: i,k,cs
    logical :: track

    string = ''
    stmp = ''
    cs = 0
    track = .false.
    k = len(trim(line))
    do i = 1,k
      digit = line(i:i)
      if (digit .eq. ' '.or.digit .eq. char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
        if (track) then
          if (cs == x) then
            string = trim(stmp)
            exit
          end if
        end if
        track = .false.
        stmp = ''
        cycle
      else
        if (.not.track) cs = cs+1
        track = .true.
        stmp = trim(stmp)//trim(digit)
        if (i .eq. k.and.x .lt. 0) then  !if x is lt. 0, the last element of the line is returned
          string = trim(stmp)
          exit
        end if
      end if
    end do
    return
  end subroutine glinex

!----------------------------------------------------------------------------
! cut the leading x words (seperated by blanks) of a string, and overwrite it

  subroutine clinex(line,x)
    implicit none
    integer,intent(in) :: x
    character(len=*),intent(inout)  :: line
    character(len=1)  :: digit
    integer :: i,j,k,cs
    logical :: track

    cs = 0
    track = .false.
    k = len(trim(line))
    do i = 1,k
      digit = line(i:i)
      if (digit .eq. ' '.or.digit .eq. char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
        track = .false.
        cycle
      else
        if (.not.track) cs = cs+1
        if (cs .gt. x) then
          j = i-1
          line = line(j:k)
          exit
        end if
        track = .true.
      end if
    end do
    return
  end subroutine clinex

!-------------------------------------------------------------
! convert a string to upper or lower case
!--------------------------------------------------------------
  subroutine to_upper(str)
    character(*),intent(in out) :: str
    integer :: i
    do i = 1,len(str)
      select case (str(i:i))
      case ("a":"z")
        str(i:i) = achar(iachar(str(i:i))-32)
      end select
    end do
  end subroutine to_upper
  function to_upper_func(str) result(out)
    character(*),intent(in) :: str
    character(len=:),allocatable :: out
    out = str
    call to_upper(out)
  end function to_upper_func
  subroutine to_lower(str)
    character(*),intent(in out) :: str
    integer :: i
    do i = 1,len(str)
      select case (str(i:i))
      case ("A":"Z")
        str(i:i) = achar(iachar(str(i:i))+32)
      end select
    end do
  end subroutine to_lower
  function to_lower_func(str) result(out)
    character(*),intent(in) :: str
    character(len=:),allocatable :: out
    out = str
    call to_lower(out)
  end function to_lower_func

!-----------------------------------------------------------------
! a simple routine to lock a feature of the code behind a password
!-----------------------------------------------------------------
  subroutine pwdlock(pwdref)
    implicit none
    character(len=*) :: pwdref    ! check for this password
    character(len=:),allocatable :: pwd
    pwd = pwdref !allocate max length of pwd
    write (stdout,'(/,1x,a)',advance='no') 'Locked feature. Enter password: '
    read (*,*) pwd
    if (pwd == pwdref) then
      write (stdout,'(1x,a,/)') 'Valid. Continue.'
    else
      error stop 'Invalid. Stop.'
    end if
    return
  end subroutine pwdlock

!====================================================================!
  function filechecker(fin,fout) result(have)
    implicit none
    logical :: have
    character(len=*),intent(in) :: fin
    character(len=*),intent(out) :: fout
    have = .false.
    inquire (file=fin,exist=have)
    if (have) then
      fout = fin
    else
      fout = ''
    end if
    return
  end function filechecker

!=========================================================================================!
!=========================================================================================!
!=========================================================================================!

!> Checking directories with Fortran's inquire is not handled by the standard.
!> The interpretation whether or not to report a directory as file is compiler
!> specific and therefore always an extension to the Fortran standard.
  function directory_exist(file) result(exist)
    character(len=*),intent(in) :: file
    logical :: exist
#ifdef __INTEL_COMPILER
    !> Intel provides the directory extension to inquire to handle this case
    inquire (directory=file,exist=exist)
#else
    !> GCC handles directories as files, to make sure we get a directory and
    !> not a file append a path separator and the current dir
    inquire (file=trim(file)//"/.",exist=exist)
#endif
  end function directory_exist

!=========================================================================================!
!=========================================================================================!
!=========================================================================================!

!> a wrapper for the intrinsic isatty function.
!> ifort only seems to work if isatty is declard as external
!> while gfortran does not want that...
  function myisatty(channel) result(term)
    implicit none
    integer,intent(in) :: channel
    logical :: term
#ifdef __INTEL_COMPILER
    logical,external :: isatty
#endif
    term = isatty(channel) 
  end function myisatty

!=========================================================================================!
!=========================================================================================!
!=========================================================================================!
!> type conversion routines
!> 8 bit integers
  subroutine i8_to_string(i_in,str_out)
    implicit none
    integer(kind=int8) :: i_in
    character(len=:),allocatable :: str_out
    character(len=4) :: stmp
    write (stmp,*) i_in
    str_out = trim(adjustl(stmp))
  end subroutine i8_to_string

!> 16 bit integers
  subroutine i16_to_string(i_in,str_out)
    implicit none
    integer(kind=int16) :: i_in
    character(len=:),allocatable :: str_out
    character(len=6) :: stmp
    write (stmp,*) i_in
    str_out = trim(adjustl(stmp))
  end subroutine i16_to_string

!> 32 bit integer
  subroutine i32_to_string(i_in,str_out)
    implicit none
    integer(kind=int32) :: i_in
    character(len=:),allocatable :: str_out
    character(len=12) :: stmp
    write (stmp,*) i_in
    str_out = trim(adjustl(stmp))
  end subroutine i32_to_string

!> 64 bit integer
  subroutine i64_to_string(i_in,str_out)
    implicit none
    integer(kind=int64) :: i_in
    character(len=:),allocatable :: str_out
    character(len=60) :: stmp
    write (stmp,*) i_in
    str_out = trim(adjustl(stmp))
  end subroutine i64_to_string

!> 32 bit float
  subroutine r32_to_string(r_in,str_out)
    implicit none
    real(kind=real32) :: r_in
    character(len=:),allocatable :: str_out
    character(len=80) :: stmp
    write (stmp,*) r_in
    str_out = trim(adjustl(stmp))
  end subroutine r32_to_string

!> 64 bit float
  subroutine r64_to_string(r_in,str_out)
    implicit none
    real(kind=wp) :: r_in
    character(len=:),allocatable :: str_out
    character(len=80) :: stmp
    write (stmp,*) r_in
    str_out = trim(adjustl(stmp))
  end subroutine r64_to_string

!> boolean
  subroutine bool_to_string(bool,str_out)
    implicit none
    logical :: bool
    character(len=:),allocatable :: str_out
    if (bool) then
      str_out = 'true'
    else
      str_out = 'false'
    end if
  end subroutine bool_to_string

!> function variants
  function to_str_i8(i8) result(str_out)
    implicit none
    integer(kind=int8) :: i8
    character(len=:),allocatable :: str_out
    call convert_to_string(i8,str_out)
  end function to_str_i8
  function to_str_i16(i16) result(str_out)
    implicit none
    integer(kind=int16) :: i16
    character(len=:),allocatable :: str_out
    call convert_to_string(i16,str_out)
  end function to_str_i16
  function to_str_i32(i32) result(str_out)
    implicit none
    integer(kind=int32) :: i32
    character(len=:),allocatable :: str_out
    call convert_to_string(i32,str_out)
  end function to_str_i32
  function to_str_i64(i64) result(str_out)
    implicit none
    integer(kind=int64) :: i64
    character(len=:),allocatable :: str_out
    call convert_to_string(i64,str_out)
  end function to_str_i64
  function to_str_r32(r32) result(str_out)
    implicit none
    real(kind=real32) :: r32
    character(len=:),allocatable :: str_out
    call convert_to_string(r32,str_out)
  end function to_str_r32
  function to_str_r64(r64) result(str_out)
    implicit none
    real(kind=wp) :: r64
    character(len=:),allocatable :: str_out
    call convert_to_string(r64,str_out)
    str_out = truncate_zeros(str_out)
  end function to_str_r64
  function to_str_bool(bool) result(str_out)
    implicit none
    logical :: bool
    character(len=:),allocatable :: str_out
    call convert_to_string(bool,str_out)
  end function to_str_bool

!> truncate zeros of a float string
  function truncate_zeros(str_in) result(str_out)
    character(len=:),allocatable :: str_in
    character(len=:),allocatable :: str_out
    integer :: i,j,k,l
    str_out = ''
    if (allocated(str_in)) then
      if (index(lowercase(str_in),'e') .ne. 0) then
        str_out = trim(adjustl(str_in))
        return
      end if
      l = len_trim(str_in)
      k = l
      do i = l,1,-1
        k = i
        if (str_in(i:i) .ne. '0') then
          exit
        end if
      end do
      if (str_in(k:k) .eq. '.') k = k+1
      str_out = str_in(1:k)
      str_out = trim(adjustl(str_out))
    end if
  end function truncate_zeros

!========================================================================================!
!========================================================================================!
!========================================================================================!

  subroutine getpath(fname,path)
    implicit none
    character(len=*) :: fname
    character(len=*) :: path
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: rcode,ich,io

    pipe = ' >/dev/null 2>/dev/null'

    checkcall = 'command -v '//trim(fname)//pipe
    call command(checkcall,rcode)

    if (rcode .eq. 0) then
      checkcall = 'command -v '//trim(fname)//' > pathout.tmp 2>/dev/null'
      call command(checkcall,io)
      open (newunit=ich,file='pathout.tmp')
      read (ich,'(a)') path
      close (ich,status='delete')
    end if
    return
  end subroutine getpath
!=========================================================================================!
  subroutine checkprog(fname,r)
    implicit none
    character(len=*) :: fname
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: rcode,r
    character(len=512) :: path

    pipe = ' >/dev/null 2>/dev/null'

    checkcall = 'command -v '//trim(fname)//pipe
    call command(checkcall,rcode)

    write (stdout,'(4x,a,a,a)') 'binary: "',trim(fname),'"'
    if (rcode .ne. 0) then
      write (stdout,'(4x,a)') 'status: not found'
      r = r+1
    else
      write (stdout,'(4x,a)') 'status: present'
      call getpath(fname,path)
      write (stdout,'(4x,a,a)') 'path  : ',trim(path)
    end if

    return
  end subroutine checkprog
!=========================================================================================!
! the same as above but only provide a printout (at all) if the program is not present!
  subroutine checkprog_silent(fname,verbose,iostat)
    implicit none
    character(len=*) :: fname
    logical,intent(in),optional :: verbose
    integer,intent(out),optional :: iostat
    character(len=:),allocatable :: checkcall
    character(len=:),allocatable :: pipe
    integer :: io

    pipe = ' >/dev/null 2>/dev/null'
    checkcall = 'command -v '//trim(fname)//pipe
    call command(trim(checkcall),exitstat=io)

    if (present(verbose)) then
      if (io .ne. 0.and.verbose) then
        write (stdout,'(4x,a,a,a)') 'binary: "',trim(fname),'"'
        write (stdout,'(4x,a)') 'status: not found!'
      end if
    end if

    if (present(iostat)) then
      iostat = io
    end if

    return
  end subroutine checkprog_silent

!========================================================================================!
!========================================================================================!
!========================================================================================!

!> For some reason the behaviour of "call system"  and "call execute_command_line"
!> differs slightly beteen Intel and GNU versions of the program that I have build.
!> I don't understand why this is the case.
!> To circumvent issues for now, I will simply include two different wrappers with
!> precompiler statements.
!> It doesn't matter too much as these are mostly relevant for legacy routines.
  subroutine command(cmd,exitstat)
    implicit none
    character(len=*),intent(in) :: cmd
    integer,intent(out),optional :: exitstat
    integer :: io
#ifdef __INTEL_COMPILER
    call execute_command_line(trim(cmd),exitstat=io)
#else
    call system(trim(cmd),io)
#endif
    if (present(exitstat)) exitstat = io
  end subroutine command

!========================================================================================!
!========================================================================================!
!========================================================================================!
end module iomod

