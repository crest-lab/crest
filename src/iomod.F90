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
module iomod
      use iso_fortran_env, only : wp => real64
      use iso_c_binding
      implicit none
      integer :: i,j,k,l,ich,och
      logical :: ex

      private :: wp
      private :: i,j,k,l,ich,och
      private :: ex

      interface
        function mkdir(path,mode) bind(c,name="mkdir")
          use iso_c_binding
          integer(c_int) :: mkdir
          character(kind=c_char,len=1) :: path(*)
          integer(c_int16_t), value :: mode
        end function mkdir
      end interface


      interface
         integer(kind=c_int) function c_setenv(c_name,c_VALUE) bind(c,name="setenv")
         use iso_c_binding
         import c_int, c_char
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
         import c_int, c_char
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


      contains
!-------------------------------------------------------------
! remove file 'fname'
! if 'fname' contains the wildcard '*', the remove is ignored
!--------------------------------------------------------------
subroutine remove(fname)
      implicit none
      character(len=*) :: fname
      if(index(fname,'*').eq.0)then
       open(newunit=ich,file=fname)
       close(ich,status='delete')
      endif
end subroutine remove

!-------------------------------------------------------------
! creates an empty file "fname", similar to the shell command "touch"
!--------------------------------------------------------------
subroutine touch(fname)
      implicit none
      character(len=*) :: fname
      open(newunit=ich,file=fname)
      close(ich)
end subroutine touch

!-------------------------------------------------------------
! prints the file "fname", similar to the shell command "cat"
!--------------------------------------------------------------
subroutine cat(fname)
      implicit none
      character(len=*) :: fname
      character(len=512) :: str
      integer :: io
      open(newunit=ich,file=fname)
       do
        read(ich,'(a)',iostat=io)str
        if (io < 0) exit
        write(*,'(a)')trim(str)
       enddo
      close(ich)
end subroutine cat

!--------------------------------------------------------------------------------------
! prints the file "fname", similar to the shell command "cat" but deletes it at the end
!--------------------------------------------------------------------------------------
subroutine catdel(fname)
      implicit none
      character(len=*) :: fname
      character(len=512) :: str
      integer :: io
       open(newunit=ich,file=fname)
        do
         read(ich,'(a)',iostat=io)str
         if (io < 0) exit
         write(*,'(a)')trim(str)
        enddo
       close(unit=ich,status='delete')
end subroutine catdel

!-------------------------------------------------------------------------
! copy a file from path "from" to path "to"
!-------------------------------------------------------------------------
subroutine copy(from,to)
      implicit none
      character(len=*) :: from
      character(len=*) :: to
      character(len=512) :: line      
      integer :: io
      inquire(file=from,exist=ex)
      if(ex) then
        open(newunit=ich,file=from)
        open(newunit=och,file=to)
        do
           read(ich,'(a)',iostat=io)line
           if (io < 0) exit
           write(och,'(a)')trim(line)
        enddo
        close(och)
        close(ich)
      endif
end subroutine copy

!-------------------------------------------------------------------------
! copy a file from path "from" to a specified sub-directory "to"
!-------------------------------------------------------------------------
subroutine copysub(from,to)
      implicit none
      character(len=*) :: from
      character(len=*) :: to
      character(len=512) :: line,dir1
      integer :: io
      inquire(file=from,exist=ex)
      if(ex)then
      call getcwd(dir1)
      inquire(file=from,exist=ex)
      if(ex) then
        open(newunit=ich,file=from)
        call chdir(to)
        open(newunit=och,file=from)
        do
           read(ich,'(a)',iostat=io)line
           if (io < 0) exit
           write(och,'(a)')trim(line)
        enddo
        close(och)
        close(ich)
      endif
      call chdir(dir1)
      endif
end subroutine copysub

!-------------------------------------------------------------------------
! move a file from path "from" to path "to". "from" file is deleted.
!-------------------------------------------------------------------------
 subroutine move(from,to)
      implicit none
      character(len=*) :: from
      character(len=*) :: to
      integer :: io
      inquire(file=from,exist=ex)
      if(ex) then
        call remove(to)
        call rename(from,to)
      else
        write(*,'(''file '',a,'' does not exist!'')')trim(from)
      endif
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
      inquire(file=f,exist=bool) !do the "inquire"
      call chdir(thispath) !switch back to original directory
end subroutine inquiresub

!-------------------------------------------------------------------------
! append content of text file "from" into text file "to", similar to "cat A >> B"
!-------------------------------------------------------------------------
subroutine appendto(from,to)
      implicit none
      integer :: io
      character(len=*) :: from
      character(len=*) :: to
      character(len=1024) :: str
      open(newunit=ich,file=to)
      open(newunit=och,file=from)
      do
        read(ich,*,iostat=io)
        if(io<0)then
            backspace(ich)
            exit
        endif
      enddo
      do
        read(och,'(a)',iostat=io)str
        if(io<0)then
            exit
        else
            write(ich,'(a)')trim(str)
        endif
      enddo
      close(och)
      close(ich)
end subroutine appendto

!-------------------------------------------------------------------------
! make a directory via iso_c_binding
!-------------------------------------------------------------------------
 function makedir(str)
      implicit none
      integer :: makedir
      character(len=*) :: str
      makedir = mkdir(str//char(0), int(o'770',c_int16_t)) !create new directory
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
      write(str,'(i0)')intval
      setenv_int = c_setenv(env//c_null_char,trim(str)//c_null_char) !create new directory
      return
 end function setenv_int
 function setenv_float(env,floatval)
      implicit none
      integer :: setenv_float
      character(len=*) :: env
      real(wp) :: floatval
      character(len=30) :: str
      write(str,'(f14.6)')floatval
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
      logical :: bool
      character(len=*) :: fil
      character(len=*) :: str
      character(len=512) :: tmp
      integer :: io
      bool=.false.
      open(newunit=ich,file=fil)
      do
        read(ich,'(a)',iostat=io) tmp
        if(io < 0)exit
        if(index(tmp,str).ne.0)then
          bool=.true.
          close(ich)
          exit
        endif
      enddo
      close(ich)
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
      logical :: bool
      character(len=*) :: fil
      character(len=*) :: str
      character(len=512) :: tmp,dum
      real(wp) :: val
      integer :: io
      bool=.false.
      val=0.0d0
      open(newunit=ich,file=fil)
      do
        read(ich,'(a)',iostat=io) tmp
        if(io < 0)exit
        io=index(tmp,str,.true.)
        if(io.ne.0)then
          io=io+len(str)
          tmp=adjustl(tmp(io+1:))
          read(tmp,*) val
          bool=.true.
          close(ich)
          exit
        endif
      enddo
      close(ich)
      !minigrep = bool
      return
end subroutine grepval

!--------------------------------------------------------------------------------------
! Write a file with a single line (INT)
!--------------------------------------------------------------------------------------
subroutine wrshort_int(fname,var)
      implicit none
      character(len=*) :: fname
      integer :: var
      integer :: io
       open(newunit=ich,file=fname)
       write(ich,*)var
       close(ich)
end subroutine wrshort_int

!--------------------------------------------------------------------------------------
! Write a file with a single line (REAL)
!--------------------------------------------------------------------------------------
subroutine wrshort_real(fname,var)
      implicit none
      character(len=*) :: fname
      real(wp) :: var
      integer :: io
       open(newunit=ich,file=fname)
       write(ich,*)var
       close(ich)
end subroutine wrshort_real

!--------------------------------------------------------------------------------------
! Write a file with a single line (STRING)
!--------------------------------------------------------------------------------------
subroutine wrshort_string(fname,var)
      implicit none
      character(len=*) :: fname
      character(len=*) :: var
      integer :: io
      open(newunit=ich,file=fname)
      write(ich,'(a)')var
      close(ich)
end subroutine wrshort_string

!--------------------------------------------------------------------------------------
! reads the first line of a file (INT)
!--------------------------------------------------------------------------------------
subroutine rdshort_int(fname,var)
      implicit none
      character(len=*) :: fname
      integer :: var
      integer :: io
       open(newunit=ich,file=fname)
       read(ich,*,iostat=io)var
       if (io .ne. 0) var = 0
       close(ich)
end subroutine rdshort_int

!--------------------------------------------------------------------------------------
! reads the first line of a file (REAL)
!--------------------------------------------------------------------------------------
subroutine rdshort_real(fname,var)
      implicit none
      character(len=*) :: fname
      real(wp) :: var
      integer :: io
       open(newunit=ich,file=fname)
       read(ich,*,iostat=io)var
       if (io .ne. 0) var = 0.0d0
       close(ich)
end subroutine rdshort_real

!--------------------------------------------------------------------------------------
! reads the first line of a file (STRING)
!--------------------------------------------------------------------------------------
subroutine rdshort_string(fname,var)
      implicit none
      character(len=*) :: fname
      character(len=*) :: var
      integer :: io
       open(newunit=ich,file=fname)
       read(ich,'(a)',iostat=io)var
       if (io .ne. 0) var = ''
       close(ich)
end subroutine rdshort_string

!----------------------------------------------------------------------------
! returns entry x of line as string, or as empty if end-of-line is reached

subroutine glinex(line,x,string)
      implicit none
      integer,intent(in) :: x
      character(len=*),intent(in)  :: line
      character(len=*),intent(out) :: string
      character(len=80) :: stmp
      character(len=1)  :: digit
      integer :: i,k,cs
      logical :: track
 
      string=''
      stmp=''
      cs=0
      track=.false.
      k=len(trim(line))
      do i=1,k
       digit=line(i:i)
       if(digit.eq.' ' .or. digit.eq.char(9))then  !should exclude tabstops and blanks, 9 is ascii code for tab
          if(track)then
            if(cs==x)then
              string=trim(stmp)
              exit
            endif
          endif
          track=.false.
          stmp=''
          cycle
       else
        if(.not.track) cs=cs+1
        track=.true.
        stmp=trim(stmp)//trim(digit)
        if(i.eq.k.and.x.lt.0)then  !if x is lt. 0, the last element of the line is returned
          string=trim(stmp)
          exit
        endif
       endif
      enddo
      return
end subroutine glinex


!----------------------------------------------------------------------------
! write a wall potential in a file used as xtb input

subroutine write_wall(env,n1,rabc1,rabc12)
  use iso_fortran_env, only : wp => real64
  use crest_data

  implicit none

  type(systemdata)     :: env 
  integer, intent(in)  :: n1
  real(wp),intent(in)  :: rabc1(3),rabc12(3)
  character (len=8)    :: flag
  
  open(unit=31,file='xcontrol')
  flag='$'
  write(31,'(a,"wall")') trim(flag)
  write(31,'(3x,"potential=polynomial")')
  write(31,'(3x,"ellipsoid:",1x,3(g0,",",1x),"all")') rabc12
  write(31,'(3x,"ellipsoid:",1x,3(g0,",",1x),"1-",i0)') rabc1,n1
  call write_cts(31,env%cts)
  call write_cts_biasext(31,env%cts)
  if(env%cts%used .eq. .true.) then !Only, if user set constrians is an $end written
     write(31,'(a)') '$end'
  end if

  close(31)

end subroutine write_wall


!----------------------------------------------------------------------------
! cut the leading x words (seperated by blanks) of a string, and overwrite it

subroutine clinex(line,x)
      implicit none
      integer,intent(in) :: x
      character(len=*),intent(inout)  :: line
      character(len=1)  :: digit
      integer :: i,j,k,cs
      logical :: track

      cs=0
      track=.false.
      k=len(trim(line))
      do i=1,k
       digit=line(i:i)
       if(digit.eq.' ' .or. digit.eq.char(9))then  !should exclude tabstops and blanks, 9 is ascii code for tab
          track=.false.
          cycle
       else
        if(.not.track) cs=cs+1
        if(cs.gt.x)then
          j=i-1
          line=line(j:k)
          exit
        endif
        track=.true.
       endif
      enddo
      return
end subroutine clinex

!-------------------------------------------------------------
! convert a string to upper or lower case
!--------------------------------------------------------------
subroutine to_upper(str)
      character(*), intent(in out) :: str
      integer :: i

       do i = 1, len(str)
         select case(str(i:i))
           case("a":"z")
             str(i:i) = achar(iachar(str(i:i))-32)
         end select
       end do
end subroutine to_upper
subroutine to_lower(str)
      character(*), intent(in out) :: str
      integer :: i

      do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
      end do
end subroutine to_lower

!-----------------------------------------------------------------
! a simple routine to lock a feature of the code behind a password
!-----------------------------------------------------------------
subroutine pwdlock(pwdref)
      implicit none
      character(len=*) :: pwdref    ! check for this password
      character(len=:),allocatable :: pwd
      pwd=pwdref !allocate max length of pwd 
      write(*,'(/,1x,a)',advance='no') 'Locked feature. Enter password: '
      read(*,*) pwd
      if(pwd == pwdref)then
        write(*,'(1x,a,/)') 'Valid. Continue.'
      else
        error stop 'Invalid. Stop.'
      endif 
      return
end subroutine pwdlock

!-----------------------------------------------------------------------------------
! copy a coord file until an $set-block is encountered
!-----------------------------------------------------------------------------------
subroutine clear_setblock(fname)
      implicit none
      character(len=*) :: fname
      character(len=512) :: atmp
      integer :: i,j,k,iost
      integer :: ich,ich2

      open(newunit=ich,file=fname)
      open(newunit=ich2,file='.setdgtmp')

      do
           read(ich,'(a)',iostat=iost)atmp
           if(iost < 0)exit
           if((index(atmp,'$set').ne.0).or.  &
           &  (index(atmp,'$end').ne.0))then
             write(ich2,'(a)')'$end'
             exit
           else
             write(ich2,'(a)')trim(atmp)
           endif
      enddo
      close(ich,status='delete')
      close(ich2)
      call rename('.setdgtmp',fname)
end subroutine clear_setblock

!====================================================================!
function filechecker(fin,fout) result(have)
    implicit none
    logical :: have
    character(len=*),intent(in) :: fin
    character(len=*),intent(out) :: fout
    have = .false.
    inquire(file=fin,exist=have)
    if(have)then
        fout=fin
    else
        fout=''
    endif
    return
end function filechecker

! Checking directories with Fortran's inquire is not handled by the standard.
! The interpretation whether or not to report a directory as file is compiler
! specific and therefore always an extension to the Fortran standard.
function directory_exist(file) result(exist)
    character(len=*), intent(in) :: file
    logical :: exist
#ifdef __INTEL_COMPILER
    ! Intel provides the directory extension to inquire to handle this case
    inquire(directory=file, exist=exist)
#else
    ! GCC handles directories as files, to make sure we get a directory and
    ! not a file append a path separator and the current dir
    inquire(file=trim(file)//"/.", exist=exist)
#endif
end function directory_exist

end module iomod
