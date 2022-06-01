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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine getname1(i,atmp)
      integer :: i
      character(len=*) :: atmp

      write(atmp,'(''scoord.'',i0)')i

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine checkname_xyz(base,fname,checkname)
      character(len=*) :: base,fname,checkname
      integer :: i
      logical :: ex

      i=0
   
      do
        !if(i.lt.1)then
        !  write(checkname,'(a,''.xyz'')')trim(base)
        !else
          write(checkname,'(a,''_'',i0,''.xyz'')')trim(base),i
        !endif
        inquire(file=trim(checkname),exist=ex)
        if(ex)then
          i=i+1
        else
         exit
        endif
      enddo

      j=max(0,i-1)
      
      !if(j.lt.1)then
      !   write(fname,'(a,''.xyz'')')trim(base)
      !else
         write(fname,'(a,''_'',i0,''.xyz'')')trim(base),j
      !endif

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine checkname_tmp(base,fname,checkname)
      character(len=*) :: base,fname,checkname
      integer :: i
      logical :: ex

      i=0

      do
        !if(i.lt.1)then
        !  write(checkname,'(a,''.xyz'')')trim(base)
        !else
          write(checkname,'(a,''_'',i0,''.tmp'')')trim(base),i
        !endif
        inquire(file=trim(checkname),exist=ex)
        if(ex)then
          i=i+1
        else
         exit
        endif
      enddo

      j=max(0,i-1)

      !if(j.lt.1)then
      !   write(fname,'(a,''.xyz'')')trim(base)
      !else
         write(fname,'(a,''_'',i0,''.tmp'')')trim(base),j
      !endif

      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine getname_dir(base,fname)
      use iomod, only : directory_exist
      character(len=*) :: base,fname
      integer :: i
      logical :: ex
      character(len=512)  :: dir

      i=0
      do
        i=i+1
        write(dir,'(a,i0)')trim(base),i
        ! Inquire is not consitent across compilers
        ex = directory_exist(dir)
        if(.not.ex)then
          fname=trim(dir) 
          exit
        endif
      enddo

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! get cregen warning
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine getcregenwarning(fname)
       implicit none
       character(len=*) :: fname
       character(len=512) :: str
       integer :: io

       open(unit=654,file=fname)
       do
          read(654,'(a)',iostat=io)str
          if(io<0)exit
          if(index(str,'WARNING:').ne.0)then
            write(*,'(a)')trim(str)
            read(654,'(a)',iostat=io)str
            write(*,'(a)')trim(str)
            read(654,'(a)',iostat=io)str
            write(*,'(a)')trim(str)
            exit
          endif
       enddo
       close(654)

       end subroutine getcregenwarning 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine delsubstr2(str,substr)

! Deletes first occurrence of substring 'substr' from string 'str' and
! shifts characters left to fill hole. Trailing spaces or blanks are
! not considered part of 'substr'.
      implicit none
      character(len=*):: str,substr
      integer :: lensubstr
      integer :: ipos

      lensubstr=len_trim(substr)
      ipos=index(str,substr)
      if(ipos==0) return
      if(ipos == 1) then
         str=str(lensubstr+1:)
      else
         str=str(:ipos-1)//str(ipos+lensubstr:)
      end if
      return

      end subroutine delsubstr2

!----------------------------------------------------------------------------
! modified readline routine

      subroutine readline3(line,floats)
      implicit none
      real*8,intent(out) :: floats(10)
      character(len=*) :: line
      character(len=80) :: strings(3)

      real*8 :: num
      character(len=80) :: stmp,str
      character(len=1)  :: digit
      integer :: i,ty,cs,cf

      stmp=''
      cs=1
      cf=1
      strings(:)=''
      floats=0.0d0
      do i=1,len(trim(line))
       digit=line(i:i)
       if(digit.ne.' '.and.digit.ne.char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
        stmp=trim(stmp)//trim(digit)
       elseif(stmp.ne.'')then
        call checktype(stmp,num,str,ty)      !get type of string, 0=number, 1=character
        if(ty.eq.0) then
         floats(cf)=num
         cf=cf+1
        elseif(ty.eq.1) then
         strings(cs)=str
         cs=cs+1
        else
          write(*,*)'Problem in checktype, must abort'
          exit
        endif
        stmp=''
       endif
       if(i.eq.len(trim(line))) then  !special case: end of line
        call checktype(stmp,num,str,ty)
        if(ty.eq.0) then
         floats(cf)=num
         cf=cf+1
        elseif(ty.eq.1) then
         strings(cs)=str
         cs=cs+1
        else
          write(*,*)'Problem in checktype, must abort'
          exit
        endif
        stmp=''
       endif
      enddo
      cs=cs-1
      cf=cf-1
      end subroutine readline3

!-------------------------------------------------------------------------
! append content of test file "from" into text file "to", similar to "cat A >> B"
! but specifically for XYZ files
!-------------------------------------------------------------------------
      subroutine XYZappendto(from,to)
      implicit none
      integer :: io
      character(len=*) :: from
      character(len=*) :: to
      integer :: i,nat
      character(len=1024) :: str
      open(unit=666,file=to)
      open(unit=777,file=from)
      do
        read(666,*,iostat=io)
        if(io<0)then
          backspace(666)  
          exit
        endif
      enddo
      read(777,*)nat
      write(666,*)nat
      do i=1,nat+1
        read(777,'(a)',iostat=io)str
        if(io<0)exit
        write(666,'(a)')trim(str)
      enddo
      close(777)
      close(666)
      end subroutine XYZappendto

!-------------------------------------------------------------------------
! append content of test file "from" into text file "to", similar to "cat A >> B"
! but specifically for TRJ files, but leaves out the first structure on "from"
!-------------------------------------------------------------------------
      subroutine TRJappendto_skipfirst(from,to)
      implicit none
      integer :: io,ich,och
      integer :: nat,i
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
      !---first structure is read, but not copied 
      read(och,'(a)',iostat=io)str
      if(io==0)then
        read(str,*)nat
        read(och,'(a)',iostat=io)str
        do i=1,nat
          read(och,'(a)',iostat=io)str
        enddo
      endif
      !---------------------------
      do
        read(och,'(a)',iostat=io)str
        if(io<0)exit
        write(ich,'(a)')trim(str)
      enddo
      close(och)
      close(ich)
      end subroutine TRJappendto_skipfirst


!-------------------------------------------------------------------------
! append content of test file "from" into text file "to", similar to "cat A >> B"
! but specifically for XYZ files
!-------------------------------------------------------------------------
      subroutine XYZappend2(from,to)
      implicit none
      integer :: io,to
      character(len=*) :: from
      integer :: fromch
      integer :: i,nat
      character(len=1024) :: str

      open(newunit=fromch,file=from)

      do
        read(to,*,iostat=io)
        if(io<0)exit
      enddo

      read(fromch,*)nat
      write(to,*)nat

      do i=1,nat+1
        read(fromch,'(a)',iostat=io)str
        if(io<0)exit
        write(to,'(a)')trim(str)
      enddo
      close(fromch)

      end subroutine XYZappend2

 
!-------------------------------------------------------------------------
! count number of arguments (seperated by space or tab) of a line
!-------------------------------------------------------------------------
      function nargums(line)
      implicit none
      integer :: nargums
      character(len=*),intent(in) :: line
      character(len=256) :: line2
      character(len=1) :: digit
      integer :: i,j
      logical :: newarg
      newarg=.true.
      nargums=0
      line2=adjustl(line)
      do i=1,len(trim(line2))
         digit=line2(i:i)
         j=i+1
         if((digit==' ').or.(digit==char(8)))then
            line2=line2(j:)
            newarg=.true.
         elseif(newarg)then
            nargums=nargums+1
            newarg=.false.
         endif
      enddo
      end function nargums

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tmp_anmrrc(atmp)
      use iomod
      implicit none
      character(len=*) :: atmp
      character(len=80) :: fname
      character(len=512) :: fname2
      logical :: ex,ex2
      !use a global .anmrrc if present (commented out for the use with orca)
      inquire(file='~/.anmrrc',exist=ex)
      if(ex)then
        atmp='~/.anmrrc'
        return
      else
      !if not present try to use a custom .anmrrc in the current directory
        fname2='.anmrrc'
        inquire(file=fname2,exist=ex2)
        if(ex2) then
           atmp=trim(fname2)
           return
        else
          !else, write a temporary .anmrrc
          fname='.anmrrc.tmp'
          !atmp='./.anmrrc.tmp'
          atmp='.anmrrc.tmp'
          !if the file exists delete it
          call remove(fname)
        open(unit=1234, file=fname)
       write(1234,'(''7 8   # XH acid atoms '')')
       write(1234,'(''elem.  calc.          exp. shift      active '')')
       write(1234,'(''1    31.538               0.0           1'')')
       write(1234,'(''6    187.34               0.0           0'')')
       write(1234,'(''9    320                  0.0           1'')')
       write(1234,'(''15   195                  0.0           1'')')
        close(1243)
         return
         endif
      endif
      end

      subroutine getanmrrc(atmp,fail)
      implicit none
      character(len=*),intent(inout) :: atmp
      logical,intent(out) :: fail
      logical :: ex
      inquire(file=atmp,exist=ex)
      if(ex)atmp=trim(atmp)
      fail=.not.ex
      return
      end subroutine getanmrrc

!-------------------------------------------------------------------------
!     read a string "str" and get the value "val" (as a string also)
!     that was assignet to the argument "arg"
!-------------------------------------------------------------------------
      subroutine rdarg(str,arg,val)
      implicit none
      character(len=*) :: str
      character(len=*) :: arg
      character(len=*) :: val
      character(len=512) :: tmp
      integer :: io
      val=''
      io=index(str,arg,.true.)
      if(io.ne.0)then
          io=io+len(arg)
          tmp=str(io:)
          val=trim(tmp)
      else
         val=''
      endif
      return
      end subroutine rdarg

!-------------------------------------------------------------------------
! minigrep: a grep subroutine that returns true or false, depending on
! if the substring "str" is present in the file "fil"
!-------------------------------------------------------------------------
      subroutine gettime(fil,secs)
      implicit none
      real*8 :: secs
      real*8 :: floats(10)
      character(len=*) :: fil
      character(len=512) :: tmp,tmp2
      integer :: io,ich
      secs=0.0d0
      open(newunit=ich,file=fil)
      do
        read(ich,'(a)',iostat=io) tmp
        if(io < 0)exit
        if(index(tmp,'finished run on').ne.0)then
          read(ich,'(a)') tmp
          read(ich,'(a)') tmp
          read(ich,'(a)') tmp
          read(ich,'(a)') tmp
          call rdarg(tmp,'time:',tmp2)
          call readline3(tmp2,floats)
        endif
      enddo
      secs=secs+floats(2)*86400.0d0  !days to seconds      
      secs=secs+floats(3)*3600.0d0   !hours to seconds
      secs=secs+floats(4)*60.0d0     !minutes to seconds
      secs=secs+floats(5)           
      return
      end subroutine gettime
