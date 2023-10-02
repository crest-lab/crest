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

subroutine rdcontrol(fname,oname,method,freqmode)
     use iso_fortran_env, wp => real64, stderr=>error_unit
     use filemod
     use iomod
     use atmasses
     use strucrd, only: i2e,e2i
     implicit none
     character(len=*) :: fname
     character(len=*) :: oname
     character(len=*) :: method
     type(filetype) :: ctrl
     integer :: i,j,k,l,m
     integer :: block(2)
     integer :: ml
     logical :: ex
     character(len=:),allocatable :: dum
     character(len=:),allocatable :: dum2 
     character(len=80) :: mline
     character(len=2) :: element(107)
     real(wp),allocatable :: amv(:)  !atom masses
     logical,allocatable :: ause(:)
     integer,allocatable :: subblock(:)
     character(len=:),allocatable :: hess_set
     integer,allocatable :: greps(:)
     logical :: freqmode

     !freqmode=.false.

     inquire(file=fname,exist=ex)

     if(.not.ex)then
        write(stderr,'(1x,3a)')'Input file ',trim(fname),' does not exist!'
        error stop 
     endif

  !--- get array with (lowercase) element symbolds
     do i=1,107
       element(i)=trim(i2e(i,'lc'))
     enddo
  !--- get atom masses
     allocate(amv(107),ause(107))  
     amv=ams(1:107)
     call freqmass(amv,method)  !get scaled masses
     ause=.true. !use all new masses since un-fitted elements are scaled anyways

  !--- open ctrl-file into memory  
     call ctrl%open(fname)

  !--- loock for the "atoms block
     call getblck(ctrl,'$atoms',block)

  !--- parse $atoms block for atom sub-blocks
     !first get all postions of the atom subblocks
     l=0
     do i=block(1),block(2)
       dum=ctrl%line(i)
       dum2=getlarg(dum,1)  
       call to_lower(dum2)
       if(any(element==trim(dum2)))then
         l=l+1  
       endif
     enddo
     allocate(subblock(l+1))
     subblock(l+1)=block(2)
     l=1
     do i=block(1),block(2)
       dum=ctrl%line(i)
       dum2=getlarg(dum,1)
       call to_lower(dum2)
       if(any(element==trim(dum2)))then
          subblock(l)=i
          l=l+1
       endif
     enddo
     !then check for each atom subblock is new masses must be written
     do j=1,l-1
        dum=ctrl%line(subblock(j))
        dum2=getlarg(dum,1)
        call to_lower(dum2)
        if(any(element==trim(dum2)))then
          k = e2i(dum2)
          if(ause(k))then !only do this if a new mass was read in
            write(mline,'(a,f12.8)')'   mass  =',amv(k)
            ml=len_trim(mline)
            ex=.false.
            do i=subblock(j),subblock(j+1)
               dum=ctrl%line(i)
               if(index(dum,'mass  =').ne.0) then
                dum(1:ml)=trim(mline)
                call ctrl%replace(i,dum)  !if "mass" is already present replace it in line
                ex=.true.
               endif
            enddo
            if(.not.ex)then !if it is not present we have to insert a new line, which requires bookkeeping
               m=subblock(j+1)-1
               dum=ctrl%line(m)
               if(len(dum).lt.81)then
                 i=max(0,(79-len_trim(dum)))
                 dum=trim(dum)//repeat(' ',i)//'\'
               else
                 dum(80:80)='\'
               endif
               call ctrl%replace(m,dum)
               call ctrl%insert(mline,subblock(j+1))
               do i=j+1,l  !update line numbers
                  subblock(i)=subblock(i)+1
               enddo
            endif
          endif
        endif
     enddo

!---- find the setting for calculating the non mass weighted hessian in the control file
     if(freqmode)then
     hess_set='$drvopt'
     call getblck(ctrl,hess_set,block)
     if(block(2).eq.0)then !no occurence, i.e., the keyword was not set in the coord yet
       call ctrl%replace(ctrl%nlines,hess_set)  !replace last line of control file (should be $end) by new keyword
       call ctrl%write('   frequency analysis only')
       call ctrl%write('$end')
     else
       greps = ctrl%findall('frequency analysis only',k)
       if(k.eq.0)then
       call ctrl%insert('   frequency analysis only',block(2))
       else
       call ctrl%replace(greps(k),'   frequency analysis only')
       endif
     endif
     endif

     ctrl%filename=trim(oname)
     call ctrl%flush
     call ctrl%close

     deallocate(ause,amv)
     return
end subroutine rdcontrol


!---- get the lines in which a block starts and ends.
subroutine getblck(ctrl,bname,b)
     use iso_fortran_env, wp => real64
     use filemod
     use iomod 
     implicit none
     type(filetype) :: ctrl
     character(len=*) :: bname
     integer :: b(2)
     integer :: i,n
     character(len=:),allocatable :: dum
     character(len=:),allocatable :: dum2

     b=0

     n=ctrl%nlines
     do i=1,n
        dum=ctrl%line(i)
        dum2=getlarg(dum,1)
        if(index(dum2,bname).ne.0)then
          b(1)=i
        endif
        if(index(dum2,'$').ne.0 .and. b(1).gt.0 &
        & .and. b(1).lt.i )then
          b(2)=i
          exit
        endif
     enddo
     return
end subroutine getblck


! get standard atmasses and fittet atmasses if they shall be uses
subroutine getatmss(mass,amv,ause)
     use iso_fortran_env, wp => real64
     use iomod
     use filemod
     use atmasses
     use strucrd, only: i2e
     implicit none
     character(len=*) :: mass
     real(wp) :: amv(107)
     logical :: ause(107)
     logical :: ex
     type(filetype) :: mfile
     integer :: i,j,k
     character(len=:),allocatable :: dum
     character(len=:),allocatable :: at
     amv=ams(1:107) !from module atmasses
     ause=.false.
     
     inquire(file=mass,exist=ex)
     if(ex)then
       write(*,'(1x,a,a)') 'File read for new atom masses:  ',trim(mass)
       call mfile%open(mass)
       ! new masses have to be read here 
        do i=1,mfile%nlines
          dum=mfile%line(i)
          if(index(dum,'element mass:').ne.0)then
            at=getlarg(dum,3)
            k=len_trim(at)
            at=at(1:k)//''
            read(at,*) j
            ause(j)=.true.
            at=getlarg(dum,4)
            read(at,*) amv(j)
            write(*,'(3x,a,a2,a,1x,f12.8)') 'element: ',i2e(j),'  mass:',amv(j)
          endif
        enddo
       call mfile%close
     endif

     return
end subroutine

!=============================================================================================
!---- add different masses to hessian calculation
!
! ON INPUT: fname - filename of the file that shall be appended
!           opt   - metadata object containing settings
!
subroutine add_mass_xtb(fname)
     use iso_fortran_env, wp => real64
     use filemod
     use atmasses
     use strucrd, only: rdnat,rdcoord
     implicit none
     character(len=*) :: fname
     type(filetype) :: coord
     real(wp),allocatable :: amv(:)
     logical,allocatable :: ause(:)
     integer :: n,i
     integer,allocatable :: at(:)
     real(wp),allocatable :: xyz(:,:)
     character(len=64) :: dum
  
     allocate(amv(107),ause(107))
     amv=ams !from module atmasses
     ause=.false.
     call freqmass(amv,'gfn2')  !get scaled masses

     call rdnat(fname,n)
     allocate(xyz(3,n),at(n))
     call rdcoord(fname,n,at,xyz)    
     do i=1,n
        ause(at(i))=.true.
     enddo
     deallocate(at,xyz)   

     call coord%open(fname)
     call coord%write('$hess')
     do i=1,107
        if(ause(i))then
          write(dum,'(a,i4,a,f10.4)')' element mass:',i,',',amv(i)
          call coord%write(trim(dum))
        endif
     enddo
     !call coord%write(' element mass: 1,    1.034')
     !call coord%write(' element mass: 6,    11.63')
     !call coord%write(' element mass: 7,    14.37')
     !call coord%write(' element mass: 8,    18.05')
     !call coord%write(' element mass: 9,    16.86')
     !call coord%write(' element mass: 16,   32.18')
     !call coord%write(' element mass: 17,   26.71')
     !call coord%write(' element mass: 35,   80.00')
     call coord%write('$end')
     call coord%flushclose

     deallocate(ause,amv)
end subroutine add_mass_xtb


