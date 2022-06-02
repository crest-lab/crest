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

!---------------------------------------------------------------------------------
! Routines for the handling of molecules in the chemoinformatical *.SDF format
!---------------------------------------------------------------------------------
subroutine inpsdf(env,fname)
      use iso_fortran_env, only: wp => real64
      use crest_data
      use strucrd
      implicit none
      !type(options) :: opt
      type(systemdata) :: env
      character(len=*) :: fname
      integer :: i
      env%sdfformat = .true.
      call checkcoordtype(fname,i)
      if(i==31)then  
        call rdsdf(env%sdf,fname)
      elseif( i==32) then
        call  rdsdf3000(env%sdf,fname)
      endif
      return
end subroutine inpsdf

!----------------------------------------------------------------------------------
subroutine rdsdf(sdf,fname)
      use iso_fortran_env, only: wp => real64
      use iomod
      use crest_data
      use strucrd
      use filemod
      implicit none
      type(sdfobj) :: sdf
      type(filetype) :: myfile
      character(len=*),intent(in) :: fname
      integer  :: nat
      integer :: i,j,k,l
      character(len=256) :: atmp,btmp

    !--- data and some allocation
      sdf%V3000 = .false.
      call rdnat(fname,nat)
      allocate(sdf%hblock(4))
      allocate(sdf%cblock(nat))
      sdf%nat = nat

    !--- read the file
      call myfile%open(fname)  

    !--- transfer header block  
      do i=1,4
         sdf%hblock(i) = ''
         sdf%hblock(i) = trim(myfile%line(i))
      enddo

    !--- transfer atom info block  
      do i=1,nat
        k =  i+4
        btmp = trim(myfile%line(k))
        atmp = getlarg(btmp,4)
        l = index(btmp,trim(atmp)) + 1
        j = len_trim(btmp)
        sdf%cblock(i) = trim(btmp(l:j))
      enddo

    !-- check number of lines for misc-block
      k = 4 + nat + 1   
      l = 0
      do i=k,myfile%nlines
        l = l + 1
        atmp = myfile%line(i)
        if(index(atmp,'M').ne.0 .and. &
        &  index(atmp,'END').ne.0)then
          exit
        endif 
      enddo  
      allocate(sdf%miscblock(l))
      sdf%nmisc = l
      j=0
      do i=k,l+k-1
        j = j+1
        sdf%miscblock(j) = trim(myfile%line(i))
      !  write(*,*) trim(sdf%miscblock(j))
      enddo
      call myfile%close()
     return
end subroutine rdsdf

subroutine rdsdf3000(sdf,fname)
      use iso_fortran_env, only: wp => real64
      use iomod
      use crest_data
      use strucrd
      use filemod
      implicit none
      type(sdfobj) :: sdf
      type(filetype) :: myfile
      character(len=*),intent(in) :: fname
      integer  :: nat
      integer :: i,j,k,l
      character(len=256) :: atmp,btmp

    !--- data and some allocation
      sdf%V3000 = .true.
      call rdnat(fname,nat)
      allocate(sdf%hblock(4))
      allocate(sdf%cblock(nat))
      sdf%nat = nat

    !--- read the file
      call myfile%open(fname)  

    !--- transfer header block  
      do i=1,4
         sdf%hblock(i) = trim(myfile%line(i))
      enddo

    !--- read counts line
    sdf%countsline  = trim(myfile%line(6))

    !--- transfer atom info block  
      do i=1,nat
        k =  i+7
        btmp = trim(myfile%line(k))
        call clinex(btmp,7)
        sdf%cblock(i) = trim(btmp)
      !  write(*,*) trim(btmp)
      enddo

    !-- check number of lines for misc-block
      k = 7 + nat + 2
      l = 0
      do i=k,myfile%nlines
        l = l + 1
        atmp = myfile%line(i)
        if(index(atmp,'M').ne.0 .and. &
        &  index(atmp,'END').ne.0 .and. &
        &  index(atmp,'V30').eq.0 )then
          exit
        endif 
      enddo  
      allocate(sdf%miscblock(l))
      sdf%nmisc = l
      j=0
      do i=k,l+k-1
        j = j+1
        sdf%miscblock(j) = trim(myfile%line(i))
      !  write(*,*) trim(sdf%miscblock(j))
      enddo
      call myfile%close()
     return
end subroutine rdsdf3000

subroutine wrsdfens(sdf,fname,oname)
      use iso_fortran_env, only: wp => real64
      use iomod
      use crest_data
      use strucrd, only: rdensembleparam,rdensemble,i2e
      implicit none
      type(sdfobj) :: sdf
      character(len=*),intent(in) :: fname
      character(len=*) :: oname
      integer  :: nat,nall
      integer,allocatable  :: at(:)
      real(wp),allocatable :: eread(:)
      real(wp),allocatable :: xyz(:,:,:)
      integer :: i,j,ich
    !---- read existing ensemble
      call rdensembleparam(fname,nat,nall)
      allocate(at(nat),eread(nat),xyz(3,nat,nall))
      call rdensemble(fname,nat,nall,at,xyz,eread)

    !---- write ensemble
      open(newunit=ich,file=trim(oname))
      do i=1,nall
         do j=1,4
           if(j.eq.1)then
              write(ich,'(a,5x,a,i0)') trim(sdf%hblock(j)),'#',i
           else
              write(ich,'(a)') trim(sdf%hblock(j))
           endif
         enddo
         if(.not.sdf%V3000)then  !--- V2000 format
           do j=1,nat
             write(ich,'(3f10.4,1x,a2,1x,a)')xyz(1:3,j,i),i2e(at(j),'nc'),trim(sdf%cblock(j))
           enddo
           do j=1,sdf%nmisc
             write(ich,'(a)') trim(sdf%miscblock(j))
           enddo
           write(ich,'(a)') '$$$$'
         else                   !--- V3000 format
           write(ich,'("M V30 BEGIN CTAB")')  
           write(ich,'(a)') trim(sdf%countsline)
           write(ich,'("M V30 BEGIN ATOM")')  
           do j=1,nat           
             write(ich,'(a,1x,i0,1x,a,3f10.4,1x,a)') 'M V30',j, &
             &     i2e(at(j),'nc'),xyz(1:3,j,i),trim(sdf%cblock(j))
           enddo
           write(ich,'("M V30 END ATOM")') 
           do j=1,sdf%nmisc
             write(ich,'(a)') trim(sdf%miscblock(j))
           enddo
           write(ich,'(a)') '$$$$'
         endif
      enddo

      deallocate(xyz,eread,at) 
      call sdf%deallocate()
      return
end subroutine wrsdfens

