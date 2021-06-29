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

!==========================================================================!
!   On request, get the bonds from the topology and write a constraints
!   file with a bond length constrain on each bond.
! 
!   A file "bondlengths" is written with a constraint on each bond.
!   The input argument is the required force constant.
!==========================================================================!
subroutine autoBondConstraint(filename,forceconstant,wbofile)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     character(len=*) :: filename
     character(len=*) :: wbofile
     real(wp) :: forceconstant
     type(zmolecule) :: zmol    
   !--- get topology
     call simpletopo_file(filename,zmol,.false.,.false.,wbofile)
   !--- get bond matrix and write "bondlengths" file
     call getbmat(zmol,1,forceconstant)
     call zmol%deallocate()
     return
end subroutine autoBondConstraint

subroutine autoBondConstraint_withEZ(filename,forceconstant,wbofile)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     character(len=*) :: filename
     character(len=*) :: wbofile
     real(wp) :: forceconstant
     type(zmolecule) :: zmol    
   !--- get topology
     call simpletopo_file(filename,zmol,.false.,.false.,wbofile)
   !--- get bond matrix and write "bondlengths" file
     call getbmat(zmol,5,forceconstant)
     call zmol%deallocate()
     return
end subroutine autoBondConstraint_withEZ

!==========================================================================!
!   On request, get the bonds from the topology and write a constraints
!   file with a bond length constrain on each bond, but only if it is
!   an M-X bond (M = transition metal atom)
! 
!   A file "bondlengths" is written with a constraint on each bond.
!   The input argument is the required force constant.
!==========================================================================!
subroutine autoMetalConstraint(filename,forceconstant,wbofile)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     character(len=*) :: filename
     character(len=*) :: wbofile
     real(wp) :: forceconstant
     type(zmolecule) :: zmol
   !--- get topology
     call simpletopo_file(filename,zmol,.false.,.false.,wbofile)
   !--- get bond matrix and write "bondlengths" file
     call getbmat(zmol,2,forceconstant)
     call zmol%deallocate()
     return
end subroutine autoMetalConstraint

!==========================================================================!
!   On request, get the bonds from the topology and write a constraints
!   file with a bond length constrain on each bond, except if the bond
!   is an X-H bond
!
!   A file "bondlengths" is written with a constraint on each bond.
!   The input argument is the required force constant.
!==========================================================================!
subroutine autoHeavyConstraint(filename,forceconstant,wbofile)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     character(len=*) :: filename
     character(len=*) :: wbofile
     real(wp) :: forceconstant
     type(zmolecule) :: zmol
   !--- get topology
     call simpletopo_file(filename,zmol,.false.,.false.,'')
   !--- get bond matrix and write "bondlengths" file
     call getbmat(zmol,3,forceconstant)
     call zmol%deallocate()
     return
end subroutine autoHeavyConstraint


!==========================================================================!
!   On request, get the bonds from the topology and write a constraints
!   file with a bond length constrain on each bond, but only if the bond
!   is an X-H bond
!
!   A file "bondlengths" is written with a constraint on each bond.
!   The input argument is the required force constant.
!==========================================================================!
subroutine autoHydrogenConstraint(filename,forceconstant,wbofile)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     character(len=*) :: filename
     character(len=*) :: wbofile
     real(wp) :: forceconstant
     type(zmolecule) :: zmol
   !--- get topology
     call simpletopo_file(filename,zmol,.false.,.false.,'')
   !--- get bond matrix and write "bondlengths" file
     call getbmat(zmol,4,forceconstant)
     call zmol%deallocate()
     return
end subroutine autoHydrogenConstraint

!===================================================================================!
!  Some routines related to the Bond matrix (BMAT)
!===================================================================================!
subroutine getbmat(zmol,r,force)
     use iso_fortran_env, wp => real64
     use zdata
     implicit none
     type(zmolecule) :: zmol
     integer :: r
     integer :: nat
     real(wp) :: force
     real(wp),allocatable :: xyz(:,:)
     real(wp),allocatable :: bmat(:,:)
     integer,allocatable :: at(:)
     integer,allocatable :: bonds(:,:)
     integer :: i,j,k,l,nb

     integer,parameter :: cbonds = 1
     integer,parameter :: cmetal = 2
     integer,parameter :: cheavy = 3
     integer,parameter :: chydro = 4
     integer,parameter :: cistrans = 5
     integer,parameter :: countbonds = 6

     real(wp),parameter :: bohr = 0.52917726_wp

     nat = zmol%nat
     at = zmol%at
     allocate(xyz(3,nat),bmat(nat,nat), source = 0.0_wp)
    
     do i=1,nat
        xyz(:,i)=zmol%zat(i)%cart(:)
     enddo

     nb = 0
     do i=1,nat
        do j=1,nat
          if(any(zmol%zat(i)%ngh(:) .eq. j))then  !only include bonds from the neighbour lists
            bmat(i,j) = zmol%distmat(i,j)
            bmat(i,j) = bmat(i,j)*bohr   !BMAT is in Angstroem!
            nb = nb +1
          else
           cycle
          endif
        enddo
     enddo

     !--- write the constrain file
     if(r == cbonds)then
     call writeBmatconstr(nat,bmat,force)
     endif
     if(r == cmetal)then
      call writeMetalconstr(nat,at,bmat,force)
     endif
     if(r == cheavy)then
      call writeHeavyconstr(nat,at,bmat,force)
     endif
     if(r == chydro)then
      call writeHydrogenconstr(nat,at,bmat,force)
     endif
     if(r == cistrans)then
      call writeBmatconstr_withEZ(zmol,nat,bmat,force)
     endif

     !--- utility
     if(r == countbonds)then
       allocate(bonds(2,nb), source=0)  
       nb=0  
       do i=1,nat
         do j=1,i
           if(bmat(i,j) > 1d-6)then
            nb=nb+1
            bonds(1,nb) = j
            bonds(2,nb) = i
           endif    
         enddo
       enddo
      deallocate(bonds)
     endif


     deallocate(bmat,xyz,at)
     return
end subroutine getbmat

!====================================================!
! Write a constraint file with all bonds
!====================================================!
subroutine writeBmatconstr(nat,bmat,force)
     use iso_fortran_env, wp => real64
     implicit none
     integer :: nat
     real(wp) :: bmat(nat,nat)
     integer :: i,j,k,l,ich
     real(wp) :: force
     character(len=20) :: dumm


     open(newunit=ich,file='bondlengths')
     write(ich,'(a)') '$constrain'
     write(dumm,'(f16.4)') force
     write(ich,'(3x,a,a)')'force constant=',adjustl(trim(dumm))
     do i=1,nat
        do j=i+1,nat
           if(bmat(i,j).gt.0.1_wp)then
             write(ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',i,',',j,',',bmat(i,j)
           endif
        enddo
     enddo
     write(ich,'(a)') '$end'
     close(ich)
     return
end subroutine writeBmatconstr

!====================================================!
! Write a constraint file with all TM bonds
!====================================================!
subroutine writeMetalconstr(nat,at,bmat,force)
     use iso_fortran_env, wp => real64
     implicit none
     integer :: nat
     integer :: at(nat)
     real(wp) :: bmat(nat,nat)
     integer :: i,j,k,l,ich
     real(wp) :: force
     character(len=20) :: dumm
     logical :: isTMetal  !this is a function

     open(newunit=ich,file='bondlengths')
     write(ich,'(a)') '$constrain'
     write(dumm,'(f16.4)') force
     write(ich,'(3x,a,a)')'force constant=',adjustl(trim(dumm))
     do i=1,nat
        do j=i+1,nat
           if(bmat(i,j).gt.0.1_wp)then
             if(isTMetal(at(i)) .or. isTMetal(at(j)))then
             write(ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',i,',',j,',',bmat(i,j)
             endif
           endif
        enddo
     enddo
     write(ich,'(a)') '$end'
     close(ich)
     return
end subroutine writeMetalconstr

logical function isTMetal(i)
     implicit none
     integer :: i
     isTMetal = .false.
     if(i.ge.21 .and. i.le.30) isTMetal = .true.
     if(i.ge.39 .and. i.le.48) isTMetal = .true.
     if(i.ge.57 .and. i.le.80) isTMetal = .true.
     return
end function isTMetal


!====================================================!
! Write a constraint file with all bonds except X-H
!====================================================!
subroutine writeHeavyconstr(nat,at,bmat,force)
     use iso_fortran_env, wp => real64
     implicit none
     integer :: nat
     integer :: at(nat)
     real(wp) :: bmat(nat,nat)
     integer :: i,j,k,l,ich
     real(wp) :: force
     character(len=20) :: dumm
     logical :: isTMetal  !this is a function

     open(newunit=ich,file='bondlengths')
     write(ich,'(a)') '$constrain'
     write(dumm,'(f16.4)') force
     write(ich,'(3x,a,a)')'force constant=',adjustl(trim(dumm))
     do i=1,nat
        do j=i+1,nat
           if(bmat(i,j).gt.0.1_wp)then
             if((at(i).ne.1) .and. (at(j).ne.1))then
             write(ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',i,',',j,',',bmat(i,j)
             endif
           endif
        enddo
     enddo
     write(ich,'(a)') '$end'
     close(ich)
     return
end subroutine writeHeavyconstr
!====================================================!
! Write a constraint file with all X-H bonds
!====================================================!
subroutine writeHydrogenconstr(nat,at,bmat,force)
     use iso_fortran_env, wp => real64
     implicit none
     integer :: nat
     integer :: at(nat)
     real(wp) :: bmat(nat,nat)
     integer :: i,j,k,l,ich
     real(wp) :: force
     character(len=20) :: dumm
     logical :: isTMetal  !this is a function

     open(newunit=ich,file='bondlengths')
     write(ich,'(a)') '$constrain'
     write(dumm,'(f16.4)') force
     write(ich,'(3x,a,a)')'force constant=',adjustl(trim(dumm))
     do i=1,nat
        do j=i+1,nat
           if(bmat(i,j).gt.0.1_wp)then
             if((at(i).eq.1) .or. (at(j).eq.1))then
             write(ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',i,',',j,',',bmat(i,j)
             endif
           endif
        enddo
     enddo
     write(ich,'(a)') '$end'
     close(ich)
     return
end subroutine writeHydrogenconstr


!====================================================!
! Write a constraint file with all bonds
!====================================================!
subroutine writeBmatconstr_withEZ(zmol,nat,bmat,force)
     use iso_fortran_env, wp => real64
     use crest_data, only : bohr
     use zdata
     implicit none
     type(zmolecule) :: zmol
     integer :: nat
     real(wp) :: bmat(nat,nat)
     integer :: i,j,k,l,ich
     real(wp) :: force
     character(len=20) :: dumm

     integer :: ni,nj
     real(wp) :: ndist
     integer :: r
     logical :: smallring

     real(wp),parameter :: springexp = 4.0_wp
     real(wp),parameter :: cclen = 1.384_wp
     integer,parameter  :: ringmax = 5
     integer,parameter  :: ringmax2 = 8

     open(newunit=ich,file='bondlengths')
     write(ich,'(a)') '$constrain'
     write(dumm,'(f16.4)') force
     write(ich,'(3x,a,a)')'force constant=',adjustl(trim(dumm))
     do i=1,nat
     do j=i+1,nat
      if(bmat(i,j).gt.0.1_wp)then

          smallring = .false.
          !if(zmol%nri > 0)then
          !  do r=1,zmol%nri
          !    if(any(zmol%zri(r)%rlist==i) .and. &
          !    &  any(zmol%zri(r)%rlist==j) .and. &
          !    &  (zmol%zri(r)%rs <= ringmax))then
          !    smallring = .true.
          !    endif
          !  enddo
          !endif

          if(.not.smallring)then
             write(ich,'(3x,a,i0,a,i0,a,f8.5)') 'distance: ',i,', ',j,', ',bmat(i,j)
          else
             write(ich,'(3x,a,i0,a,i0,a,f8.5,a,f8.5)') 'DISTANCE: ',i,', ',j,', ',bmat(i,j), &
            &       ', ',springexp 
          endif
      endif
     enddo
     enddo
     !write(ich,*) '#E/Z constraints'
     ILOOP : do i=1,nat
     JLOOP : do j=i+1,nat
      if(bmat(i,j).gt.0.1_wp)then
       if((zmol%at(i)==6 .and. zmol%at(j)==6).and.  &                     ! ij are Carbon atoms?
       &  (nint(zmol%zat(i)%cn)==3 .and. nint(zmol%zat(j)%cn)==3).and. &  ! ij are sp2 C?
       &  (bmat(i,j)<cclen)   )then                                       ! r_ij could be C=C bond?

          ! check if ij define a C=C bond within a small ring? (cycle, if so)
          do r=1,zmol%nri
            if(any(zmol%zri(r)%rlist==i) .and. &
            &  any(zmol%zri(r)%rlist==j) .and. &
            &  (zmol%zri(r)%rs < ringmax2))then
            cycle JLOOP
            endif
          enddo
          
          do k=1,zmol%zat(i)%nei
            ni = zmol%zat(i)%ngh(k)
            if(ni == j) cycle
          do l=1,zmol%zat(j)%nei
            nj = zmol%zat(j)%ngh(l)
            if( nj == i ) cycle
            ndist = zmol%dist(ni,nj)*bohr
            !write(ich,'(3x,a,i0,a,i0,a,f8.5)') 'distance: ',ni,', ',nj,', ',ndist
            !--- WARNING: The constraint has to be case-sensitive in order to have an adjustable springexp.
            write(ich,'(3x,a,i0,a,i0,a,f8.5,a,f8.5)') 'DISTANCE: ',ni,', ',nj,', ',ndist, &
            &       ', ',springexp
          enddo
          enddo
       endif
      endif
     enddo JLOOP
     enddo ILOOP
     write(ich,'(a)') '$end'
     close(ich)

     return
end subroutine writeBmatconstr_withEZ
