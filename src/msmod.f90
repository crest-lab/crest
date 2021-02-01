!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht
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

module msmod

      use iso_fortran_env, wp => real64
      implicit none

      
      !-- storage for a single mol
      type msmol
          integer :: nat
          integer,allocatable :: at(:)
          real(wp),allocatable :: xyz(:,:)
          real(wp) :: etot
          integer :: chrg
          integer :: nfrag
          integer,allocatable :: frag(:)
          integer :: gen
          contains
              procedure :: newmol => msmod_newmol
              procedure :: dist => msmoldist
      end type msmol

      public :: msmol

      !-- isomer list
      type ilist
          integer :: nmol
          type(msmol),allocatable :: mol(:)
          logical,allocatable :: new(:)
          contains
              procedure :: dealloc => deallocate_ilist
              procedure :: append => ilist_append
      end type ilist

      !-- fragmentized structure list
      type flist
          integer :: nmol
          type(msmol),allocatable :: mol(:)
          logical,allocatable :: new(:)
          contains
              procedure :: dealloc => deallocate_flist
      end type flist


      !-- global data for msreact
      type msobj

          real(wp) :: ethr    = 1.0_wp  !cov. struc. chekc. thr.
          real(wp) :: wbormsd = 0.5_wp  !wbo rmsd comp. thr.
         
          real(wp) :: T  = 3000.0_wp    !start temperature
          real(wp) :: fc = 0.05_wp      !start fc
          real(wp) :: cdist = 1.5_wp    !constraing distance scaling factor of rcov
          integer  :: maxc = 30         !max optimization cycle

          type(ilist) :: il   
          type(flist) :: fl
          
      end type msobj

      public :: msobj

      private
contains

subroutine deallocate_ilist(self)
     class(ilist) :: self
     if(allocated(self%mol)) deallocate(self%mol)
     if(allocated(self%new)) deallocate(self%new)
     return
end subroutine deallocate_ilist

subroutine deallocate_flist(self)
     class(flist) :: self
     if(allocated(self%mol)) deallocate(self%mol)
     if(allocated(self%new)) deallocate(self%new)
     return
end subroutine deallocate_flist

subroutine msmod_newmol(self,nat,xyz,at,etot,chrg,gen)
     implicit none
     class(msmol) :: self
     integer :: nat
     real(wp) :: xyz(3,nat)
     integer :: at(nat)
     real(wp) :: etot
     integer :: chrg
     integer :: gen
     self%nat = nat
     allocate(self%xyz(3,nat))
     self%xyz = xyz
     allocate(self%at(nat))
     self%at = at
     self%etot = etot
     self%chrg = chrg
     self%gen = gen
     return    
end subroutine msmod_newmol

!calculate the distance between two atoms i and j
function msmoldist(self,i,j) result(dist)
     implicit none
     class(msmol) :: self
     integer :: i,j
     real(wp) :: dist
     dist = 0.0_wp
     if(allocated(self%xyz))then
      dist=(self%xyz(1,i)-self%xyz(1,j))**2 + &
     &     (self%xyz(2,i)-self%xyz(2,j))**2 + &
     &     (self%xyz(3,i)-self%xyz(3,j))**2
      dist = sqrt(dist)
     endif
     return
end function msmoldist



subroutine ilist_append(self,nat,at,xyz,etot,chrg,gen)
     implicit none
     class(ilist) :: self
     type(msmol) :: mol
     integer :: nat
     real(wp) :: xyz(3,nat)
     integer :: at(nat)
     real(wp) :: etot
     integer :: chrg
     integer :: gen

     type(msmol),allocatable :: dummy(:)
     logical,allocatable     :: btmp(:)
     integer :: n
     call mol%newmol(nat,xyz,at,etot,chrg,gen)
     if(.not.allocated(self%mol))then
         self%nmol = 1
         allocate(self%mol(1))
         allocate(self%new(1))
         self%mol(1) = mol
         self%new(1) = .true.
     else
         n = self%nmol + 1
         allocate(dummy(n))
         allocate(btmp(n))
         dummy(1:self%nmol) = self%mol(1:self%nmol)
         dummy(n) = mol
         btmp(1:self%nmol) = self%new(1:self%nmol)
         btmp(n) = .true.
         self%nmol = n
         call move_alloc(btmp, self%new)
         call move_alloc(dummy, self%mol)
     endif
     return
end subroutine ilist_append


end module msmod
