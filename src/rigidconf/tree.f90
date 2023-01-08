!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht, Christopher Zurek, Christoph Bannwarth
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

!========================================================================================!
!========================================================================================!
subroutine rigidconf_tree(env,mol)
!************************************************************
!* Generate conformers in a straight-forward approach.
!* Selected dihedral angles may have a set of possible values.
!* All combinations of these values will be generated.
!* The generation process follows a simple recursive strategy,
!* which can be interpreted as setting up a tree graph
!* e.g., A,B,C,... being dihedral angles:       
!*
!*  ref.
!*   └──┐ 
!*      ├──A1 
!*      │  ├──B1
!*      ┊  │  ├──C1
!*      ┊  │     ...
!*      ┊  ├──B2
!*      ┊  ┊  ├──C1  
!*      ┊        ...
!*      ├─A2 
!*      │  ├──B1
!*      ┊  ┊  ├─┄...
!*      ┊     ...
!*      etc.
!*
!* Input:
!*    env  - CREST's systemdata
!*    mol  - Input (reference) geometry
!*
!************************************************************
  use crest_parameters
  use crest_data
  use strucrd
  use zdata, only: readwbo
  use adjacency
  implicit none
  !> INPUT/OUTPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol  !> by convention mol is in Bohrs
  !> LOCAL
  real(wp),allocatable :: wbo(:,:)
  integer,allocatable  :: na(:),nb(:),nc(:)
  real(wp),allocatable :: zmat(:,:)
  integer,allocatable  :: Amat(:,:)
  integer :: i,j,k,l

  !> number of dihedral angles that are adjusted
  !> this defines the overall number of conformers
  integer :: ndieder 
  !> the number of values each of the dihedral angles
  !> can assume (equally spaced in -180°< ϕ <180° )
  !> i.e., if dvalues(i)=3, the dihedral i will have 3 defined values,
  !> ϕ (as in the input geometry), ϕ+120°, and ϕ+240°
  integer,allocatable :: dvalues(:)
  !> the step size for each dihedral angle, in Radians
  real(wp),allocatable :: dstep(:)
  !> mapping zmat entry to selected dihedral angles. multiple zmat entries 
  !> may coincide with the same bond, making this mapping necessary.
  integer,allocatable :: ztod(:) 
  !> saved combinations of dihedral angles in RAM. We first generate
  !> these "fingerprints" using the recursive strategy, and only
  !> later reconstruct the molecular structures. This allows running
  !> the structure checks in parallel
  integer :: ncombi
  real(wp) :: ncombi_f
  integer(int8),allocatable :: combi(:,:)

!========================================================================================!
!>--- generate the topology information
  allocate(Amat(mol%nat,mol%nat), source=0) !> adjacency matrix

  select case( env%rigidconf_toposource ) 
  !case ( 1 ) !> 
  !  
  case default !> do a GFN0-xTB singlepoint and read WBOs
   call crest_xtbsp(env,0,mol) !> writes file "wbo"
   allocate( wbo(mol%nat,mol%nat), source=0.0_wp)
   call readwbo("wbo",mol%nat,wbo) !> reads file "wbo"
   call wbo2adjacency(mol%nat, wbo, Amat, 0.02_wp) !> setup Amat
  end select

!========================================================================================!
!>--- generate zmatrix based on adjacency
  allocate(na(mol%nat), nb(mol%nat), nc(mol%nat), source=0)
  allocate(zmat(3,mol%nat), source=0.0_wp)
  call BETTER_XYZINT( mol%nat, mol%xyz, Amat, na, nb, nc, zmat)
  call smallhead('Internal coordinates:')
  call print_zmat(stdout,mol%nat,zmat,na,nb,nc)
!--- Note: zmat angles (columns 2 and 3) are in Radians

!========================================================================================!
!>--- analyze dihedral angles encoded in the zmatrix to select the ones to vary
!>    and how to vary them

   !=================================================!
   !> IMPLEMENTATION (or subroutine call) GOES HERE <!
   !=================================================!

!>--- fallback implementation for testing: All single-bonds with a corresponding
!>    entry in the zmatrix (this excludes terminal atoms, e.g. H)
   if(.true.)then
     call rigidconf_count_fallback(mol%nat,na,nb,nc,wbo,ndieder)
     if(ndieder < 1) stop 'no dihedral angles selected!'
     allocate(dvalues(ndieder), source=0)
     allocate(dstep(ndieder), source=0.0_wp)
     allocate(ztod(mol%nat), source=0)
     call rigidconf_analyze_fallback(env,mol,zmat,na,nb,nc,wbo, &
     &                               ndieder,dvalues,dstep,ztod )
   endif 

!========================================================================================!
!>--- For now, I am limiting the max. number of dihedral angles to 10.
!>    At some point we will have to think about what to allow.  
   if( ndieder > 10 )then
     error stop 'Too many dihedral angles. '
   endif
   
!========================================================================================!
!>--- Allocate space for combination "fingerprints"
   ncombi_f = float(dvalues(1))
   do i=2,ndieder   
     ncombi_f = ncombi_f * float(dvalues(i))
   enddo
   if(ncombi_f .gt. float(huge(ncombi)))then
     error stop 'cannot allocate combi()'
   else
     ncombi = nint(ncombi_f)
     allocate( combi(ndieder,ncombi), source = 1_int8)
     write(stdout,*)
     write(stdout,'(">",1x,i0,a,i0,a)') ndieder,' dihedral angles resulting in ', &
     &     ncombi,' possible conformers.'
   endif

!>--- Generate combinations recursively



!========================================================================================!
  if(allocated(combi)) deallocate(combi)
  if(allocated(dvalues)) deallocate(dvalues)
  if(allocated(dstep)) deallocate(dstep)
  if(allocated(ztod)) deallocate(ztod)
  if(allocated(wbo)) deallocate(wbo)
  if(allocated(na)) deallocate(na)
  if(allocated(nb)) deallocate(nb)
  if(allocated(nc)) deallocate(nc)
  if(allocated(zmat)) deallocate(zmat)
  if(allocated(Amat)) deallocate(Amat)
  return
end subroutine rigidconf_tree
!========================================================================================!
!========================================================================================!
