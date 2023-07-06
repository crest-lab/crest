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
subroutine discretize_trj(env,mol)
  use crest_parameters
  use crest_data
  use strucrd
  use zdata,only:readwbo
  use adjacency
  implicit none
  !> INPUT/OUTPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol  !> by convention mol is in Bohrs
  !> LOCAL
  real(wp),allocatable :: wbo(:,:)
  integer,allocatable  :: na(:),nb(:),nc(:)
  real(wp),allocatable :: zmat(:,:),zmat_new(:,:)
  integer,allocatable  :: Amat(:,:)
  integer :: i,j,k,l,ich
  character(len=124) :: stmp
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
  !> work space array for the recursive routine and iteration variables
  integer(int8),allocatable :: wcombi(:)
  integer :: iteratord, iteratork

  !> data for new structure and reference checks
  logical,allocatable :: sane(:)
  type(coord) :: newmol
  real(wp),allocatable :: rcov(:),cnref(:)
  real(wp) :: cthr,p 
  integer :: nremain

  character(len=19),parameter :: outputfile = "crest_rigidconf.xyz"

!========================================================================================!
   cthr = 0.3d0 !> CN clash threshold

!========================================================================================!
!>--- generate the topology information
  allocate (Amat(mol%nat,mol%nat),source=0) !> adjacency matrix

  select case (env%rigidconf_toposource)
    !case ( 1 ) !>
    !
  case default !> do a GFN0-xTB singlepoint and read WBOs

    call crest_xtbsp(env,0,mol) !> writes file "wbo"
    allocate (wbo(mol%nat,mol%nat),source=0.0_wp)
    call readwbo("wbo",mol%nat,wbo) !> reads file "wbo"
    call wbo2adjacency(mol%nat,wbo,Amat,0.02_wp) !> setup Amat

  end select

!========================================================================================!
!>--- generate zmatrix based on adjacency
  allocate (na(mol%nat),nb(mol%nat),nc(mol%nat),source=0)
  allocate (zmat(3,mol%nat),source=0.0_wp)
  call BETTER_XYZINT(mol%nat,mol%xyz,Amat,na,nb,nc,zmat)
  call smallhead('Internal coordinates:')
  call print_zmat(stdout,mol%nat,mol%at,zmat,na,nb,nc,.true.)
!--- Note: zmat angles (columns 2 and 3) are in Radians

!========================================================================================!
!>--- analyze dihedral angles encoded in the zmatrix to select the ones to vary
!>    and how to vary them

  !=================================================!
  !> IMPLEMENTATION (or subroutine call) GOES HERE <!
  !=================================================!

!>--- fallback implementation for testing: All single-bonds with a corresponding
!>    entry in the zmatrix (this excludes terminal atoms, e.g. H)
  if (.true.) then
    call rigidconf_count_fallback(mol%nat,na,nb,nc,wbo,ndieder)
    if (ndieder < 1) stop 'no dihedral angles selected!'
    allocate (dvalues(ndieder),source=0)
    allocate (dstep(ndieder),source=0.0_wp)
    allocate (ztod(mol%nat),source=0)
    call rigidconf_analyze_fallback(env,mol,zmat,na,nb,nc,wbo, &
    &                               ndieder,dvalues,dstep,ztod)
    
    !call prune_zmat_dihedrals(mol%nat, mol%xyz, zmat, na,nb,nc, ztod )
    !call smallhead('New internal coordinates:')
    !call print_zmat(stdout,mol%nat,mol%at,zmat,na,nb,nc,.true.)
   end if

!========================================================================================!
  write (stdout,*)
  call smallhead('Discretization of dihedral angles')




!========================================================================================!
  if (allocated(sane)) deallocate(sane)
  if (allocated(zmat_new)) deallocate(zmat_new) 
  if (allocated(wcombi)) deallocate(wcombi)
  if (allocated(combi)) deallocate (combi)
  if (allocated(dvalues)) deallocate (dvalues)
  if (allocated(dstep)) deallocate (dstep)
  if (allocated(ztod)) deallocate (ztod)
  if (allocated(wbo)) deallocate (wbo)
  if (allocated(na)) deallocate (na)
  if (allocated(nb)) deallocate (nb)
  if (allocated(nc)) deallocate (nc)
  if (allocated(zmat)) deallocate (zmat)
  if (allocated(Amat)) deallocate (Amat)
  return
end subroutine discretize_trj
!========================================================================================!

subroutine test_vonMises(env,kappa,n,mu)
   use crest_parameters
   use crest_data
   use probabilities_module
   use discretize_module
   implicit none
   type(systemdata),intent(inout) :: env
   integer,intent(in) :: n
   real(wp),intent(in) :: mu(n)
   real(wp),intent(in) ::  kappa
   real(wp) :: theta,p,dpdt,dpdt2,numdpdt,pp,pm,tmp
   real(wp),parameter :: d=1.0d-1
   integer :: i,j,k,l, NM
   integer :: ich

   !> Plot an example von Mises distribution

   theta = 0.0_wp
   open(newunit=ich,file='vonmises.txt')
   do i=1,360
     theta = theta + degtorad
    call vonMises(theta,kappa,mu,p,dpdt,dpdt2)
    !write(ich,'(4f16.8)') theta*radtodeg,p,dpdt,dpdt2
     write(ich,'(3f16.8)') p,dpdt,dpdt2
   enddo
   close(ich)

  call probability_count_minima(0.0_wp,2.0*pi,kappa,mu,NM) 
  write(*,*)
  write(*,*) 'Distirbution has ',NM,' minima'
  write(*,*)  

end subroutine test_vonMises

!========================================================================================!
!========================================================================================!

