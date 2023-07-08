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
subroutine discretize_trj(env)
  use crest_parameters
  use crest_data
  use strucrd
  use zdata,only:readwbo
  use adjacency
  use discretize_module
  implicit none
  !> INPUT/OUTPUT
  type(systemdata),intent(inout) :: env
  !> LOCAL
  type(coord) :: mol  !> by convention mol is in Bohrs
  real(wp),allocatable :: wbo(:,:)
  integer,allocatable  :: na(:),nb(:),nc(:)
  real(wp),allocatable :: zmat(:,:),zmat_new(:,:)
  integer,allocatable  :: Amat(:,:)
  integer :: i,j,k,l,ich
  character(len=124) :: stmp
  !> number of dihedral angles that are adjusted
  !> this defines the overall number of conformers
  integer :: ndieder
  integer :: ncombi
  real(wp) :: ncombi_f
  !> data for new structure and reference checks
  real(wp) :: cthr,p 
  integer :: nremain
  !> trajectory data
  integer :: nall,nat
  integer,allocatable :: trjat(:)
  real(wp),allocatable :: trjxyz(:,:,:) 
  real(wp),allocatable :: trje(:)
  real(wp),allocatable :: trjzmat(:,:,:)

  !> discretization data  
  type(struc_info),allocatable :: discinfo(:)

!========================================================================================!
   cthr = 0.3d0 !> CN clash threshold

!>--- get reference mol
   call env%ref%to( mol ) 

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
!>--- analyze dihedral angles encoded in the zmatrix
   !TODO this is also where to select dihedral angles to put into the
   !discretization routines. For now I just take any
   if(.true.)then
    call rigidconf_count_fallback(mol%nat,na,nb,nc,wbo,ndieder)
    if (ndieder < 1) stop 'no dihedral angles detected!'
   end if
  

!========================================================================================!
  write (stdout,*)
  call smallhead('Reading trajectory')
  !>---- read the input ensemble
  call rdensembleparam(env%ensemblename,nat,nall)
  if (nall .lt. 1) return
  allocate (trjxyz(3,nat,nall),trjat(nat),trje(nall))
  call rdensemble(env%ensemblename,nat,nall,trjat,trjxyz,trje)
  write (stdout,'(a,3(i0,a))') '> ',nall,' structures with ',nat,' atoms and ', &
  & ndieder,' dihedral angles'

  allocate(trjzmat(3,nat,nall),source=0.0_wp)
  write (stdout,'(a)',advance='no') '> converting to internal coordinates ... '
  flush(stdout)
  do i=1,nall
     call XYZGEO2(nat,trjxyz,NA,NB,NC,1.0_wp,zmat)  
     trjzmat(:,:,i) = zmat(:,:)
  enddo
  write (stdout,'(a)') 'done.'

!========================================================================================!
  write (stdout,*)
  call smallhead('Analyzing dihedral angle distributions')
  allocate(discinfo(ndieder))
  do j=1,ndieder
   
  enddo

!========================================================================================!
  if (allocated(trjat)) deallocate(trjat)
  if (allocated(trjxyz)) deallocate(trjxyz)
  if(allocated(trjzmat)) deallocate(trjzmat)
  if (allocated(zmat_new)) deallocate(zmat_new) 
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
   type(struc_info) :: disctest   

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

  !call probability_count_minima(0.0_wp,2.0*pi,kappa,mu,NM) 
  call disctest%setup(kappa,mu)
  NM = disctest%nmax
  write(*,*)
  write(*,*) 'Distirbution has ',NM,' minima'
  write(*,*)  

end subroutine test_vonMises

!========================================================================================!
!========================================================================================!

