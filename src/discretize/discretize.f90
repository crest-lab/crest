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
  use INTERNALS_mod
  use rigidconf_analyze
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
  integer,allocatable  :: ztod(:)
  integer :: zmax
  integer :: i,j,k,l,m,ich
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
  integer,allocatable  :: drep(:)              !> mapping zmat entry to unique dihedral
  type(struc_info),allocatable :: discinfo(:)  !> discretization data (position min/max, prob)
  real(wp),allocatable :: mutmp(:)             !> temporary data array
  integer,allocatable :: datmat(:,:)           !> data matrix, discretized values
  integer,allocatable :: dummat(:,:)
  real(wp) :: dtmp

!========================================================================================!
  cthr = 0.3d0 !> CN clash threshold

!>--- get reference mol
  call env%ref%to(mol)

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
  call smallhead('Internal coordinates of reference structure:')
  call print_zmat(stdout,mol%nat,mol%at,zmat,na,nb,nc,.true.)
!--- Note: zmat angles (columns 2 and 3) are in Radians

!========================================================================================!
!>--- analyze dihedral angles encoded in the zmatrix
  !TODO this is also where to select dihedral angles to put into the
  !discretization routines. For now I just take any
  if (.true.) then
    zmax = mol%nat
    allocate (ztod(zmax),source=0)
    call rigidconf_count_fallback(mol%nat,na,nb,nc,wbo,ndieder,ztod)
    if (ndieder < 1) stop 'no dihedral angles detected!'

    allocate (drep(ndieder),source=0)
    do j = 1,ndieder
      do i = 1,zmax
        if (ztod(i) == j) then
          drep(j) = i
          exit
        end if
      end do
    end do
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

  allocate (trjzmat(3,nat,nall),source=0.0_wp)
  write (stdout,'(a)',advance='no') '> converting to internal coordinates ... '
  flush (stdout)
  do i = 1,nall
    call XYZGEO2(nat,trjxyz(:,:,i),NA,NB,NC,1.0_wp,zmat)
    trjzmat(:,:,i) = zmat(:,:)
  end do
  write (stdout,'(a)') 'done.'

!========================================================================================!
  write (stdout,*)
  call smallhead('Analyzing dihedral angle distributions')
  allocate(mutmp(nall))
  allocate (discinfo(ndieder))
  do j = 1,ndieder
    m = drep(j)
    mutmp(:) = trjzmat(3,m,:)
    call discinfo(j)%setup(env%kappa,mutmp)
    call discinfo(j)%info(j,radtodeg)
  end do
  deallocate(mutmp)

  write(stdout,*)
  write(stdout,'(a,a)') '> ',repeat('*',73)
  write(stdout,'(a,f18.8,a)') '> Total marginal (first order) entropy :',sum(discinfo(:)%S),' cal mol⁻¹ K⁻¹'
  write(stdout,'(a,a)') '> ',repeat('*',73)

!========================================================================================!
  allocate(datmat(ndieder,nall), source = 0)
  write (stdout,*)
  write (stdout,'(a)',advance='no') '> discretizing ensemble data ... '
  flush (stdout)
  do i=1,nall
    do j=1,ndieder
       m = drep(j)
       dtmp = trjzmat(3,m,i)
       datmat(j,i) = discinfo(j)%discrete( dtmp )
    enddo
  enddo 
  write (stdout,'(a)') 'done.'
  
  call write_data_matrix_bin(datmat)
  write (stdout,'(a)') '> written discretized data to binary file ddata.bin'
  !call read_data_matrix_bin(i,j,dummat)
  !write(*,*) i,j,all(datmat==dummat)

!========================================================================================!
  if (allocated(datmat)) deallocate(datmat)
  if (allocated(drep)) deallocate(drep)
  if (allocated(trjat)) deallocate (trjat)
  if (allocated(trjxyz)) deallocate (trjxyz)
  if (allocated(trjzmat)) deallocate (trjzmat)
  if (allocated(zmat_new)) deallocate (zmat_new)
  if (allocated(wbo)) deallocate (wbo)
  if (allocated(na)) deallocate (na)
  if (allocated(nb)) deallocate (nb)
  if (allocated(nc)) deallocate (nc)
  if (allocated(zmat)) deallocate (zmat)
  if (allocated(Amat)) deallocate (Amat)
  return
end subroutine discretize_trj

!========================================================================================!
!========================================================================================!
