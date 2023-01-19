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
!* which can be interpreted as setting up an undirected tree graph
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
    !call print_zmat(stdout,mol%nat,zmat,na,nb,nc)
   end if


!========================================================================================!
!>--- For now, I am limiting the max. number of dihedral angles to 10.
!>    At some point we will have to think about what to allow.
  if (ndieder > 10) then
    error stop 'Too many dihedral angles. '
  end if

!========================================================================================!
  write (stdout,*)
  call smallhead('Rule-based conformer genration')

!>--- Allocate space for combination "fingerprints"
  ncombi_f = float(dvalues(1))
  do i = 2,ndieder
    ncombi_f = ncombi_f*float(dvalues(i))
  end do
  if (ncombi_f .gt. float(huge(ncombi))) then
    error stop 'cannot allocate combi()'
  else
    ncombi = nint(ncombi_f)
    allocate (combi(ndieder,ncombi),source=1_int8)
    write (stdout,'(">",1x,i0,a,i0,a)') ndieder,' dihedral angles resulting in ', &
    &     ncombi,' possible conformers.'
  end if

!>--- Generate combinations recursively
  write(stdout,'(">",1x,a)',advance='no') 'Generating combination references ...'
  flush(stdout)
  allocate (wcombi(ndieder),source=1_int8)
  iteratord = 1
  iteratork = 0
  call generate_dihedral_combinations(ndieder,wcombi,ncombi,combi, &
  &                    dvalues,iteratord,iteratork)
  write(stdout,*) 'done.'


!========================================================================================!
!>--- Calculate reference CNs
  allocate (rcov(94),cnref(mol%nat),source=0.0_wp)
  call setrcov(rcov)
  call ycoord(mol%nat,rcov,mol%at,mol%xyz,cnref,100.0d0) !> refernce CNs
  
!>--- Generate all the conformers and check for CN clashes
  allocate(zmat_new(3,mol%nat), source=0.0_wp)
  allocate(sane(ncombi), source=.false.)
  
  open(newunit=ich, file=outputfile)
  write(stdout,'(">",1x,a)',advance='no') 'Reconstructing structures and checking for CN clashes ...'
  flush(stdout)
  do i=1,ncombi

!>--- Generate new zmat from zmat and combi entry
    call construct_new_zmat(mol%nat,zmat,combi(:,i),ndieder,dvalues,dstep,ztod,zmat_new) 

!>--- Reconstruct Cartesian coordinates for new zmat
    call reconstruct_zmat_to_mol(mol%nat,mol%at,zmat_new,na,nb,nc,newmol)

!>--- Check new structure for CN clashes w.r.t. reference and cthr
    call ycoord2(newmol%nat,rcov,newmol%at,newmol%xyz,cnref,100.d0,cthr,sane(i)) !> CN clashes
    sane(i) = .not.sane(i)

!>--- Dump to file
    if(sane(i))then
    call newmol%append(ich)
    endif
    
  enddo
  write(stdout,*) 'done.'
  close(ich) 

  nremain = count(sane,1)
  p = float(nremain)/float(ncombi) * 100.0_wp
  write(stdout,'(">",1x,i0,1x,a,1x,i0,1x,a,f5.1,a)') nremain,'of',ncombi, &
  & 'structures remaining (',p,'%). Clashes discarded.'
  write(stdout,'(">",1x,a,1x,a)') 'Written to file',outputfile

!========================================================================================!
!>--- Post-processing

  !=================================================!
  !> IMPLEMENTATION (or subroutine call) GOES HERE <!
  !=================================================!

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
!========================================================================================!
contains
!========================================================================================!
  recursive subroutine generate_dihedral_combinations(ndieder,wcombi,ncombi,combi, &
  &                    dvalues,iteratord,iteratork)
    implicit none
    integer,intent(in) :: ndieder
    integer(int8),intent(inout) :: wcombi(ndieder)
    integer,intent(in) :: ncombi
    integer(int8),intent(inout) :: combi(ndieder,ncombi)
    integer,intent(in) :: dvalues(ndieder)
    integer,intent(inout) :: iteratord,iteratork

    integer :: i,j,k

    if (iteratork > ncombi) return
!>-- if we reached the last dihedral in the list, dump and return
    if (iteratord > ndieder) then
      iteratork = iteratork+1
      combi(:,iteratork) = wcombi(:)
      return
    end if

!>-- iterate the current dihedral angle (j)
    j = iteratord
    k = iteratord+1
    do i = 1,dvalues(j)
!>--- set the combination value to
      wcombi(j) = int(i,int8)
!>--- go to the next dihedral angle (k)
      call generate_dihedral_combinations(ndieder,wcombi,ncombi,combi, &
 &                    dvalues,k,iteratork)
!>--- reset after the last value for the current dihedral angle
      if (i == dvalues(j)) then
        wcombi(j) = 1_int8
      end if
    end do

  end subroutine generate_dihedral_combinations
!========================================================================================!

end subroutine rigidconf_tree
!========================================================================================!
!========================================================================================!

