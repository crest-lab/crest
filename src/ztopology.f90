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

!================================================================================!
!  This file contains many subroutines related to topology setup of a molecule.
!  It requires the polymorphic datatype "zmolecule", which can be found in the
!  "zdata.f90" module.
!
!  Most routines are based on a representation where the molecule is seen as a
!  "graph" and the topology is analyzed by navigating through the different
!  paths of the graph.
!
!  P. Pracht, 2019/20
!================================================================================!
!================================================================================!
!  SIMPLETOPO: set up topology.
!  All the information is stored in the complex datatype object "zmol"
!  Although the name implies it, there is nothing much "simple" for some
!  of the routines ...
!
!  simpeltopo_file version is a file wrapper
!  simpeltopo_mol version is a wrapper for a coord object
!================================================================================!
subroutine simpletopo_file(fname,zmol,verbose,getrings,wbofile)
  use iso_fortran_env,wp => real64
  use zdata
  use strucrd,only:rdnat,rdcoord
  implicit none
  character(len=*) :: fname
  type(zmolecule)  :: zmol
  logical          :: verbose
  logical          :: getrings
  character(len=*),optional :: wbofile
  integer :: n
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:)
  logical :: ex
  interface
    subroutine simpletopo(n,at,xyz,zmol,verbose,getrings,wbofile)
      import :: zmolecule,wp
      implicit none
      type(zmolecule)  :: zmol
      logical          :: verbose
      logical          :: getrings
      integer,intent(in)  :: n
      integer,intent(in)  :: at(n)
      real(wp),intent(in) :: xyz(3,n) !in Bohrs
      character(len=*),intent(in),optional :: wbofile
    end subroutine simpletopo
  end interface

  call rdnat(fname,n)
  allocate (at(n),xyz(3,n))
  call rdcoord(fname,n,at,xyz)

  ex = .false.
  if (present(wbofile)) then
    inquire (file=wbofile,exist=ex)
  end if
  if (.not.ex) then
    call simpletopo(n,at,xyz,zmol,verbose,getrings,'')
  else
    call simpletopo(n,at,xyz,zmol,verbose,getrings,wbofile)
  end if

  deallocate (xyz,at)
  return
end subroutine simpletopo_file

subroutine simpletopo_mol(mol,zmol,verbose,getrings)
  use iso_fortran_env,wp => real64
  use zdata
  use strucrd
  implicit none
  type(coord)      :: mol    !in
  type(zmolecule)  :: zmol   !out
  logical          :: verbose
  logical          :: getrings
  interface
    subroutine simpletopo(n,at,xyz,zmol,verbose,getrings,wbofile)
      import :: zmolecule,wp
      implicit none
      type(zmolecule)  :: zmol
      logical          :: verbose
      logical          :: getrings
      integer,intent(in)  :: n
      integer,intent(in)  :: at(n)
      real(wp),intent(in) :: xyz(3,n) !in Bohrs
      character(len=*),intent(in),optional :: wbofile
    end subroutine simpletopo
  end interface
  call simpletopo(mol%nat,mol%at,mol%xyz,zmol,verbose,getrings,'')
  return
end subroutine simpletopo_mol

!================================================================================!
!  SIMPLETOPO: set up topology.
!  n   - number of atoms
!  at  - atom types (as atom number)
!  xyz - cartesian coordinates (in Bohrs)
!  zmol - polymorphic datatype containing the molecule data
!  verbose - boolean to activate printouts
!  wbofile - (optional) name of the file containing WBOs
!================================================================================!
subroutine simpletopo(n,at,xyz,zmol,verbose,getrings,wbofile)
  use crest_parameters
  use zdata
  use miscdata,only:rcov
  implicit none
  type(zmolecule)  :: zmol
  type(zmolecule)  :: zfrag
  type(zring)      :: newring
  logical          :: verbose
  logical          :: getrings
  integer,intent(in)  :: n
  integer,intent(in)  :: at(n)
  real(wp),intent(in) :: xyz(3,n) !in Bohrs
  character(len=*),intent(in),optional :: wbofile
  type(zatom),allocatable :: zat(:)
  real(wp),allocatable :: cn(:)
  real(wp),allocatable :: bond(:,:)
  real(wp),allocatable :: wbo(:,:)
  integer :: i,j,k
  integer :: ntopo
  logical,allocatable :: neighmat(:,:)
  integer,allocatable :: topovec(:)
  integer :: nrings

  logical :: ex,useWBO

!>--- header
  if (verbose) then
    write (*,*)
    call smallhead('TOPOLOGY ANALYSIS')
  end if

!>--- set covalent radii and calculate coordination numbers
  allocate (cn(n),bond(n,n))
  call xcoord2(n,at,xyz,rcov,cn,900.0_wp,bond)
!>--- read in WBOs if required
  ex = .false.
  if (present(wbofile)) then
    inquire (file=wbofile,exist=ex)
  end if
  if (ex) then
    allocate (wbo(n,n))
    call readwbo("wbo",n,wbo)
    !write(*,*) wbo
    !stop
    useWBO = .true.
  else
    useWBO = .false.
  end if

  allocate (zat(n))
  if (.not.useWBO) then
!--- get individual neighbour lists and set up molecule object "zat"
    ntopo = n*(n+1)/2
    allocate (neighmat(n,n),source=.false.)
    allocate (topovec(ntopo))
    call bondtotopo(n,at,bond,cn,ntopo,topovec,neighmat)
    call neighbourset(zmol,n,at,xyz,cn,ntopo,topovec)
    deallocate (topovec,neighmat)
  else
!--- neighbour list could also be set up from WBOs, which is not necessarily better
    do i = 1,n
      zat(i)%cn = cn(i)
      call wboneighbours(i,n,at,xyz,wbo,0.7d0,zat(i))
      call wborepaireta(i,n,at,xyz,cn,wbo,zat(i))  !handle eta-coordinated atoms
    end do

    !--- transfer information to zmol class
    zmol%nat = n      !number of atoms
    zmol%zat = zat    !list of atoms, including data
    zmol%at = at      !list of atom types (integer)
  end if
  do i = 1,n
    call zat(i)%deallocate()
  end do
  !if(allocated(zat))deallocate(zat)

!--- analyze system for fragments
  call zmol%mrec()

  allocate (zmol%distmat(n,n))
  do i = 1,n
    do j = 1,n
      zmol%distmat(i,j) = (xyz(1,i)-xyz(1,j))**2+ &
    &                   (xyz(2,i)-xyz(2,j))**2+ &
    &                   (xyz(3,i)-xyz(3,j))**2
      zmol%distmat(i,j) = sqrt(zmol%distmat(i,j))
    end do
  end do

!--- printouts?
  if (verbose) then
    write (*,'(1x,a)') 'Coordinates (Bohr):'
    do i = 1,n
      call zmol%zat(i)%wrtmline(6)
    end do
    write (*,*)
    call zmol%wrtable(6)
    write (*,*)

    if (allocated(zmol%molvec)) then
      write (*,'(1x,a,i0)') 'Number of fragments in the system:  ',zmol%nfrag
      write (*,*)
    end if
  end if

!>--- identify rings
  if (getrings) then
    do i = 1,zmol%nfrag
      call zmol%fragment(i,zfrag)
      call newgetrings(zfrag,.false.)
      nrings = zfrag%nri
      do j = 1,nrings
        newring = zfrag%zri(j)
        do k = 1,newring%rs
          newring%rlist(k) = zfrag%revmap(newring%rlist(k))
        end do
        call zmol%addring(newring)
        call newring%deallocate()
      end do
      !end if
      if (verbose.and.nrings > 0) then
        if (zmol%nfrag > 1) then
          write (stdout,'(1x,a,i0)') 'Fragment ',i
          write (stdout,'(1x,a,i0,/)') 'Total number of rings in the fragment: ',nrings
        else
          write (stdout,'(1x,a,i0,/)') 'Total number of rings in the system: ',nrings
        end if
      end if
      call zfrag%deallocate()
    end do
    if (verbose) then
      call zmol%prrings(stdout)
      if (zmol%nri > 0) then
        write (stdout,'(/,1x,a,i0,/)') 'Total number of rings in the system: ',zmol%nri
      end if
    end if
  end if

!>--- deallocation of memory
  if (allocated(wbo)) deallocate (wbo)
  deallocate (bond,cn)
  return
end subroutine simpletopo

!=======================================================================!
!> compute coordination numbers by adding an inverse damping function
!=======================================================================!
subroutine xcoord2(nat,iz,xyz,rcov,cn,cn_thr,bond)
  use iso_fortran_env,wp => real64
  implicit none
  integer,intent(in) :: nat
  integer,intent(in) :: iz(nat)
  real(wp),intent(in) :: xyz(3,nat)
  real(wp),intent(out) :: cn(nat)
  real(wp),intent(in)  :: cn_thr
  real(wp),intent(in)  :: rcov(*)
  real(wp),intent(out) :: bond(nat,nat)
  integer :: i,k1
  integer :: iat
  real(wp) :: dx,dy,dz,r,damp,xn,rr,rco,r2,rcovi,rcovj
  k1 = 16
  bond = 0.0d0
  cn = 0.0d0
  do i = 1,nat
    xn = 0.0d0
    rcovi = rcov(iz(i))
    do iat = 1,nat
      if (iat .ne. i) then
        dx = xyz(1,iat)-xyz(1,i)
        dy = xyz(2,iat)-xyz(2,i)
        dz = xyz(3,iat)-xyz(3,i)
        r2 = dx*dx+dy*dy+dz*dz
        r = sqrt(r2)
        if (r2 .gt. cn_thr) cycle
        rcovj = rcov(iz(iat))
!> covalent distance in Bohr
        rco = (rcovi+rcovj)*1.0  ! this scaling reduces the size of the clusters
        rr = rco/r
!> counting function exponential has a better long-range behavior than MHGs inverse damping
        damp = 1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
        bond(iat,i) = damp
        xn = xn+damp
      end if
    end do
    cn(i) = xn
  end do
  return
end subroutine xcoord2

!===================================================!
!> generate the topo array for a given structure
!===================================================!
subroutine bondtotopo(nat,at,bond,cn,ntopo,topo,neighbourmat)
  use iso_fortran_env,only:wp => real64
  use utilities
  integer,intent(in)  :: nat
  integer,intent(in) :: at(nat)
  real(wp),intent(inout) :: bond(nat,nat)
  real(wp),intent(in) :: cn(nat)
  integer,intent(in)  :: ntopo
  integer,intent(out) :: topo(ntopo)
  real(wp),allocatable :: cn2(:)
  logical,intent(inout) :: neighbourmat(nat,nat)
  integer :: i,j,k,l
  integer :: icn,rcn
  allocate (cn2(nat),source=0.0_wp)
  topo = 0
  neighbourmat = .false.
  !>--- some heuristic rules and CN array setup
  do i = 1,nat
    cn2(i) = cn(i)
    rcn = floor(cn(i))
    select case (at(i)) !additional empirical topology rules
      ! case( 5 ) !B
      !   if( nint(cn(i)) > 4) cn2(i)=4.0_wp
      ! case( 9,17,35,53 ) !F,Cl,Br,I
      !   cn2(i) = min(cn(i),1.0_wp)
    case (6) !C
      if ((cn(i)-rcn) < 0.7_wp) then
        cn2(i) = rcn
      end if
    end select
    !-- extreme CN cases
    if (nint(cn(i)) > 8) cn2(i) = 8.0_wp
    !empirical: rounding down up to .6 is better for topo setup
    if ((cn(i)-rcn) < 0.6_wp) then
      cn2(i) = rcn
    end if
  end do
  !>--- build the topology
  do i = 1,nat
    icn = nint(cn2(i))
    do k = 1,icn
      j = maxloc(bond(:,i),1)
      bond(j,i) = 0.0d0
      if (i .eq. j) cycle
      neighbourmat(i,j) = .true. !--important: not automatically (i,j)=(j,i)
    end do
  end do
  do i = 1,nat
    do j = 1,nat
      if (i == j) cycle
      l = lin(i,j)
      !>-- only save matching topology --> prevent high CN failures
      if (neighbourmat(i,j).and.neighbourmat(j,i)) then
        topo(l) = 1
      else
        !> special case for carbon (because the carbon CN is typically correct)
        !> this helps, e.g. with eta-coordination in ferrocene
        !> (used, except if both are carbon)
        if (.not. (at(i) == 6.and.at(j) == 6)) then
          if (at(i) == 6.and.neighbourmat(i,j)) topo(l) = 1
          if (at(j) == 6.and.neighbourmat(j,i)) topo(l) = 1
        end if
      end if
    end do
  end do
  deallocate (cn2)
  return
end subroutine bondtotopo

subroutine bondtotopo_excl(nat,at,bond,cn,ntopo,topo,neighbourmat,excl)
  use iso_fortran_env,only:wp => real64
  use utilities
  integer,intent(in)  :: nat
  integer,intent(in) :: at(nat)
  real(wp),intent(inout) :: bond(nat,nat)
  real(wp),intent(in) :: cn(nat)
  integer,intent(in)  :: ntopo
  integer,intent(out) :: topo(ntopo)
  real(wp),allocatable :: cn2(:)
  logical,intent(inout) :: neighbourmat(nat,nat)
  logical,intent(in) :: excl(nat)
  integer :: i,j,k,l
  integer :: icn,rcn
  allocate (cn2(nat),source=0.0_wp)
  topo = 0
  neighbourmat = .false.
  !--- some heuristic rules and CN array setup
  do i = 1,nat
    cn2(i) = cn(i)
    rcn = floor(cn(i))
    select case (at(i)) !additional empirical topology rules
      ! case( 5 ) !B
      !   if( nint(cn(i)) > 4) cn2(i)=4.0_wp
      ! case( 9,17,35,53 ) !F,Cl,Br,I
      !   cn2(i) = min(cn(i),1.0_wp)
    case (6) !C
      if ((cn(i)-rcn) < 0.7_wp) then
        cn2(i) = rcn
      end if
    end select
    !-- extreme CN cases
    if (nint(cn(i)) > 8) cn2(i) = 8.0_wp
    !empirical: rounding down up to .6 is better for topo setup
    if ((cn(i)-rcn) < 0.6_wp) then
      cn2(i) = rcn
    end if
  end do
  !--- build the topology
  do i = 1,nat
    icn = nint(cn2(i))
    do k = 1,icn
      j = maxloc(bond(:,i),1)
      bond(j,i) = 0.0d0
      if (i .eq. j) cycle
      neighbourmat(i,j) = .true. !--important: not automatically (i,j)=(j,i)
      if (excl(i).or.excl(j)) neighbourmat(i,j) = .false.
    end do
  end do
  do i = 1,nat
    do j = 1,nat
      if (i == j) cycle
      l = lin(i,j)
      !-- only save matching topology --> prevent high CN failures
      if (neighbourmat(i,j).and.neighbourmat(j,i)) then
        topo(l) = 1
      else
        ! special case for carbon (because the carbon CN is typically correct)
        ! this helps, e.g. with eta-coordination in ferrocene
        ! (used, except if both are carbon)
        if (.not. (at(i) == 6.and.at(j) == 6)) then
          if (at(i) == 6.and.neighbourmat(i,j)) topo(l) = 1
          if (at(j) == 6.and.neighbourmat(j,i)) topo(l) = 1
        end if
      end if
    end do
  end do
  deallocate (cn2)
  return
end subroutine bondtotopo_excl
subroutine quicktopo(nat,at,xyz,ntopo,topovec)
  use iso_fortran_env,only:wp => real64
  use miscdata,only:rcov
  implicit none
  integer :: nat
  integer :: at(nat)
  real(wp) :: xyz(3,nat) !must be in Bohrs
  integer :: ntopo
  integer :: topovec(ntopo)
  real(wp),allocatable :: cn(:),bond(:,:)
  logical,allocatable :: neighmat(:,:)
  allocate (bond(nat,nat),cn(nat),source=0.0_wp)
  allocate (neighmat(nat,nat),source=.false.)
  cn = 0.0d0
  bond = 0.0d0
  call xcoord2(nat,at,xyz,rcov,cn,900.0_wp,bond)
  call bondtotopo(nat,at,bond,cn,ntopo,topovec,neighmat)
  deallocate (neighmat,cn,bond)
  return
end subroutine quicktopo

!-- transfer topology data to the zmol object
subroutine neighbourset(zmol,nat,at,xyz,cn,ntopo,topovec)
  use iso_fortran_env,only:wp => real64
  use zdata
  use strucrd,only:i2e
  use utilities
  implicit none
  type(zmolecule) :: zmol
  integer,intent(in)  :: nat
  real(wp),intent(in) :: xyz(3,nat)
  integer,intent(in)  :: at(nat)
  real(wp),intent(in) :: cn(nat)
  integer,intent(in)     :: ntopo
  integer,intent(in) :: topovec(ntopo)
  integer :: i,j,k,l
  integer :: inei

  allocate (zmol%zat(nat))
  zmol%nat = nat    !number of atoms
  zmol%at = at      !list of atom types (integer)
  do i = 1,nat
    zmol%zat(i)%atype = at(i)
    zmol%zat(i)%cart(:) = xyz(:,i)
    zmol%zat(i)%el = i2e(at(i),'nc')
    zmol%zat(i)%cn = cn(i)
    zmol%zat(i)%pos = i
    inei = 0
    do j = 1,nat
      if (i == j) cycle
      l = lin(i,j)
      if (topovec(l) == 1) inei = inei+1
    end do
    allocate (zmol%zat(i)%ngh(inei))
    allocate (zmol%zat(i)%ngt(inei))
    k = 0
    do j = 1,nat
      if (i == j) cycle
      l = lin(i,j)
      if (topovec(l) == 1) then
        k = k+1
        zmol%zat(i)%ngh(k) = j      ! the k-th neighbour of atom i is atom j
        zmol%zat(i)%ngt(k) = at(j)  ! atom j has this atom type
      end if
    end do
    zmol%zat(i)%nei = inei
    call quicksort(inei,zmol%zat(i)%ngh)
    !call quicksort(inei,zmol%zat(i)%ngt)
    do j = 1,inei
      zmol%zat(i)%ngt(j) = at(zmol%zat(i)%ngh(j))
    end do
  end do

  return
end subroutine neighbourset

!=======================================================================!
!  compute Zmat (internal coordinates) for the zmol struct
!=======================================================================!
subroutine ztopozmat(zmol,pr)
  use iso_fortran_env,wp => real64
  use zdata
  implicit none
  type(zmolecule)  :: zmol
  logical          :: pr

  integer :: nat
  real(wp),allocatable :: xyz(:,:)
  real(wp),allocatable :: geo(:,:)
  integer,allocatable :: na(:),nb(:),nc(:)
  integer :: i
  real(wp),parameter :: pi = 3.14159265358979D0
  real(wp),parameter :: rad = 180.0d0/pi

  nat = zmol%nat
  allocate (geo(3,nat),xyz(3,nat),source=0.0d0)
  allocate (na(nat),nb(nat),nc(nat),source=0)

  do i = 1,nat
    xyz(1:3,i) = zmol%zat(i)%cart(1:3)
  end do

  !call xyzint(xyz,nat,na,nb,nc,rad,geo)
  call xyzint(xyz,nat,na,nb,nc,1.0d0,geo)
  !call xyzgeo(xyz,nat,na,nb,nc,rad,geo)

  if (pr) then
    call smallhead('INTERNAL COORDINATES (ZMATRIX)')
    do i = 1,nat
      write (*,'(2x,a,4x,3f10.4,1x,3i5)') zmol%zat(i)%el,geo(1:3,i),na(i),nb(i),nc(i)
    end do
  end if

  call move_alloc(geo,zmol%zmat)
  call move_alloc(na,zmol%zna)
  call move_alloc(nb,zmol%znb)
  call move_alloc(nc,zmol%znc)

  deallocate (xyz)
  return
end subroutine ztopozmat

!=======================================================================!
! With the CN-based neighbourlist it can happen that two atoms
! are shared as a neighbour only in the list of one of the atoms.
! This needs correction and is one of the reasons why WBO-based
! neighbourlists are better.
!=======================================================================!
subroutine crosscheckCNnei(zmol)
  use iso_fortran_env,wp => real64
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer :: i,j,k,l
  integer :: nei,newnei
  integer,allocatable :: newngh(:)
  integer,allocatable :: newngt(:)
  do i = 1,zmol%nat
    nei = zmol%zat(i)%nei
    newnei = nei
    allocate (newngh(nei),source=0)
    allocate (newngt(nei),source=0)
    do j = 1,nei
      k = zmol%zat(i)%ngh(j)
      if (any(zmol%zat(k)%ngh(:) .eq. i)) then
        newngh(j) = k
        newngt(j) = zmol%at(k)
        cycle !atom i is also the neighbour of its own neighbour
      else
        !else, the atoms seem to be "artificial" neighbours and
        !has to be removed from the neighbourlist
        newngh(j) = 0
        newnei = newnei-1
      end if
    end do
    !--- if the neighbour list was changed we update it
    if (newnei .lt. nei) then
      deallocate (zmol%zat(i)%ngh,zmol%zat(i)%ngt)
      allocate (zmol%zat(i)%ngh(newnei),zmol%zat(i)%ngt(newnei))
      l = 0
      do j = 1,nei
        if (newngh(j) .ne. 0) then
          l = l+1
          zmol%zat(i)%ngh(l) = newngh(j)
          zmol%zat(i)%ngt(l) = newngt(j)
        end if
      end do
      zmol%zat(i)%nei = newnei
    end if
    deallocate (newngt,newngh)
  end do

  return
end subroutine crosscheckCNnei

!=======================================================================!
!  get the neighbours of an atom based on the WBOs
!=======================================================================!
subroutine wboneighbours(i,nat,at,xyz,wbo,wbothr,zat)
  use iso_fortran_env,only:wp => real64
  use zdata
  use strucrd,only:i2e
  implicit none
  integer,intent(in)  :: i
  real(wp),intent(in) :: xyz(3,nat)
  real(wp),intent(in) :: wbo(nat,nat)
  real(wp),intent(in) :: wbothr
  integer,intent(in)  :: at(nat)
  integer,intent(in)  :: nat
  type(zatom)         :: zat    !--- "zat" is the complex datatype object for atom i
  integer :: icn
  integer :: k,l,j

  zat%atype = at(i)
  zat%cart(:) = xyz(:,i)
  zat%el = i2e(at(i),'nc')
  zat%pos = i

  icn = 0
  do k = 1,nat
    if (wbo(i,k) .ge. wbothr) then
      icn = icn+1
    end if
  end do
  zat%nei = icn
  allocate (zat%ngh(icn),zat%ngt(icn))
  l = 0
  do j = 1,nat
    if (wbo(i,j) .ge. wbothr) then
      l = l+1
      zat%ngh(l) = j      ! the l-th neighbour of atom i is atom j
      zat%ngt(l) = at(j)  ! atom j has this atom type
    end if
  end do
  call quicksort(icn,zat%ngh)
  !call quicksort(icn,zat%ngt)
  do j = 1,icn
    zat%ngt(j) = at(zat%ngh(j))
  end do

  return
end subroutine wboneighbours
subroutine wborepaireta(i,nat,at,xyz,cn,wbo,zat)
  use iso_fortran_env,only:wp => real64
  use zdata
  implicit none
  integer,intent(in)  :: i
  real(wp),intent(in) :: xyz(3,nat)
  real(wp),intent(in) :: cn(nat)
  real(wp),intent(in) :: wbo(nat,nat)
  integer,intent(in)  :: at(nat)
  integer,intent(in)  :: nat
  type(zatom)         :: zat

  if (cn(i) .gt. 1.and.zat%nei .lt. 1) then
    if (allocated(zat%ngt)) deallocate (zat%ngt)
    if (allocated(zat%ngh)) deallocate (zat%ngh)
    write (*,*) i,at(i),cn(i)
    call wboneighbours(i,nat,at,xyz,wbo,0.1d0,zat)
  end if
  return
end subroutine wborepaireta

!=================================================================!
! subroutine wbomrec
! molcount: number of total fragments (increased during search)
! nat: overall number of atoms
! wbo: bond orders
! wbothr: threshold for when a bond is present
! molvec: assignment vector of atom to fragment
!=================================================================!
subroutine wbomrec(molcount,nat,wbo,wbothr,molvec)
  use iso_fortran_env,wp => real64
  implicit none
  real(wp) :: wbo(nat,nat)
  real(wp) :: wbothr
  integer :: nat,molvec(nat),i,molcount
  logical :: taken(nat)
  molvec = 0
  molcount = 1
  taken = .false.
  do i = 1,nat
    if (.not.taken(i)) then
      molvec(i) = molcount
      taken(i) = .true.
      call wbofrags(i,nat,wbo,wbothr,taken,molvec,molcount)
      molcount = molcount+1
    end if
  end do
  molcount = molcount-1
  return
end subroutine wbomrec
!- depth-first search through the WBO list to determine all connected atoms
recursive subroutine wbofrags(i,nat,wbo,wbothr,taken,molvec,molcnt)
  use iso_fortran_env,wp => real64
  implicit none
  integer :: i,nat
  real(wp) :: wbo(nat,nat)
  real(wp) :: wbothr
  integer :: molcnt,molvec(nat),j
  logical taken(nat)
  do j = 1,nat
    if (i .eq. j) cycle
    if (wbo(i,j) .ge. wbothr) then
      if (.not.taken(j)) then
        molvec(j) = molcnt
        taken(j) = .true.
        call wbofrags(j,nat,wbo,wbothr,taken,molvec,molcnt)
      end if
    end if
  end do
  return
end subroutine wbofrags

!=======================================================================!
!  compare two arrays and check if their content is identical
!=======================================================================!
subroutine arrcomp(n,narr,m,marr,equi)
  implicit none
  integer :: n
  integer :: narr(*)
  integer :: m
  integer :: marr(*)
  logical :: equi
  logical,allocatable :: mask(:)
  if (n .ne. m) then
    equi = .false.
    return
  end if
  allocate (mask(n))
  mask(:) = narr(1:n) .eq. marr(1:n)
  if (any(.not.mask(:))) then
    equi = .false.
  else
    equi = .true.
  end if
  deallocate (mask)
  return
end subroutine arrcomp

!=======================================================================!
! compare two arrays and check if they have exactly one element in common
!=======================================================================!
function arruniqel(n,narr,m,marr,el)
  implicit none
  logical :: arruniqel
  integer :: n
  integer :: narr(n)
  integer :: m
  integer :: marr(m)
  integer :: el
  integer :: incr
  integer :: i
  arruniqel = .false.
  el = 0
  incr = 0
  do i = 1,n
    if (any(marr(:) .eq. narr(i))) then
      incr = incr+1
      el = narr(i)
    end if
  end do
  if (incr .eq. 1) then
    arruniqel = .true.
    !el = el
  else
    arruniqel = .false.
    el = 0
  end if
  return
end function arruniqel

!=========================================================================!
! check if all atoms contained in narr are in the neighbour list of atom i
!=========================================================================!
logical function allneigh(zmol,i,n,narr)
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer :: i
  integer :: n
  integer :: narr(n)
  integer :: j
  logical :: dum
  dum = .true.
  do j = 1,n
    dum = dum.and.any(zmol%zat(i)%ngh(:) == narr(j))
  end do
  allneigh = dum
  return
end function allneigh

!=======================================================================!
! get information about the side chain
!=======================================================================!
subroutine getsideweight(zmol,i,j,mside)
  use iso_fortran_env,only:wp => real64
  use zdata
  use atmasses
  implicit none
  integer,intent(in) :: i  !central atom
  integer,intent(in) :: j  !first atom of the side chain
  type(zmolecule) :: zmol
  real(wp),intent(out) :: mside

  integer :: k,l
  integer :: atm

  integer,allocatable :: path(:)

  mside = 0.0_wp

  allocate (path(zmol%nat),source=0)
  l = 0
  call recside(zmol,i,j,path,l)

  do k = 1,zmol%nat
    if (path(k) .eq. 0) cycle
    l = path(k)
    atm = zmol%zat(l)%atype
    mside = mside+ams(atm)
  end do

  deallocate (path)
  return
end subroutine getsideweight

!=======================================================================!
!  count the number of rings in the system
!  (just a call on the same "newgetrings" routine,
!   but some of the data has to be reset)
!=======================================================================!
subroutine countrings(zmol,nrings)
  use iso_fortran_env,only:wp => real64
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer,intent(out) :: nrings
  integer :: i
  nrings = 0
  call newgetrings(zmol,.false.)
  nrings = zmol%nri
  do i = 1,zmol%nat
    zmol%zat(i)%ring = .false.
  end do
  return
end subroutine countrings

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!================================================================================================!
! new "getring" routine.
! requires a FULLY SET UP TOPOLOGY in the "zmol" datatype
! Utilizes graph theory to detect rings via Dijkstra's algorithm.
!
!   DO. NOT. TOUCH.
!
!================================================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
subroutine newgetrings(zmol,verbose)
!> newest implementation of the getrings routine
!> that employs clean rewrites of graph and shortest-path routines
  use crest_parameters
  use zdata
  use adjacency
  implicit none
  !> INPUT/OUTPUT
  type(zmolecule),intent(inout) :: zmol
  logical,intent(in) :: verbose
  !> LOCAL
  type(zring) :: zri
  integer :: V,nrings,nmem
  integer,allocatable :: A(:,:) !> adjacency matrix
  logical,allocatable :: rings(:,:)
  logical,allocatable :: ringtracker(:,:)
  integer,allocatable :: tmppath(:)
  logical,allocatable :: tmpmem(:)
  logical :: duplicate
  integer :: i,j,k

  !>-- some printout
  if (verbose) then
    write (stdout,'(1x,a,1x,i0,1x,a)') 'ring analysis on',zmol%nat,'atoms:'
  end if

  !>-- reset number of rings
  nrings = 0

  !>-- atoms as vertices
  V = zmol%nat
  allocate (A(V,V),source=0)
  allocate (rings(V,V),ringtracker(V,V),source=.false.)

  !>-- get the adjacency matrix
  call zmol%adjacency(A)

  !>-- check possible rings
  call check_rings_min(V,A,rings)

  allocate (tmppath(V),source=0)
  allocate (tmpmem(V),source=.false.)
  !>-- sweeps
  do i = 1,V
    do j = 1,i-1
      nmem = 0
      !>-- if a possible ring was detected, get it
      if (rings(j,i)) then
        call get_ring_min(V,A,i,j,tmppath,nmem)
      end if
      !>-- if we were successful and now the number or ring members
      if (nmem > 0) then
        tmpmem(:) = .false.
        do k = 1,V
          if (tmppath(k) > 0) tmpmem(tmppath(k)) = .true.
        end do
        !>-- add the first ring
        if (nrings == 0) then
          nrings = nrings+1
          ringtracker(:,nrings) = tmpmem(:)
          !>-- or check duplicates
        else
          duplicate = .false.
          do k = 1,nrings
            if (all(ringtracker(:,k).eqv.tmpmem(:))) then
              duplicate = .true.
              exit
            end if
          end do
          if (.not.duplicate) then
            nrings = nrings+1
            ringtracker(:,nrings) = tmpmem(:)
          end if
        end if
      end if
    end do !> j loop
  end do !> i loop
  deallocate (tmpmem,tmppath)

  if (verbose) then
    if (nrings > 0) then
      write (stdout,'(3x,i0,1x,a)') nrings,'unique rings detected'
    else
      write (stdout,'(3x,a)') 'system contains no rings'
    end if
  end if

  !>--- put rings into zmol
  if (allocated(zmol%zri)) deallocate (zmol%zri)
  zmol%nri = 0
  zmol%zat(:)%ring = .false.
  do i = 1,nrings
    call zri%deallocate !> clear space (if allocated)
    nmem = count(ringtracker(:,i))
    allocate (tmppath(nmem),source=0)
    k = 0
    do j = 1,V
      if (ringtracker(j,i)) then
        k = k+1
        tmppath(k) = j
        zmol%zat(i)%ring = .true.
      end if
    end do
    zri%rs = nmem
    call move_alloc(tmppath,zri%rlist)
    call zmol%addring(zri)
    if (allocated(tmppath)) deallocate (tmppath)
  end do

  if (verbose) then
    call zmol%prrings(stdout)
  end if

  deallocate (ringtracker,rings,A)
  return
end subroutine newgetrings

!===============================================================!
!  get all atoms of the side chain l of k and write them to path
!===============================================================!
recursive subroutine recside(zmol,k,l,path,j)
  use iso_fortran_env,only:wp => real64
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer,intent(in) :: k          !the starting atom
  integer,intent(in) :: l          !current atom
  integer,intent(inout) :: j       !postion counter
  integer,intent(inout) :: path(zmol%nat)  !entire path
  integer :: i,p

  if (k .eq. l) return
  if (any(path(:) == l)) return
  if (j+1 .gt. zmol%nat) return

  j = j+1
  path(j) = l

  if (zmol%zat(l)%nei .gt. 1) then
    do i = 1,zmol%zat(l)%nei
      p = zmol%zat(l)%ngh(i)
      call recside(zmol,k,p,path,j)
    end do
  end if

  return
end subroutine recside
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!================================================================================================!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!

!==========================================================!
! "commonring" is a function (return value .true./.false.)
! that is used to check is a given list of atoms is part
! of the same ring.
!
! On Input:  zmol - molecule and topology
!            n    - # of atoms to be checked
!            atms - array containing the atoms
!
! On Output: function value
!            r - the common ring (if function value = .true.)
!==========================================================!
logical function commonring(zmol,n,atms,r)
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer,intent(in) :: n
  integer,intent(in) :: atms(n)
  integer,intent(out) :: r
  integer :: i,j
  logical :: l1,l2
  commonring = .false.
  r = 0
  !--- first check if the atoms are in a ring at all.
  !    if not, we can return immediatly.
  l1 = .true.
  do i = 1,n
    l1 = l1.and.zmol%zat(atms(i))%ring
  end do
  if (.not.l1) return
  !--- loop over all rings in the molecule and check if
  !    the atoms are part of it
  OUTER: do i = 1,zmol%nri
    l2 = .true.
    INNER: do j = 1,n
      l2 = l2.and.(any(zmol%zri(i)%rlist(:) .eq. atms(j)))
      if (.not.l2) cycle OUTER
    end do INNER
    if (l2) then
      commonring = .true.
      r = i
    end if
  end do OUTER
  return
end function commonring

!==========================================================!
! "ringside" is a subroutine that is used to check
! in which direction of which atom a ring extends,
! starting from atom i.
!
! On Input:  zmol - molecule and topology
!            n    - starting atom
!            ring - which ring to look at
!
! On Output: atms - array containing the atoms
!
!==========================================================!
subroutine ringside(zmol,i,ring,atms)
  use zdata
  implicit none
  type(zmolecule) :: zmol
  integer,intent(in) :: i
  integer,intent(in) :: ring
  integer,intent(inout) :: atms(zmol%nat)
  integer :: j,k,n,l
  integer,allocatable :: path(:)
  n = zmol%nat
  allocate (path(n),source=0)
  !--- loop over neighbour list of i
  do j = 1,zmol%zat(i)%nei
    k = zmol%zat(i)%ngh(j)
    !--- if the neighbour and i are in the selected ring, do nothing
    if (any(zmol%zri(ring)%rlist .eq. k)) cycle
    !--- else, get everything attached to  the neighbour
    path = 0
    l = 0
    call recside(zmol,i,k,path,l)
    call mergearr(n,atms,n,path)  !merge path into atms
  end do
  deallocate (path)

  return
end subroutine ringside
