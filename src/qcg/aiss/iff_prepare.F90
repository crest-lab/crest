!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2026 Christoph Plett, Philipp Pracht
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
! Ported from xtb src/iff/iff_prepare.f90 and src/local.f90
! (github.com/grimme-lab/xtb), LGPL-3.0-or-later.
!================================================================================!

!> Electronic-structure precomputation (LMOs, CT data) of the xTB-IFF fragments
module aiss_iff_prepare
  use crest_parameters,only:wp,stdout,autoev
  use strucrd,only:coord
  use aiss_param,only:iff_data
  use mo_localize,only:local_iff,Dints,cao2saop
#ifdef WITH_TBLITE
  use tblite_api,only:tblite_data,tblite_setup,tblite_singlepoint,xtblvl
  use mctc_io,only:structure_type,new
  use tblite_integral_multipole,only:get_multipole_integrals
  use tblite_basis_type,only:get_cutoff
#endif
#ifdef WITH_GFN0
  use gfn0_api,only:gfn0_data,gfn0_setup,gfn0_addsettings,gfn0_sp
#endif
  implicit none
  private

  public :: iff_precomp
  public :: lmo_basic_atoms

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine iff_precomp(mol,lmolvl,etemp,d,imol,etot,io)
    !***************************************************************
    !* Gas-phase singlepoint of fragment imol (1 or 2) and IFF-mode
    !* localization; fills charges, LMO centers/types, HOMO/LUMO
    !* data and dipole of that fragment in d (allocated by caller).
    !* lmolvl : 'gfn2'/'gfn1' (tblite) or 'gfn0'; etemp in K
    !***************************************************************
    type(coord),intent(in) :: mol
    character(len=*),intent(in) :: lmolvl
    real(wp),intent(in) :: etemp
    type(iff_data),intent(inout) :: d
    integer,intent(in) :: imol
    real(wp),intent(out) :: etot
    integer,intent(out) :: io

    integer :: nat,nlmo
    integer,allocatable  :: lmo(:)
    real(wp),allocatable :: q(:),rlmo(:,:),qct(:,:)
    real(wp) :: ehomo,elumo,dipol

    nat = mol%nat
    allocate (q(nat),qct(nat,2),source=0.0_wp)
    allocate (rlmo(4,10*nat),source=0.0_wp)
    allocate (lmo(10*nat),source=0)
    etot = 0.0_wp
    nlmo = 0
    ehomo = 0.0_wp
    elumo = 0.0_wp
    dipol = 0.0_wp

    select case (lmolvl)
    case ('gfn2','gfn1')
      call precomp_tblite(mol,lmolvl,etemp,etot,q,nlmo,lmo,rlmo,qct,ehomo,elumo,dipol,io)
    case ('gfn0')
      call precomp_gfn0(mol,etemp,etot,q,nlmo,lmo,rlmo,qct,ehomo,elumo,dipol,io)
    case default
      write (stdout,'(a,a)') 'aISS: unknown level for the IFF precomputation: ',lmolvl
      io = 1
    end select
    if (io /= 0) return

    select case (imol)
    case (1)
      d%at1 = mol%at
      d%xyz1 = mol%xyz
      d%q1 = q
      d%nlmo1 = nlmo
      d%lmo1 = lmo
      d%rlmo1 = rlmo
      d%qct1 = qct
    case (2)
      d%at2 = mol%at
      d%xyz2 = mol%xyz
      d%q2 = q
      d%nlmo2 = nlmo
      d%lmo2 = lmo
      d%rlmo2 = rlmo
      d%qct2 = qct
    end select
    d%ehomo(imol) = ehomo
    d%elumo(imol) = elumo
    d%dipol(imol) = dipol
    d%etot(imol) = etot
  end subroutine iff_precomp

!========================================================================================!

  subroutine lmo_basic_atoms(mol,lmolvl,etemp,basat,io)
    !***************************************************************
    !* Atoms with a lone pair or pi LMO (xtb --lmo), e.g. as
    !* targets for H shifts; mol%chrg/uhf are used
    !***************************************************************
    type(coord),intent(in) :: mol
    character(len=*),intent(in) :: lmolvl
    real(wp),intent(in) :: etemp
    logical,intent(out) :: basat(mol%nat)
    integer,intent(out) :: io
    integer :: nlmo
    integer,allocatable  :: lmo(:)
    real(wp),allocatable :: q(:),rlmo(:,:),qct(:,:)
    real(wp) :: etot,ehomo,elumo,dipol

    basat = .false.
    allocate (q(mol%nat),qct(mol%nat,2),source=0.0_wp)
    allocate (rlmo(4,10*mol%nat),source=0.0_wp)
    allocate (lmo(10*mol%nat),source=0)
    select case (lmolvl)
    case ('gfn2','gfn1')
      call precomp_tblite(mol,lmolvl,etemp,etot,q,nlmo,lmo,rlmo,qct,ehomo,elumo,dipol,io,basat)
    case ('gfn0')
      call precomp_gfn0(mol,etemp,etot,q,nlmo,lmo,rlmo,qct,ehomo,elumo,dipol,io,basat)
    case default
      write (stdout,'(a,a)') 'unknown level for the LMO calculation: ',lmolvl
      io = 1
    end select
  end subroutine lmo_basic_atoms

!========================================================================================!

  subroutine precomp_tblite(mol,lmolvl,etemp,etot,q,nlmo,lmo,rlmo,qct,ehomo,elumo,dipol,io,basat)
    !***************************************************************
    !* GFN2/GFN1-xTB singlepoint via tblite, then local_iff
    !***************************************************************
    type(coord),intent(in) :: mol
    character(len=*),intent(in) :: lmolvl
    real(wp),intent(in) :: etemp
    real(wp),intent(out) :: etot,q(:),rlmo(:,:),qct(:,:),ehomo,elumo,dipol
    integer,intent(out) :: nlmo,lmo(:),io
    logical,intent(out),optional :: basat(:)
#ifdef WITH_TBLITE
    type(tblite_data),allocatable :: tb
    type(structure_type) :: mctcmol
    real(wp),allocatable :: grd(:,:),s(:,:),dpint(:,:,:),qpint(:,:,:),dip(:,:),z(:)
    real(wp) :: trans(3,1)
    integer :: nao,lvl,ihomo,i,j,k

    lvl = xtblvl%gfn2
    if (lmolvl == 'gfn1') lvl = xtblvl%gfn1
    allocate (tb)
    tb%ctx%verbosity = 0
    call tblite_setup(mol,mol%chrg,mol%uhf,lvl,etemp,tb,.false.)
    allocate (grd(3,mol%nat))
    call tblite_singlepoint(mol,mol%chrg,mol%uhf,tb,etot,grd,io)
    if (io /= 0) then
      write (stdout,'(a)') 'aISS: tblite singlepoint of the IFF precomputation failed'
      return
    end if

!>-- overlap and dipole integrals; tblite dipoles are <mu|r-R_B|nu> with B the atom of nu
    nao = tb%calc%bas%nao
    allocate (s(nao,nao),dpint(3,nao,nao),qpint(6,nao,nao))
    call new(mctcmol,mol%at,mol%xyz,charge=real(mol%chrg,wp),uhf=mol%uhf)
    trans = 0.0_wp
    call get_multipole_integrals(mctcmol,trans,get_cutoff(tb%calc%bas),tb%calc%bas, &
    &                            s,dpint,qpint)
    allocate (dip(nao*(nao+1)/2,3))
    k = 0
    do i = 1,nao
      do j = 1,i
        k = k+1
        dip(k,1:3) = -(dpint(1:3,j,i)+mol%xyz(1:3,tb%calc%bas%ao2at(i))*s(j,i))
      end do
    end do
    deallocate (dpint,qpint)

    call mol%get_z(z)
    ihomo = nint(tb%wfn%nel(1))
    q(:) = tb%wfn%qat(:,1)
    call local_iff(mol%nat,mol%at,nao,ihomo,mol%xyz,z,tb%wfn%focc(:,1),s, &
    &              tb%wfn%density(:,:,1),tb%wfn%coeff(:,:,1),tb%wfn%emo(:,1), &
    &              tb%calc%bas%ao2at,dip,nlmo,lmo,rlmo,qct,ehomo,elumo,dipol,basat)
    io = 0
#else
    etot = 0.0_wp
    nlmo = 0
    dipol = 0.0_wp
    ehomo = 0.0_wp
    elumo = 0.0_wp
    write (stdout,'(a)') 'aISS: the '//lmolvl//' IFF precomputation requires tblite (compiled without)'
    io = 1
#endif
  end subroutine precomp_tblite

!========================================================================================!

  subroutine precomp_gfn0(mol,etemp,etot,q,nlmo,lmo,rlmo,qct,ehomo,elumo,dipol,io,basat)
    !***************************************************************
    !* GFN0-xTB singlepoint via the gfn0 library, then local_iff;
    !* q are the EEQ charges of the non-SCC model, emo eV -> Eh
    !***************************************************************
    type(coord),intent(in) :: mol
    real(wp),intent(in) :: etemp
    real(wp),intent(out) :: etot,q(:),rlmo(:,:),qct(:,:),ehomo,elumo,dipol
    integer,intent(out) :: nlmo,lmo(:),io
    logical,intent(out),optional :: basat(:)
#ifdef WITH_GFN0
    type(gfn0_data),allocatable :: g0
    real(wp),allocatable :: grd(:,:),dip(:,:),z(:),emo(:)
    integer :: nbf,nao

    allocate (g0)
    call gfn0_setup(mol,mol%chrg,mol%uhf,g0)
    call gfn0_addsettings(mol,g0,loadwbo=.true.)
    call g0%wfn%refresh_occu(g0%wfn%nel,mol%uhf)
    call gfn0_addsettings(mol,g0,etemp=etemp)
    allocate (grd(3,mol%nat))
    call gfn0_sp(mol,mol%chrg,mol%uhf,g0,etot,grd,io)
    if (io /= 0) then
      write (stdout,'(a)') 'aISS: GFN0 singlepoint of the IFF precomputation failed'
      return
    end if

    nbf = g0%basis%nbf
    nao = g0%basis%nao
    allocate (dip(nbf*(nbf+1)/2,3))
    call Dints(mol%nat,nbf,mol%xyz,dip(:,1),dip(:,2),dip(:,3), &
    &          g0%basis%aoat,g0%basis%nprim,g0%basis%alp,g0%basis%lao,g0%basis%cont)
    call cao2saop(nbf,nao,dip(:,1),g0%basis%lao)
    call cao2saop(nbf,nao,dip(:,2),g0%basis%lao)
    call cao2saop(nbf,nao,dip(:,3),g0%basis%lao)

    call mol%get_z(z)
    emo = g0%wfn%emo/autoev
    q(:) = g0%wfn%q(:)
    call local_iff(mol%nat,mol%at,nao,g0%wfn%ihomo,mol%xyz,z,g0%wfn%focc,g0%wfn%S, &
    &              g0%wfn%P,g0%wfn%C,emo,g0%basis%aoat2,dip(1:nao*(nao+1)/2,:), &
    &              nlmo,lmo,rlmo,qct,ehomo,elumo,dipol,basat)
    io = 0
#else
    etot = 0.0_wp
    nlmo = 0
    dipol = 0.0_wp
    ehomo = 0.0_wp
    elumo = 0.0_wp
    write (stdout,'(a)') 'aISS: the gfn0 IFF precomputation requires the gfn0 library (compiled without)'
    io = 1
#endif
  end subroutine precomp_gfn0

!========================================================================================!
end module aiss_iff_prepare
