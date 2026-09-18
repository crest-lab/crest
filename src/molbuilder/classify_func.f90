!=============================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
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
!=============================================================================!

module molbuilder_classify_func
  use crest_parameters,only:wp,stdout
  use strucrd,only:coord,i2e
  use adjacency
  use canonical_mod
  use molbuilder_classify_type
  use quicksort_interface,only:qqsorti
  implicit none
  private

  public :: functional_group_classify !> try to determine some functional groups

!==============================================================================!
contains  !> MODULE PROCEDURES START HERE
!==============================================================================!

!> FUNCTIONAL GROUP CLASSIFICATION ROUTINES

  subroutine functional_group_classify(molc)
    implicit none
    type(coord_classify),intent(inout) :: molc
    type(functional_group) :: fg
    integer :: ii,jj,nn,nfunc,nfunc2
    logical :: success,updated,duplicate

    if (.not.allocated(molc%atinfo)) then
      call atinfo_classify(molc)
    end if

    !> basic functional groups from single atoms
    do ii = 1,molc%nat
      call functional_group_classify_simple(molc,ii,fg,success)
      if (success) then
        call molc%add(fg)
      end if
    end do

    updated = .true.
    nfunc = molc%nfuncs
    do while (updated)

      do ii = 1,nfunc
        if (molc%funcgroups(ii)%seeded) cycle
        call functional_group_classify_extended(molc,ii,fg,success)
        if (success) then
          !> check for duplicates, only add new ones
          duplicate = .false.
          nfunc2 = molc%nfuncs
          do jj = 1,nfunc2
            if (molc%funcgroups(jj)%natms .ne. fg%natms) cycle
            if (all(molc%funcgroups(jj)%ids .eq. fg%ids)) duplicate = .true.
          end do
          if (.not.duplicate) call molc%add(fg)
        end if
        molc%funcgroups(ii)%seeded = .true.
      end do

      nn = molc%nfuncs
      if (nn == nfunc) then
        updated = .false.
      else
        nfunc = nn
      end if
    end do

  end subroutine functional_group_classify

  subroutine functional_group_classify_simple(molc,ii,fg,success)
    implicit none
    type(coord_classify),intent(in) :: molc
    type(functional_group),intent(out) :: fg
    integer,intent(in) :: ii
    logical,intent(out) :: success
    integer :: jj,kk
    success = .false.
    call fg%clear()
    select case (trim(molc%atinfo(ii)))
    case ('CH3')
      fg%name = 'methyl'
      fg%natms = 4
      allocate (fg%ids(4),source=0)
      fg%ids(1) = ii
      kk = 1
      do jj = 1,molc%nat
        if (molc%bond(jj,ii) == 1) then
          if (molc%at(jj) == 1) then
            kk = kk+1
            fg%ids(kk) = jj
          else
            fg%attached_to = jj
          end if
        end if
      end do
      success = .true.

    case ('OH')
      fg%name = "hydroxy"
      fg%natms = 2
      allocate (fg%ids(2),source=0)
      fg%ids(1) = ii
      kk = 1
      do jj = 1,molc%nat
        if (molc%bond(jj,ii) == 1) then
          if (molc%at(jj) == 1) then
            kk = kk+1
            fg%ids(kk) = jj
          else
            fg%attached_to = jj
          end if
        end if
      end do
      success = .true.

    case ('SH')
      fg%name = "thiol"
      fg%natms = 2
      allocate (fg%ids(2),source=0)
      fg%ids(1) = ii
      kk = 1
      do jj = 1,molc%nat
        if (molc%bond(jj,ii) == 1) then
          if (molc%at(jj) == 1) then
            kk = kk+1
            fg%ids(kk) = jj
          else
            fg%attached_to = jj
          end if
        end if
      end do
      success = .true.

    case ('NH2','NR2','NR3')
      fg%name = "amine"
      if (trim(molc%atinfo(ii)) == 'NH2') fg%name = fg%name//' (1°)'
      if (trim(molc%atinfo(ii)) == 'NR2') fg%name = fg%name//' (2°)'
      if (trim(molc%atinfo(ii)) == 'NR3') fg%name = fg%name//' (3°)'
      fg%natms = 3
      allocate (fg%ids(3),source=0)
      fg%ids(1) = ii
      kk = 1
      do jj = 1,molc%nat
        if (molc%bond(jj,ii) == 1) then
          if (molc%at(jj) == 1) then
            kk = kk+1
            fg%ids(kk) = jj
          else
            !TODO: logic for NR2 and NR3
            fg%attached_to = jj
          end if
        end if
      end do
      success = .true.

    case ('o')
      fg%natms = 2
      allocate (fg%ids(2),source=0)
      fg%ids(1) = ii
      kk = 1
      do jj = 1,molc%nat
        if (molc%bond(jj,ii) == 1) then
          if (molc%at(jj) == 6) then
            kk = kk+1
            fg%ids(kk) = jj
            fg%attached_to = jj
          end if
        end if
      end do
      if (kk == 2) then
        fg%name = 'carbonyl'
      end if
      success = .true.

    case ('F','Cl','Br','I')
      fg%name = 'halide'
      fg%natms = 1
      allocate (fg%ids(1),source=ii)
      fg%attached_to = maxloc(molc%bond(:,ii),1)
      success = .true.

    end select
  end subroutine functional_group_classify_simple

  subroutine functional_group_classify_extended(molc,ii,fg,success)
    implicit none
    type(coord_classify),intent(inout) :: molc
    type(functional_group),intent(out) :: fg
    integer,intent(in) :: ii
    logical,intent(out) :: success
    integer :: jj,kk
    success = .false.
    call fg%clear()
    select case (trim(molc%funcgroups(ii)%name))
    case ('methyl')

      !> n-alkyl chains always start from methyl
      call check_alkyl(molc,ii,fg,success)
      if (success) then
        fg%natms = count(molc%lwork)
        allocate (fg%ids(fg%natms),source=0)
        kk = 0
        do jj = 1,molc%nat
          if (molc%lwork(jj)) then
            kk = kk+1
            fg%ids(kk) = jj
          end if
        end do
        call qqsorti(fg%ids,1,fg%natms)
      end if

    end select
  end subroutine functional_group_classify_extended

!=============================================================================!
!#############################################################################!
!=============================================================================!

!> routines to check for specific functional groups

  subroutine check_alkyl(molc,istart,fg,success)
    implicit none
    type(coord_classify),intent(inout) :: molc
    type(functional_group),intent(out) :: fg
    integer,intent(in) :: istart
    logical,intent(out) :: success
    logical :: contin
    integer :: ii,jj,nat,ati,atii

    success = .false.

    !> prepare
    nat = molc%nat
    if (.not.allocated(molc%lwork)) allocate (molc%lwork(nat))

    molc%lwork(:) = .false.
    ati = molc%funcgroups(istart)%ids(1)
    do ii = 1,molc%funcgroups(istart)%natms
      molc%lwork(molc%funcgroups(istart)%ids(ii)) = .true.
    end do
    !> walk
    contin = .true.
    atii = ati
    do while (contin)
      do ii = 1,nat
        if (molc%lwork(ii)) cycle !> skip alreaty iterated atoms
        if (molc%bond(ii,ati) == 1) then
          if (molc%at(ii) == 1) then
            molc%lwork(ii) = .true. !> H's simply set to true
          else if (molc%at(ii) == 6) then
            if (trim(molc%atinfo(ii)) == 'CH2') then
              !> continue chain
              fg%name = 'alkyl'
              success = .true. !> at the first occurence of CH2 we have at least ethyl
              atii = ii
              molc%lwork(ii) = .true.
            else if (trim(molc%atinfo(ii)) == 'CH3') then
              !> continue chain
              fg%name = 'alkane'
              success = .true. !> at the first occurence of CH2 we have at least ethyl
              atii = ii
              molc%lwork(ii) = .true.
            else
              !> terminate chain on non-CH2
              if (success) fg%attached_to = ii
              contin = .false.
            end if
          else
            !> terminate chain for hetero atoms
            if (success) fg%attached_to = ii
            contin = .false.
          end if
        end if
      end do
      if (atii == ati) contin = .false.
      ati = atii
    end do

  end subroutine check_alkyl

!=============================================================================!
!#############################################################################!
!=============================================================================!
end module molbuilder_classify_func

