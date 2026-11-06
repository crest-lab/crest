!================================================================================!
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
!================================================================================!

module molecule_io
  use iso_c_binding
  use molecule_parameters
  use molecule_type_components
!> simple geomerty and vector operations
  use geo
!> element symbols
  use crest_cn_module,only:calculate_cn
  implicit none

! ──────────────────────────────────────────────────────────────────────────────
!>--- private module variables and parameters
  private

!>--- private utility subroutines
  private :: upperCase,lowerCase
  private :: convertlable,fextension

! ──────────────────────────────────────────────────────────────────────────────
!>--- public subroutines
  public :: i2e          !> function to convert atomic number to element symbol
  public :: asym         !> "
  interface asym         !> "
    module procedure i2e !> "
  end interface asym
  public :: e2i          !> function to convert element symbol into atomic number
  public :: grepenergy
  public :: checkcoordtype

  public :: rdnat       !-- procedure to read number of atoms Nat
  public :: rdcoord     !-- read an input file, determine format automatically
  public :: rdxmol      !-- read a file in the Xmol (.xyz) format specifically
  public :: rdxmolselec !-- read only a certain structure in Xmol file
  public :: rdPDB
  public :: read_extxyz_frame
  public :: get_at_from_ext,get_xyz_from_ext,get_grad_from_ext

  !>--- write a TM coord file
  public :: wrc0
  interface wrc0
    module procedure wrc0_file
    module procedure wrc0_channel
  end interface wrc0
  public :: wrcoord
  interface wrcoord
    module procedure wrc0_file
    module procedure wrc0_channel
  end interface wrcoord

  !>--- write a XYZ coord file
  public :: wrxyz
  interface wrxyz
    module procedure wrxyz_file
    module procedure wrxyz_file_mask
    module procedure wrxyz_channel_energy
    module procedure wrxyz_channel
  end interface wrxyz

  !>--- write a sdf molfile
  public :: wrsdf
  interface wrsdf
    module procedure wrsdf_channel
  end interface wrsdf
  public :: wrsdfV2000
  interface wrsdfV2000
    module procedure wrsdf_channel
  end interface wrsdfV2000
  interface wrsdfV3000
    module procedure wrsdfV3000_channel
  end interface wrsdfV3000
  public :: wrsdfV3000

  public :: coordline
  public :: sgrep
  public :: get_atlist
  public :: sumform

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════
!  ROUTINES FOR READING SINGLE STRUCTURES (COORDS)
! ──────────────────────────────────────────────────────────────────────────────

  subroutine checkcoordtype(fname,typint)
!*****************************************************
!* subroutine checkcoordtype                         *
!* try to identify the filetype of the coord type.   *
!* first based on file extension, if that fails by   *
!* a keyword within the file.                        *
!*****************************************************
    implicit none
    character(len=*) :: fname
    integer,intent(out) :: typint
    typint = coordtype%unknown
    !-- check file extension first
    select case (fextension(fname))
    case ('.coord','.COORD')
      typint = coordtype%turbomole
    case ('.xyz','.XYZ', &
        & '.trj','.TRJ','.sorted')
      typint = coordtype%xyz
      if (sgrep(fname,'Properties=',casesensitive=.false.)) then
        typint = coordtype%extxyz
      end if
    case ('.extxyz','.EXTXYZ')
      typint = coordtype%extxyz
    case ('.sd','.sdf','.SDF','.mol','.MOL')
      typint = coordtype%sdf
      if (sgrep(fname,'V2000')) then
        typint = coordtype%sdfV2000
      end if
      if (sgrep(fname,'V3000')) then
        typint = coordtype%sdfV3000
      end if
    case ('.pdb','.PDB')
      typint = coordtype%PDB
    case default
      typint = coordtype%unknown
    end select

    if (typint .ne. coordtype%unknown) return !-- file extension was recognized
    !-- grep for keywords otherwise
    if (sgrep(fname,'$coord')) then
      typint = coordtype%turbomole
    else !--no match found
      typint = coordtype%unknown
    end if
    return
  end subroutine checkcoordtype

! ──────────────────────────────────────────────────────────────────────────────

  subroutine rdnat(fname,nat,ftype)
!*******************************************************************
!* subroutine rdnat                                                *
!* read number of atoms "nat" form file                            *
!*                                                                 *
!* On Input: fname  - name of the coord file                       *
!*           ftype  - (OPTIONAL) format of the input coord file    *
!*                    if ftype is not present, it is determined    *
!* On Output: nat   - number of atoms                              *
!*******************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(out) :: nat
    integer,optional :: ftype
    integer :: ftypedum
    integer :: ich,i,j,io,k
    logical :: ex
    character(len=256) :: atmp
    nat = 0
    inquire (file=fname,exist=ex)
    if (.not.ex) then
      write (stdout,'(a)') '**ERROR** could not find coord file '//trim(fname)
      call exit(1)
    end if
    if (present(ftype)) then
      ftypedum = ftype
    else
      call checkcoordtype(fname,ftypedum)
    end if
    open (newunit=ich,file=fname)
    select case (ftypedum)

    case (coordtype%xyz,coordtype%extxyz)  !--- *.xyz files
      read (ich,*,iostat=io) nat

    case (coordtype%turbomole)      !--- TM coord file
      do
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit
        atmp = adjustl(atmp)
        if (index(atmp,"$coord") .eq. 1) exit
      end do
      do
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit
        atmp = adjustl(atmp)
        if (atmp(1:1) == '$') exit
        nat = nat+1
      end do

    case (coordtype%sdfV2000)      !--- sdf V2000 (or *.mol) file
      do i = 1,3 !-- first three comment lines
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit
      end do
      read (ich,'(a)',iostat=io) atmp
      if (index(atmp,'V2000') .ne. 0) then
        read (atmp,'(i3)') nat !- first argument is nat
      end if

    case (coordtype%sdfV3000)      !--- sdf V3000 file
      do
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit
        if ((index(atmp,'V30') .ne. 0).and. &
        &  (index(atmp,'COUNTS') .ne. 0)) then
          j = index(atmp,'COUNTS')+6
          k = len_trim(atmp)
          atmp = atmp(j:k)
          atmp = adjustl(atmp)
          read (atmp,*) nat
        end if
      end do

    case (coordtype%PDB)      !--- pdb file
      nat = 0
      do
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit
        if ((index(atmp,'ATOM') .eq. 1).or. &
        &  (index(atmp,'HETATM') .eq. 1)) then
          nat = nat+1
        end if
      end do

    case default
      continue
    end select
    close (ich)
    return
  end subroutine rdnat

! ──────────────────────────────────────────────────────────────────────────────

  subroutine rdcoord(fname,nat,at,xyz,energy,ftype)
!*****************************************************************
!* subroutine rdcoord                                            *
!* read in a structure. The format is determined automatically   *
!*                                                               *
!* On Input: fname  - name of the coord file                     *
!*           nat    - number of atoms                            *
!*           ftype  - coord file type (optional)                 *
!*                                                               *
!* On Output: at   - atom number as integer                      *
!*            xyz  - coordinates (always in Bohr)                *
!*            energy - (OPTIONAL) if present, try to get energy  *
!*                      mainly from xyz files                    *
!*****************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    real(wp),optional :: energy
    integer,intent(in),optional :: ftype
    character(len=256) :: atmp
    integer :: ftypedum
    type(pdbdata) :: pdbdummy

    if (present(ftype)) then
      ftypedum = ftype
    else
      call checkcoordtype(fname,ftypedum)
    end if

    select case (ftypedum)
    case (coordtype%turbomole)  !-- TM coord file, always retruns coords in Bohr
      call rdtmcoord(fname,nat,at,xyz)

    case (coordtype%xyz)     !-- XYZ file, is Angström, needs conversion
      if (present(energy)) then
        call rdxmol(fname,nat,at,xyz,atmp)
        energy = grepenergy(atmp)
      else
        call rdxmol(fname,nat,at,xyz)
      end if
      xyz = xyz/bohr

    case (coordtype%sdfV2000)      !-- SDF/MOL V2000 file, also Angström
      call rdsdf(fname,nat,at,xyz)
      xyz = xyz/bohr

    case (coordtype%sdfV3000)     !-- SDF V3000 file, Angström
      call rdsdfV3000(fname,nat,at,xyz)
      xyz = xyz/bohr

    case (coordtype%PDB)          !-- PDB file, Angström
      call rdPDB(fname,nat,at,xyz,pdbdummy)
      xyz = xyz/bohr
      call pdbdummy%deallocate()

    case default
      continue
    end select

    return
  end subroutine rdcoord

! ──────────────────────────────────────────────────────────────────────────────

  subroutine rdtmcoord(fname,nat,at,xyz)
!**************************************************************
!* subroutine rdtmcoord                                       *
!* read a struncture in the TM coord style.                   *
!*                                                            *
!* On Input: fname  - name of the coord file                  *
!*           nat    - number of atoms                         *
!*                                                            *
!* On Output: at   - atom number as integer                   *
!*            xyz  - coordinates (always in Bohr)             *
!**************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    character(len=6) :: sym
    integer :: ich,io,i
    real(wp) :: convert
    character(len=256) :: atmp
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      if (index(atmp,"$coord") .eq. 1) exit
    end do
    if (index(atmp,'ang') .ne. 0) then
      !> coord files allow explicit specification in Angström
      convert = aatoau
    else
      convert = 1.0_wp
    end if
    do i = 1,nat
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      if (atmp(1:1) == '$') exit
      call coordline(atmp,sym,xyz(1:3,i),io)
      if (io < 0) then
        write (*,*) 'error while reading coord line. EOF'
        exit
      end if
      at(i) = e2i(sym)
    end do
    close (ich)
    xyz = xyz*convert
    return
  end subroutine rdtmcoord

! ──────────────────────────────────────────────────────────────────────────────

  subroutine rdxmol(fname,nat,at,xyz,comment)
!***************************************************************
!* subroutine rdxmol                                           *
!* read a struncture in the *.xyz (Xmol) style.                *
!* The commentary (second) line is ignored                     *
!*                                                             *
!* On Input: fname  - name of the coord file                   *
!*           nat    - number of atoms                          *
!*                                                             *
!* On Output: at   - atom number as integer                    *
!*            xyz  - coordinates (in Angström)                 *
!*            comment - (OPTIONAL) commentary line of the file *
!***************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    character(len=*),optional :: comment
    character(len=6) :: sym
    integer :: ich,io,i
    integer :: dum
    character(len=256) :: atmp
    open (newunit=ich,file=fname)
    read (ich,*,iostat=io) dum
    if (nat .ne. dum) then
      write (stdout,'(a)') '**ERROR** Mismatch in expected atom number for file '//trim(fname)
      write (stdout,'(a,i0,a,i0)') '          Expected ',nat,' got ',dum
      call exit(1)
    end if
    read (ich,'(a)') atmp !--commentary line
    if (present(comment)) comment = trim(adjustl(atmp))
    do i = 1,nat
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      call coordline(atmp,sym,xyz(1:3,i),io)
      if (io < 0) then
        write (stdout,'(a)') '**ERROR** Unexpected EOF while reading file '//trim(fname)
        exit
      end if
      at(i) = e2i(sym)
    end do
    close (ich)
    return
  end subroutine rdxmol

! ──────────────────────────────────────────────────────────────────────────────

  subroutine rdsdf(fname,nat,at,xyz,comment)
!***************************************************************
!* subroutine rdsdf                                            *
!* read a struncture in the .sdf/.mol V2000 style.             *
!*                                                             *
!* On Input: fname  - name of the coord file                   *
!*           nat    - number of atoms                          *
!*                                                             *
!* On Output: at   - atom number as integer                    *
!*            xyz  - coordinates (in Angström)                 *
!*            comment - (OPTIONAL) commentary line of the file *
!***************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    character(len=*),optional :: comment
    character(len=6) :: sym
    integer :: ich,io,i
    integer :: dum
    character(len=256) :: atmp
    open (newunit=ich,file=fname)
    read (ich,'(a)',iostat=io) atmp
    read (ich,'(a)',iostat=io) atmp
    read (ich,'(a)',iostat=io) atmp
    if (present(comment)) comment = trim(adjustl(atmp))
    read (ich,'(i3)',iostat=io) dum
    if (nat .ne. dum) then
      write (stdout,'(a)') '**ERROR** Mismatch in expected atom number for file '//trim(fname)
      write (stdout,'(a,i0,a,i0)') '          Expected ',nat,' got ',dum
      call exit(1)
    end if
    do i = 1,nat
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      call coordline(atmp,sym,xyz(1:3,i),io)
      if (io < 0) then
        write (*,*) 'error while reading coord line. EOF'
        exit
      end if
      at(i) = e2i(sym)
    end do
    close (ich)
    return
  end subroutine rdsdf

! ──────────────────────────────────────────────────────────────────────────────

  subroutine rdsdfV3000(fname,nat,at,xyz,comment)
!***************************************************************
!* subroutine rdsdfV3000                                       *
!* read a struncture in the .sdf/.mol V3000 style.             *
!*                                                             *
!* On Input: fname  - name of the coord file                   *
!*           nat    - number of atoms                          *
!*                                                             *
!* On Output: at   - atom number as integer                    *
!*            xyz  - coordinates (in Angström)                 *
!*            comment - (OPTIONAL) commentary line of the file *
!***************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    character(len=*),optional :: comment
    character(len=6) :: sym
    integer :: ich,io,i,j,k,l
    integer :: dum
    character(len=256) :: atmp
    character(len=32) :: btmp
    open (newunit=ich,file=fname)
    read (ich,'(a)',iostat=io) atmp
    read (ich,'(a)',iostat=io) atmp
    read (ich,'(a)',iostat=io) atmp
    if (present(comment)) comment = trim(adjustl(atmp))
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      if ((index(atmp,'V30') .ne. 0).and. &
      &  (index(atmp,'COUNTS') .ne. 0)) then
        j = index(atmp,'COUNTS')+6
        k = len_trim(atmp)
        atmp = atmp(j:k)
        atmp = adjustl(atmp)
        read (atmp,*) dum
      end if
      if ((index(atmp,'V30') .ne. 0).and. &
      &  (index(atmp,'ATOM') .ne. 0)) then
        exit
      end if
    end do
    if (nat .ne. dum) then
      write (stdout,'(a)') '**ERROR** Mismatch in expected atom number for file '//trim(fname)
      write (stdout,'(a,i0,a,i0)') '          Expected ',nat,' got ',dum
      call exit(1)
    end if
    do i = 1,nat
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      write (btmp,'(i0)') i
      l = len_trim(btmp)+1
      j = index(atmp,'V30')+3
      k = len_trim(atmp)
      atmp = atmp(j:k)
      atmp = adjustl(atmp)
      atmp = atmp(l:k)
      call coordline(atmp,sym,xyz(1:3,i),io)
      if (io < 0) then
        write (*,*) 'error while reading coord line. EOF'
        exit
      end if
      at(i) = e2i(sym)
    end do
    close (ich)
    return
  end subroutine rdsdfV3000

! ──────────────────────────────────────────────────────────────────────────────

  subroutine rdPDB(fname,nat,at,xyz,pdb)
!***********************************************
!* subroutine rdPDB                            *
!* read a struncture in the .PDB style.        *
!*                                             *
!* On Input: fname  - name of the coord file   *
!*           nat    - number of atoms          *
!*                                             *
!* On Output: at   - atom number as integer    *
!*            xyz  - coordinates (in Angström) *
!*            pdb  - pdbdata object            *
!***********************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    type(pdbdata) :: pdb
    character(len=2) :: sym
    integer :: ich,io,i,j,k
    character(len=256) :: atmp
    character(len=6) :: dum1
    character(len=1) :: dum2,dum3,pdbgp
    character(len=3) :: pdbas
    character(len=2) :: dum4
    character(len=4) :: pdbat
    real(wp) :: r1,r2
    call pdb%allocate(nat)
    open (newunit=ich,file=fname)
    k = 0
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      if ((index(atmp,'ATOM') .eq. 1).or. &
      &  (index(atmp,'HETATM') .eq. 1)) then
        k = k+1
        read (atmp,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)') &
        &  dum1,i,pdbat,dum2,pdbas,pdbgp,j,dum3,xyz(1:3,k),r1,r2,sym,dum4
        at(k) = e2i(sym)
        pdb%pdbat(k) = pdbat
        pdb%pdbas(k) = pdbas
        pdb%pdbgrp(k) = pdbgp
        pdb%pdbfrag(k) = j
        pdb%pdbocc(k) = r1
        pdb%pdbtf(k) = r2
      end if
    end do
    close (ich)
    return
  end subroutine rdPDB

! ──────────────────────────────────────────────────────────────────────────────

  subroutine rdxmolselec(fname,m,nat,at,xyz,comment)
!*******************************************************************
!* subroutine rdxmolselec                                          *
!* Read a file with multiple structures in the *.xyz (Xmol) style. *
!* Picks one structure.                                            *
!* The commentary (second) line is ignored                         *
!*                                                                 *
!* On Input: fname  - name of the coord file                       *
!*           m      - position of the desired structure            *
!*           nat    - number of atoms                              *
!*                                                                 *
!* On Output: at   - atom number as integer                        *
!*            xyz  - coordinates (in Bohr)                         *
!*******************************************************************
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat,m
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    character(len=*),optional :: comment
    character(len=6) :: sym
    integer :: ich,io,i,j
    integer :: dum
    character(len=256) :: atmp

    open (newunit=ich,file=fname)

    do j = 1,m
      read (ich,*,iostat=io) dum
      if (nat .ne. dum) then
        write (stdout,'(a)') '**ERROR** Mismatch in expected atom number for file '//trim(fname)
        write (stdout,'(a,i0,a,i0)') '          Expected ',nat,' got ',dum
        call exit(1)
      end if
      read (ich,'(a)') atmp !--commentary line
      if (present(comment)) comment = trim(adjustl(atmp))
      do i = 1,nat
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit
        atmp = adjustl(atmp)
        call coordline(atmp,sym,xyz(1:3,i),io)
        if (io < 0) then
          write (*,*) 'error while reading coord line. EOF'
          exit
        end if
        at(i) = e2i(sym)
      end do
    end do
    close (ich)
    xyz = xyz/bohr
    return
  end subroutine rdxmolselec

!=========================================================================================!
!=========================================================================================!
!  3. ROUTINES FOR WRITING STRUCTURES AND CONVERTING THEM
!=========================================================================================!
!=========================================================================================!

!============================================================!
! subroutine wrc0_file
! this is the typical quick write routine for TM coord files
! version for writing directly to a new file
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Bohr)
!
! On Output: file written to "fname"
!============================================================!
  subroutine wrc0_file(fname,nat,at,xyz)
    implicit none
    character(len=*) :: fname
    integer :: nat
    integer :: at(nat)
    real(wp) ::  xyz(3,nat)
    integer :: i,j,k,ich,io
    logical :: ex
    open (newunit=ich,file=fname,status='replace')
    write (ich,'(''$coord'')')
    do j = 1,nat
      write (ich,'(3F24.12,5x,a2)') xyz(1:3,j),i2e(at(j),'lc')
    end do
    write (ich,'(''$end'')')
    close (ich)
    return
  end subroutine wrc0_file

!============================================================!
! subroutine wrc0_channel
! this is the typical quick write routine for TM coord files
! version for writing to a output channel
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Bohr)
!
! On Output: file written to "fname"
!============================================================!
  subroutine wrc0_channel(ch,nat,at,xyz)
    implicit none
    integer :: ch
    integer :: nat
    integer :: at(nat)
    real(wp) ::  xyz(3,nat)
    integer :: i,j,k,ich,io
    logical :: ex
    write (ch,'(''$coord'')')
    do j = 1,nat
      write (ch,'(3F24.12,5x,a2)') xyz(1:3,j),i2e(at(j),'lc')
    end do
    write (ch,'(''$end'')')
    return
  end subroutine wrc0_channel

!============================================================!
! subroutine wrxyz_file
! this is the typical quick write routine for TM coord files
! version for writing directly to a new file
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           comment - (OPTIONAL) comment line
!
! On Output: file written to "fname"
!============================================================!
  subroutine wrxyz_file(fname,nat,at,xyz,comment)
    implicit none
    character(len=*) :: fname
    integer :: nat
    integer :: at(nat)
    real(wp) ::  xyz(3,nat)
    character(len=*),optional :: comment
    integer :: i,j,k,ich,io
    logical :: ex
    open (newunit=ich,file=fname,status='replace')
    write (ich,'(2x,i0)') nat
    if (present(comment)) then
      write (ich,'(a)') trim(comment)
    else
      write (ich,*)
    end if
    do j = 1,nat
      write (ich,'(1x,a2,1x,3f20.10)') i2e(at(j),'nc'),xyz(1:3,j)
    end do
    close (ich)
    return
  end subroutine wrxyz_file

!============================================================!
! subroutine wrxyz_file_mask
! this is the typical quick write routine for TM coord files
! version for writing directly to a new file
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           mask - a mask to determine to write which atoms
!           comment - (OPTIONAL) comment line
!
! On Output: file written to "fname"
!============================================================!
  subroutine wrxyz_file_mask(fname,nat,at,xyz,mask,comment)
    implicit none
    character(len=*) :: fname
    integer :: nat
    integer :: at(nat)
    real(wp) ::  xyz(3,nat)
    logical :: mask(nat)
    integer :: maskednat
    character(len=*),optional :: comment
    integer :: i,j,k,ich,io
    logical :: ex
    open (newunit=ich,file=fname,status='replace')
    maskednat = count(mask(:))
    write (ich,'(2x,i0)') maskednat
    if (present(comment)) then
      write (ich,'(a)') trim(comment)
    else
      write (ich,*)
    end if
    do j = 1,nat
      if (mask(j)) then
        write (ich,'(1x,a2,1x,3f20.10)') i2e(at(j),'nc'),xyz(1:3,j)
      end if
    end do
    close (ich)
    return
  end subroutine wrxyz_file_mask

!============================================================!
! subroutine wrxyz_channel
! this is the typical quick write routine for xyz files
! version for writing to a output channel
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           comment - (OPTIONAL) the comment line
!
! On Output: file written to "fname"
!============================================================!
  subroutine wrxyz_channel(ch,nat,at,xyz,comment)
    implicit none
    integer :: ch
    integer :: nat
    integer :: at(nat)
    real(wp) ::  xyz(3,nat)
    character(len=*),optional :: comment
    integer :: i,j,k,ich,io
    logical :: ex
    write (ch,'(2x,i0)') nat
    if (present(comment)) then
      write (ch,'(a)') trim(comment)
    else
      write (ch,*)
    end if
    do j = 1,nat
      write (ch,'(1x,a2,1x,3f20.10)') i2e(at(j),'nc'),xyz(1:3,j)
    end do
    return
  end subroutine wrxyz_channel

!============================================================!
! subroutine wrxyz_channel
! this is the typical quick write routine for xyz files
! version for writing to a output channel
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           er   - energy
!
! On Output: file written to "fname"
!============================================================!
  subroutine wrxyz_channel_energy(ch,nat,at,xyz,er)
    implicit none
    integer :: ch
    integer :: nat
    integer :: at(nat)
    real(wp) ::  xyz(3,nat)
    real(wp) :: er
    integer :: i,j,k,ich,io
    logical :: ex
    character(len=30) :: etmp
    write (ch,'(2x,i0)') nat
    write (etmp,'(f20.10)') er
    write (ch,'(2x,a,a)') "energy=",adjustl(etmp)
    do j = 1,nat
      write (ch,'(1x,a2,1x,3f20.10)') i2e(at(j),'nc'),xyz(1:3,j)
    end do
    return
  end subroutine wrxyz_channel_energy

!============================================================!
! subroutine wrsdf_channel
! this is the quick write routine for sdf files
! version for writing to a output channel
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           er   - energy
!           wbo  - bond order matrix
!
! On Output: written to channel "ch"
!============================================================!
  subroutine wrsdf_channel(ch,nat,at,xyz,er,chrg,wbo,comment,icharges)
    implicit none
    integer,intent(in) :: ch
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) ::  xyz(3,nat)
    real(wp),intent(in) :: er
    integer,intent(in) :: chrg
    real(wp),intent(in) :: wbo(nat,nat)
    character(len=*),intent(in) :: comment
    real(wp),intent(in),optional :: icharges(nat)
    character(len=8)  :: date
    character(len=10) :: time
    integer :: list12(12),nbd
    integer,parameter :: list4(4) = 0
    integer,parameter :: list8(8) = 0
    character(len=*),parameter :: countsfmt = '(3i3, 8i3, 1x, a5)'
    character(len=*),parameter :: atmfmt = '(3f10.4, 1x, a2, 12i3)'
    character(len=*),parameter :: bndfmt = '(7i3)'
    integer :: i,j,k,ich,io
    logical :: ex

    !>--- generate data
    call date_and_time(date,time)
    nbd = countbonds(nat,wbo)
    list12 = 0
    !>--- comment lines
    call date_and_time(date,time)
    write (ch,'(a)') trim(comment)
    write (ch,'(1x,a, 3a2, a4, "3D",1x,a,f18.8,5x)') &
    & 'crest',date(5:6),date(7:8),date(3:4),time(:4),'Energy =',er
    write (ch,'(a)')
    !>--- counts line
    write (ch,countsfmt) nat,nbd,list8,999,'V2000'
    !>--- atom block
    do j = 1,nat
      write (ch,atmfmt) xyz(1:3,j),i2e(at(j),'nc'),list12
    end do
    !>--- bonds block
    do i = 1,nat
      do j = i+1,nat
        k = nint(wbo(j,i))
        if (k > 0) then
          write (ch,bndfmt) i,j,k,list4
        end if
      end do
    end do
    !>--- other
    if (present(icharges)) then
      do i = 1,nat
        if (abs(nint(icharges(i))) /= 0) then
          write (ch,'(a, *(i3, 1x, i3, 1x, i3))') "M  CHG",1,i,nint(icharges(i))
        end if
      end do
    else if (chrg .ne. 0) then
      write (ch,'(a, *(i3, 1x, i3, 1x, i3))') "M  CHG",1,1,chrg
    end if
    write (ch,'(a)') 'M  END'
    write (ch,'(a)') '$$$$'
    return
  end subroutine wrsdf_channel

!============================================================!
! subroutine wrsdfV3000_channel
! this is the quick write routine for sdf files
! version for writing to a output channel
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!           at   - atom number as integer
!           xyz  - coordinates (in Angström)
!           er   - energy
!           wbo  - bond order matrix
!
! On Output: written to channel "ch"
!============================================================!
  subroutine wrsdfV3000_channel(ch,nat,at,xyz,er,chrg,wbo,comment)
    implicit none
    integer,intent(in) :: ch
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) ::  xyz(3,nat)
    real(wp),intent(in) :: er
    real(wp),intent(in) :: chrg
    real(wp),intent(in) :: wbo(nat,nat)
    character(len=*),intent(in),optional :: comment
    character(len=8)  :: date
    character(len=10) :: time
    integer :: list12(12),nbd,b
    integer,parameter :: list4(4) = 0
    character(len=*),parameter :: countsfmt = '(3i3, 8i3, 1x, a5)'
    character(len=*),parameter :: countsfmt2 = '(a,2i3, 3i3)'
    character(len=*),parameter :: atmfmt = '(a,1x,i0,1x, a,3f10.4, i2, 11i3)'
    character(len=*),parameter :: bndfmt = '(a,1x,i0,1x,7i3)'
    integer :: i,j,k,ich,io
    logical :: ex

    !>--- generate data
    call date_and_time(date,time)
    nbd = countbonds(nat,wbo)
    !>--- comment lines
    call date_and_time(date,time)
    if (present(comment)) then
      write (ch,'(1x,a)') comment
    else
      write (ch,'(1x,a)') 'structure written by crest'
    end if
    write (ch,'(1x,a,f18.8,5x, 3a2, a4, "3D")') &
    & 'Energy =',er,date(5:6),date(7:8),date(3:4),time(:4)
    write (ch,'(a)')
    !>--- counts line
    write (ch,countsfmt) nat,nbd,0,0,0,999,'V2000'
    write (ch,'("M V30 BEGIN CTAB")')
    write (ch,countsfmt2) "M V30 COUNTS",nat,nbd,0,0,0
    !>--- atom block
    write (ch,'("M V30 BEGIN ATOM")')
    do j = 1,nat
      write (ch,atmfmt) 'M V30',j, &
      &     i2e(at(j),'nc'),xyz(1:3,j),list12
    end do
    write (ch,'("M V30 END ATOM")')
    !>--- bonds block
    write (ch,'("M V30 BEGIN BOND")')
    b = 0
    do i = 1,nat
      do j = i+1,nat
        k = nint(wbo(j,i))
        if (k > 0) then
          b = b+1
          write (ch,bndfmt) "M V30",b,i,j,k,list4
        end if
      end do
    end do
    write (ch,'("M V30 END BOND")')
    !>--- other
    if (chrg .ne. 0) then
      write (ch,'(a, *(i3, 1x, i3, 1x, i3))') "M V30 CHG",1,1,chrg
    end if
    write (ch,'(a)') 'M V30 END CTAB'
    write (ch,'(a)') 'M  END'
    write (ch,'(a)') '$$$$'
    return
  end subroutine wrsdfV3000_channel

! ──────────────────────────────────────────────────────────────────────────────

  subroutine read_extxyz_frame(iunit,ext_sigs,ext_props,nat,energy,lat,success, &
  &                            energy_units,forces_units)
    implicit none

    ! Formal Arguments
    integer,intent(in)          :: iunit
    type(extxyz_signatures),intent(out) :: ext_sigs
    type(extxyz_properties),intent(out) :: ext_props
    integer,intent(out)  :: nat
    real(wp),intent(out) :: energy
    real(wp),intent(out),allocatable :: lat(:,:)
    logical,intent(out)         :: success
    character(len=32),intent(out),optional :: energy_units
    character(len=32),intent(out),optional :: forces_units

    ! Internal variables
    integer                      :: i,ierr,total_fields
    character(len=5000)          :: comment_line
    character(len=2000)          :: val_str
    logical                      :: found
    real(wp)                     :: lattice(3,3)
    real(wp)                     :: lat_raw(9)
    character(len=128),allocatable :: line_fields(:)
    character(len=2000)          :: current_line

    success = .true.

    ! 1. Read Number of Atoms (nat)
    read (iunit,*,iostat=ierr) nat
    if (ierr /= 0) then
      success = .false.
      return
    end if

    ! 2. Read the long comment line
    read (iunit,'(A)',iostat=ierr) comment_line
    if (ierr /= 0) then
      success = .false.
      return
    end if

    ! 3. Extract Key-Value Pairs
    ! Extract Energy
    call get_key_value(comment_line,"energy",val_str,found,case_sensitive=.false.)
    if (found) read (val_str,*) energy

    ! Extract optional unit specifications (defaults: eV for energy, eV/Ang for forces)
    ! Values are stored lowercase so callers can do case-insensitive comparisons.
    if (present(energy_units)) then
      energy_units = 'ev'
      call get_key_value(comment_line,"energy_units",val_str,found,case_sensitive=.false.)
      if (found) energy_units = lowerCase(trim(adjustl(val_str)))
    end if
    if (present(forces_units)) then
      forces_units = 'ev/ang'
      call get_key_value(comment_line,"forces_units",val_str,found,case_sensitive=.false.)
      if (found) forces_units = lowerCase(trim(adjustl(val_str)))
    end if

    ! Extract Lattice
    call get_key_value(comment_line,"lattice",val_str,found,case_sensitive=.false.)
    if (found) then
      read (val_str,*) lat_raw
      lattice = reshape(lat_raw, (/3,3/))
      allocate (lat(3,3))
      lat = lattice*aatoau
    end if

    ! Extract and Parse Properties Signature
    call get_extxyz_value(comment_line,"Properties",val_str,found)
    call get_key_value(comment_line,'Properties',val_str,found,case_sensitive=.false.)
    if (found) then
      call parse_properties_tag(val_str,ext_sigs)
    else
      success = .false.
      return
    end if

    ! 4. Allocate extxyz_properties based on signatures
    call allocate_extxyz_properties_from_sigs(nat,ext_sigs,ext_props)

    ! 5. Read Atom Data Lines
    total_fields = ext_sigs%total_fields
    allocate (line_fields(total_fields))

    do i = 1,nat
      read (iunit,'(A)',iostat=ierr) current_line
      if (ierr /= 0) then
        success = .false.
        exit
      end if

      ! 6. Fill entries in extxyz_properties
      call fill_atom_properties(current_line,ext_props,i)
    end do
  end subroutine read_extxyz_frame

! ──────────────────────────────────────────────────────────────────────────────

  subroutine fill_atom_properties(current_line,ext_props,i)
    implicit none
    character(len=*),intent(in) :: current_line
    type(extxyz_properties),intent(inout) :: ext_props
    integer,intent(in) :: i
    integer :: ii,jj,kk,ierr
    character(len=32),allocatable :: line_fields(:)

    allocate (line_fields(ext_props%total_fields),source=repeat(' ',32))

    read (current_line,*,iostat=ierr) line_fields
    if (ierr /= 0) then
      write (stdout,*) '**ERROR** unexpected line fromat in extxyz parsing for atom',i
    end if

    kk = 0
    do ii = 1,ext_props%n_props
      associate (prop => ext_props%props(ii))
        do jj = 1,prop%signat%n_fields
          kk = kk+1
          select case (prop%signat%p_type)
          case ('S')
            prop%S(jj,i) = trim(line_fields(kk))
          case ('I')
            read (line_fields(kk),*,iostat=ierr) prop%I(jj,i)
          case ('R')
            read (line_fields(kk),*,iostat=ierr) prop%R(jj,i)
          end select
          if (ierr /= 0) then
            write (stdout,*) '**ERROR** unexpected line fromat in extxyz parsing for element',jj,'of atom',i
            return
          end if
        end do
      end associate
    end do

  end subroutine fill_atom_properties

! ──────────────────────────────────────────────────────────────────────────────

  subroutine get_at_from_ext(ext_props,at)
    implicit none
    type(extxyz_properties) :: ext_props
    integer,intent(out),allocatable :: at(:)
    integer :: ii,jj,nat
    do ii = 1,ext_props%n_props
      associate (prop => ext_props%props(ii))
        select case (trim(prop%signat%name))
        case ('species')
          nat = prop%natoms
          allocate (at(nat),source=0)
          do jj = 1,nat
            at(jj) = e2i(trim(prop%S(1,jj)))
          end do
        end select
      end associate
    end do
  end subroutine get_at_from_ext

  subroutine get_xyz_from_ext(ext_props,xyz)
    implicit none
    type(extxyz_properties) :: ext_props
    real(wp),intent(out),allocatable :: xyz(:,:)
    integer :: ii,jj,nat
    do ii = 1,ext_props%n_props
      associate (prop => ext_props%props(ii))
        select case (trim(prop%signat%name))
        case ('pos')
          nat = prop%natoms
          allocate (xyz(3,nat),source=0.0_wp)
          do jj = 1,nat
            xyz(:,jj) = prop%R(:,jj)*aatoau
          end do
        end select
      end associate
    end do
  end subroutine get_xyz_from_ext

  subroutine get_grad_from_ext(ext_props,grad,forces_units)
    implicit none
    type(extxyz_properties) :: ext_props
    real(wp),intent(out),allocatable :: grad(:,:)
    character(len=32),intent(in),optional :: forces_units
    integer :: ii,jj,nat
    character(len=32) :: units_loc
    units_loc = 'eV/Ang'
    if (present(forces_units)) units_loc = forces_units
    do ii = 1,ext_props%n_props
      associate (prop => ext_props%props(ii))
        select case (trim(prop%signat%name))
        case ('forces')
          nat = prop%natoms
          allocate (grad(3,nat),source=0.0_wp)
          do jj = 1,nat
            select case (trim(units_loc))
            case ('ha/bohr','hartree/bohr','au')
              grad(:,jj) = prop%R(:,jj)*(-1.0_wp)  ! forces → gradient (sign flip only)
            case default                             ! 'ev/ang' and anything unrecognised
              grad(:,jj) = prop%R(:,jj)*(-autoaa/autoeV)
            end select
          end do
        end select
      end associate
    end do
  end subroutine get_grad_from_ext

!=========================================================================================!
!=========================================================================================!
!  4. GENERAL UTILITY ROUTINES
!=========================================================================================!
!=========================================================================================!

!============================================================!
! read a line of coordinates and determine by itself
! if the format is x,y,z,at or at,x,y,z
!============================================================!
  subroutine coordline(line,sym,xyz,io)
    implicit none
    character(len=*) :: line
    character(len=*) :: sym
    real(wp) :: xyz(3)
    integer,intent(out) :: io

    io = 0
    read (line,*,iostat=io) xyz(1:3),sym
    if (io .ne. 0) then
      read (line,*,iostat=io) sym,xyz(1:3)
    end if

    return
  end subroutine coordline

!============================================================!
! convert a string into uppercase
!============================================================!
  function upperCase(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: upperCase
    integer :: ic,i
    character(26),Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,LEN_TRIM(s)
      ic = INDEX(low,s(i:i))
      if (ic > 0) sout(i:i) = high(ic:ic)
    end do
    call move_alloc(sout,upperCase)
  end function upperCase

!============================================================!
! convert a string into lowercase
!============================================================!
  function lowerCase(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: lowerCase
    integer :: ic,i
    character(26),Parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,LEN_TRIM(s)
      ic = INDEX(high,s(i:i))
      if (ic > 0) sout(i:i) = low(ic:ic)
    end do
    call move_alloc(sout,lowerCase)
  end function lowerCase

!============================================================!
! split element lable if some isotope indicator was given
! and convert to uppercase
!============================================================!
  function convertlable(s)
    implicit none
    character(len=*),intent(in) :: s
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: convertlable
    integer :: ic,i
    character(14),parameter :: lab = '0123456789*_+-'
    character(26),parameter :: high = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26),parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
    sout = s
    do i = 1,len_trim(s)
      ic = index(lab,s(i:i))
      if (ic > 0) sout(i:i) = ' '
      ic = index(low,s(i:i))
      if (ic > 0) sout(i:i) = high(ic:ic)
    end do
    sout = trim(adjustl(sout))
    if (len_trim(sout) .gt. 1) then
      sout(2:2) = lowerCase(sout(2:2))
    else
      sout = sout//' '
    end if
    call move_alloc(sout,convertlable)
  end function convertlable

!============================================================!
! e2i is used to map the element (as a string) to integer
!============================================================!
  integer function e2i(cin)
    implicit none
    character(len=*),intent(in) :: cin
    character(len=:),allocatable :: c
    integer :: iout
    integer :: i,j,k,ich,io,Z
    logical :: ex
    c = trim(convertlable(cin))
    read (cin,*,iostat=io) j
    if (io == 0) Z = j
    if (any(PSE(:) .eq. c)) then
      do i = 1,118
        if (trim(PSE(i)) .eq. c) then
          iout = i
          exit
        end if
      end do
    else if (io == 0.and.Z <= 118) then
      iout = Z
    else !> special cases
      select case (trim(c))
      case ('D'); iout = 1
      case ('T'); iout = 1
      case default; iout = 0
      end select
    end if
    e2i = iout
  end function e2i

!============================================================!
! i2e is used to map the element (as a integer) to a string
!============================================================!
  character(len=2) function i2e(iin,oformat)
    implicit none
    integer,intent(in) :: iin
    character(len=:),allocatable :: c
    character(len=*),optional :: oformat
    if (iin <= 118) then
      c = uppercase(PSE(iin))
    else
      c = 'XX'
    end if
    i2e = trim(c)
    if (present(oformat)) then
      select case (oformat)
      case ('lc','lowercase')
        i2e = lowerCase(trim(c))
      case ('nc','nicecase')
        if (len_trim(c) .gt. 1) then
          c(2:2) = lowerCase(c(2:2))
          i2e = trim(c)
        end if
      case default
        continue
      end select
    end if
  end function i2e

!============================================================!
! get the file extension
!============================================================!
  function fextension(s)
    implicit none
    character(len=*),intent(in) :: s !filename
    character(len=:),allocatable :: sout
    character(len=:),allocatable :: fextension !output
    integer :: ic,i
    sout = trim(adjustl(s))
    i = len_trim(sout)
    ic = index(sout,'.',.true.)
    if (ic .ne. 0) then
      fextension = sout(ic:i)
    else
      fextension = 'none'
    end if
    return
  end function fextension

!============================================================!
! grep for a keyword within the file
!============================================================!
  function sgrep(fname,key,casesensitive)
    implicit none
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) :: key
    logical,intent(in),optional :: casesensitive
    logical :: sgrep,ex
    logical :: convert
    character(len=5000) :: atmp
    character(len=:),allocatable :: kkey
    integer :: ic,io
    sgrep = .false.
    convert = .false.
    inquire (file=fname,exist=ex)
    if (.not.ex) return
    kkey = trim(key)
    if (present(casesensitive)) then
      if (.not.casesensitive)then
        kkey = lowercase(key)
        convert = .true.
      endif
    end if
    open (newunit=ic,file=fname)
    do
      read (ic,'(a)',iostat=io) atmp
      if (io < 0) exit !EOF
      if (convert) then
        if (index(lowercase(atmp),kkey) .ne. 0) then
          sgrep = .true.
          exit
        end if
      else
        if (index(atmp,kkey) .ne. 0) then
          sgrep = .true.
          exit
        end if
      end if
    end do
    close (ic)
    return
  end function sgrep

!============================================================!
! grep the energy from a line of strings
!============================================================!
  function grepenergy(line)
    implicit none
    real(wp) :: grepenergy
    character(len=*),intent(in) :: line
    real(wp) :: energy
    character(len=:),allocatable :: atmp
    integer :: i,io,k
    atmp = trim(line)
    energy = 0.0_wp
    if (index(atmp,'energy=') .ne. 0) then
      k = index(atmp,'energy=')
      atmp = atmp(k+7:)
      read (atmp,*,iostat=io) energy
      if (io .ne. 0) energy = 0.0_wp
    else if (index(atmp,'energy:') .ne. 0) then
      k = index(atmp,'energy:')
      atmp = atmp(k+7:)
      read (atmp,*,iostat=io) energy
      if (io .ne. 0) energy = 0.0_wp
    else
      !> assumes that the first float in the line is the energy
      do i = 1,len_trim(atmp)
        if (len_trim(atmp) .lt. 1) exit
        read (atmp,*,iostat=io) energy
        if (io > 0) then
          atmp = atmp(2:)
          atmp = adjustl(atmp)
          cycle
        else
          exit
        end if
      end do
    end if
    grepenergy = energy
    return
  end function grepenergy
! ──────────────────────────────────────────────────────────────────────────────

  subroutine get_extxyz_value(comment_line,key,value,found)
!*************************************************************************
!* subroutine get_extxyz_value                                           *
!* grep a key-value-pair from the comment line of an extended XYZ file   *
!* On input:                                                             *
!*      comment_line - the comment line                                  *
!*      key          - the key to look for (case INSENSITIVE)            *
!*                                                                       *
!* On output:                                                            *
!*      value - the value as raw string                                  *
!*      found - success logical, did we find the key?                    *
!*************************************************************************
    implicit none
    character(len=*),intent(in)  :: comment_line
    character(len=*),intent(in)  :: key
    character(len=*),intent(out) :: value
    logical,intent(out)          :: found

    integer :: key_start,val_start,val_end,line_len
    character(len=:),allocatable :: search_key

    found = .false.
    value = ""
    line_len = len_trim(comment_line)

    search_key = lowercase(key)//"="
    key_start = index(comment_line,trim(search_key))

    if (key_start > 0) then
      val_start = key_start+len_trim(search_key)

      ! --- Skip any spaces between '=' and the value
      do while (val_start <= line_len.and.comment_line(val_start:val_start) == " ")
        val_start = val_start+1
      end do

      ! If we hit the end of the line, the key had no value
      if (val_start > line_len) return
      found = .true.

      ! Check for quotes
      if (comment_line(val_start:val_start) == '"'.or. &
          comment_line(val_start:val_start) == "'") then

        val_start = val_start+1
        val_end = val_start+index(comment_line(val_start:),comment_line(val_start-1:val_start-1))-2
      else
        ! Bare value: find next space
        val_end = val_start+index(comment_line(val_start:)," ")-2
        if (val_end < val_start) val_end = line_len
      end if

      value = comment_line(val_start:val_end)
    end if
  end subroutine get_extxyz_value

  subroutine get_key_value(input_str,key,value,success,case_sensitive)
    implicit none

    ! Arguments
    character(len=*),intent(in)          :: input_str     ! The full string to search
    character(len=*),intent(in)          :: key           ! The key to look for
    character(len=*),intent(out)         :: value         ! The extracted value
    logical,intent(out)         :: success       ! True if found and parsed
    logical,intent(in),optional :: case_sensitive ! Toggle case sensitivity

    ! Internal variables
    integer :: key_len,str_len,start_pos,val_start,val_end
    character(len=len(input_str)) :: search_str,search_key
    character(len=1) :: quote_char
    logical :: sensitive

    ! Initialize
    success = .false.
    value = ""
    sensitive = .true.
    if (present(case_sensitive)) sensitive = case_sensitive

    key_len = len_trim(key)
    str_len = len_trim(input_str)

    ! Prepare strings for case-insensitive search if requested
    if (.not.sensitive) then
      search_str = lowercase(input_str)
      search_key = lowercase(key)
    else
      search_str = input_str
      search_key = key
    end if

    ! Find key followed immediately by '='
    ! Note: Searching for 'key=' to satisfy the "no whitespace" requirement
    start_pos = index(search_str,trim(search_key)//'=')

    if (start_pos > 0) then
      ! Value starts right after 'key='
      val_start = start_pos+key_len+1

      ! Check for quotes
      quote_char = input_str(val_start:val_start)

      if (quote_char == '"'.or.quote_char == "'") then
        ! Handle quoted value
        val_start = val_start+1
        ! Find the closing quote starting from the next character
        val_end = index(input_str(val_start:),quote_char)

        if (val_end > 0) then
          val_end = val_start+val_end-2
          value = input_str(val_start:val_end)
          success = .true.
        end if
      else
        ! Handle unquoted value (extract until next space or end of string)
        val_end = index(input_str(val_start:)," ")
        if (val_end == 0) then
          val_end = str_len
        else
          val_end = val_start+val_end-2
        end if
        value = input_str(val_start:val_end)
        success = .true.
      end if
    end if
  end subroutine get_key_value

! ──────────────────────────────────────────────────────────────────────────────

  function count_extxyz_pairs(comment_line) result(num_pairs)
    implicit none
    character(len=*),intent(in) :: comment_line
    integer :: num_pairs
    integer :: i,line_len
    logical :: in_quotes
    character :: quote_char

    num_pairs = 0
    in_quotes = .false.
    line_len = len_trim(comment_line)
    quote_char = ' '

    do i = 1,line_len
      ! Check if we are entering or leaving a quoted section
      if (.not.in_quotes) then
        if (comment_line(i:i) == '"'.or.comment_line(i:i) == "'") then
          in_quotes = .true.
          quote_char = comment_line(i:i)
        end if
      else
        ! If we are in quotes, look for the matching closing quote
        if (comment_line(i:i) == quote_char) then
          in_quotes = .false.
        end if
      end if

      ! If we find an '=' while NOT in quotes, it's a new key-value pair
      if (.not.in_quotes.and.comment_line(i:i) == '=') then
        num_pairs = num_pairs+1
      end if
    end do
  end function count_extxyz_pairs

! ──────────────────────────────────────────────────────────────────────────────

!============================================================!
! count number of bonds from an wbo matrix
!============================================================!
  function countbonds(nat,wbo) result(nbd)
    implicit none
    integer,intent(in)  :: nat
    real(wp),intent(in) :: wbo(nat,nat)
    integer :: nbd
    integer :: i,j,k
    nbd = 0
    do i = 1,nat
      do j = 1,i-1
        k = nint(wbo(i,j))
        if (k > 0) nbd = nbd+1
      end do
    end do
    return
  end function countbonds

!=========================================================================================!

  subroutine get_atlist(nat,atlist,line,at)
!******************************************************
!* Analyze a string containing atom specifications.
!* "atlist" is a array of booleans for each atom,
!* which is set to .true. should the atom be contained
!* in atlist.
!******************************************************
    implicit none
    integer,intent(in) :: nat
    logical,intent(out),allocatable :: atlist(:)
    character(len=*),intent(in) :: line
    integer,intent(in),optional :: at(nat)
    character(len=:),allocatable :: substr(:)
    integer :: i,j,k,l,io,ns,ll,i1,i2,io1,io2,i3,i4
    character(len=:),allocatable :: atmp,btmp

    allocate (atlist(nat),source=.false.)
!>-- count stuff
    ll = len_trim(line)
    ns = 1
    do i = 1,ll
      if (line(i:i) .eq. ',') ns = ns+1
    end do
    allocate (substr(ns),source=repeat(' ',ll))
!>-- cut stuff
    if (ns > 1) then
      j = 1
      k = 1
      do i = 1,ll
        if (k == ns) then
          substr(k) = lowercase(adjustl(line(j:)))
          exit
        end if
        if (line(i:i) .eq. ',') then
          substr(k) = lowercase(adjustl(line(j:i-1)))
          k = k+1
          j = i+1
        end if
      end do
    else
      substr(1) = trim(line)
    end if
!>--- analyze stuff
    do i = 1,ns
      atmp = trim(substr(i))
      if (atmp .eq. 'all') then
        atlist(:) = .true.
        exit
      end if
      if (index(atmp,'.') .ne. 0) cycle !> exclude floats
      l = index(atmp,'-')
      if (l .eq. 0) then
        read (atmp,*,iostat=io) i1
        !> check if it is an element symbol
        if (io /= 0) then
          if (len_trim(atmp) > 2) then
            if (index(trim(atmp),'heavy') .ne. 0) then !> all heavy atoms
              if (present(at)) then
                do j = 1,nat
                  if (at(j) > 1) atlist(j) = .true.
                end do
              end if
            end if
          else !> element symbols
            k = e2i(atmp)
            if (present(at)) then
              do j = 1,nat
                if (at(j) == k) atlist(j) = .true.
              end do
            end if
          end if
        else
          atlist(i1) = .true.
        end if
      else
        btmp = atmp(:l-1)
        read (btmp,*,iostat=io1) i1
        btmp = atmp(l+1:)
        read (btmp,*,iostat=io2) i2
        if (io1 .eq. 0.and.io2 .eq. 0) then
          i4 = max(i1,i2)
          i3 = min(i1,i2)
          do j = 1,nat
            if (i3 <= j.and.j <= i4) atlist(j) = .true.
          end do
        end if
      end if
    end do
    deallocate (substr)
  end subroutine get_atlist

!=========================================================================================!
  function sumform(nat,at) result(sumformula)
!************************************************
!* get sumformula as a string from the AT array
!************************************************
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    character(len=:),allocatable :: sumformula
    integer :: sumat(118)
    integer :: i
    character(len=6) :: str
    sumformula = ''
    sumat = 0
    do i = 1,nat
      sumat(at(i)) = sumat(at(i))+1
    end do
    !> carbon always first
    if (sumat(6) > 0) then
      if (sumat(6) > 1) then
        write (str,'(a,i0)') trim(adjustl(i2e(6,'nc'))),sumat(6)
      else
        str = 'C'
      end if
      sumformula = trim(sumformula)//trim(str)
    end if
    do i = 2,118
      if (i == 6) cycle
      if (sumat(i) .lt. 1) cycle
      if (sumat(i) > 1) then
        write (str,'(a,i0)') trim(adjustl(i2e(i,'nc'))),sumat(i)
      else
        str = trim(i2e(i,'nc'))
      end if
      sumformula = trim(sumformula)//trim(str)
    end do
    !> hydrogen always last
    if (sumat(1) > 0) then
      if (sumat(1) > 1) then
        write (str,'(a,i0)') trim(adjustl(i2e(1,'nc'))),sumat(1)
      else
        str = 'H'
      end if
      sumformula = trim(sumformula)//trim(str)
    end if
    return
  end function sumform

! ══════════════════════════════════════════════════════════════════════════════
! end of the module
! ══════════════════════════════════════════════════════════════════════════════
end module molecule_io
