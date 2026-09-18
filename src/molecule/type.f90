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

module molecule_type
  use iso_c_binding
  use molecule_parameters
  use molecule_io
  use molecule_type_components
!> simple geomerty and vector operations
  use geo
!> element symbols
  use crest_cn_module,only:calculate_cn
  implicit none
  private
! ══════════════════════════════════════════════════════════════════════════════
  !> EXPORTS
  public :: coord
  public :: coord2xyz,xyz2coord
! ══════════════════════════════════════════════════════════════════════════════

  type :: coord
    !> coord class. contains a single structure
    !> by convention coordinates are in atomic units (Bohr) for a single structure!

    !********************************************!
    !> data that's typically used in coord type <!
    !********************************************!
    !>-- number of atoms
    integer :: nat = 0
    !>-- atom types as integer, dimension will be at(nat)
    integer,allocatable  :: at(:)
    !>-- atomic coordinates, by convention in Bohrs
    real(wp),allocatable :: xyz(:,:)

    !**************************************!
    !> (optional) data, often not present <!
    !**************************************!
    !>-- energy
    real(wp) :: energy = 0.0_wp
    !>-- gradient
    real(wp),allocatable :: gradient(:,:)
    !>-- a comment line
    character(len=:),allocatable :: comment
    !>-- "origin" tag
    character(len=:),allocatable :: origin
    !>-- molecular charge
    integer :: chrg = 0
    !>-- multiplicity information
    integer :: uhf = 0
    !>-- number of bonds
    integer :: nbd = 0
    !>-- bond info
    integer,allocatable :: bond(:,:)
    !>-- lattice vectors
    real(wp),allocatable :: lat(:,:)

    !>-- atomic charges
    real(wp),allocatable :: qat(:)

    !>-- atom ids (for iRMSDs)
    integer,allocatable :: id(:)

    !>-- (optional) PDB data
    type(pdbdata) :: pdb

    !>-- extxyz signature
    logical :: wrextxyz = .false.
    type(extxyz_signatures),allocatable :: extxyz

  contains
    procedure :: deallocate => deallocate_coord !> clear memory space
    procedure :: open => opencoord              !> read an coord file
    procedure :: write => writecoord            !> write (detected from file extension)
    procedure :: writeextxyz => write_extxyz    !> write extxyz file to a given iunit
    procedure :: append => appendcoord          !> append
    procedure :: get => getcoord                !> allocate & fill with data
    procedure :: appendlog => appendcoord       !> append .log file with coordinates and energy
    procedure :: dist => coord_getdistance      !> calculate distance between two atoms
    procedure :: angle => coord_getangle        !> calculate angle between three atoms
    procedure :: dihedral => coord_getdihedral  !> calculate dihedral angle between four atoms
    procedure :: cutout => coord_getcutout      !> create a substructure
    procedure :: get_CN => coord_get_CN         !> calculate coordination number
    procedure :: cellvol => coord_cellvol       !> cell volume |det(lat)| in Bohr³
    procedure :: get_z => coord_get_z           !> calculate nuclear charge
    procedure :: cn_to_bond => coord_cn_to_bond !> generate neighbour matrix from CN
    procedure :: swap => atswp                  !> swap two atoms coordinates and their at() entries
    procedure :: sumform => coord_sumform       !> generate a string with the sum formula
    procedure :: copy => coord_copy             !> deep copy from another coord object
  end type coord

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════
!  ROUTINES FOR READING SINGLE STRUCTURES (COORDS)
! ──────────────────────────────────────────────────────────────────────────────

  subroutine deallocate_coord(self)
!**********************************************
!* subroutine deallocate_coord                *
!* is used to clear memory for the coord type *
!**********************************************
    implicit none
    class(coord) :: self
    self%nat = 0
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    call self%pdb%deallocate()
    return
  end subroutine deallocate_coord

! ──────────────────────────────────────────────────────────────────────────────

  subroutine coord_copy(self,src)
!*************************************************************
!* Deep copy of a coord object from src to self.            *
!*                                                          *
!* On Input:  src  - source coord object                    *
!* On Output: self - destination, populated with src data   *
!*************************************************************
    implicit none
    class(coord),intent(out) :: self
    type(coord),intent(in)   :: src

    ! ── scalar fields ────────────────────────────────────────────────────────
    self%nat      = src%nat
    self%energy   = src%energy
    self%chrg     = src%chrg
    self%uhf      = src%uhf
    self%nbd      = src%nbd
    self%wrextxyz = src%wrextxyz

    ! ── mandatory allocatable arrays ─────────────────────────────────────────
    if (allocated(src%at))  self%at  = src%at
    if (allocated(src%xyz)) self%xyz = src%xyz

    ! ── optional allocatable arrays ──────────────────────────────────────────
    if (allocated(src%gradient)) self%gradient = src%gradient
    if (allocated(src%bond))     self%bond     = src%bond
    if (allocated(src%lat))      self%lat      = src%lat
    if (allocated(src%qat))      self%qat      = src%qat

    ! ── optional character fields ─────────────────────────────────────────────
    if (allocated(src%comment)) self%comment = src%comment
    if (allocated(src%origin))  self%origin  = src%origin

    ! ── derived-type components ───────────────────────────────────────────────
    self%pdb = src%pdb
    if (allocated(src%extxyz)) self%extxyz = src%extxyz

    return
  end subroutine coord_copy

! ──────────────────────────────────────────────────────────────────────────────

  subroutine opencoord(self,fname)
!************************************************
!* subroutine opencoord                         *
!* is the open procedure for the "coord" class. *
!************************************************
    implicit none
    class(coord) :: self
    character(len=*),intent(in) :: fname
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp),allocatable :: grad(:,:)
    real(wp),allocatable :: lat(:,:)
    integer :: ftype
    integer :: i,j,k,ich,io,iunit
    logical :: ex,success
    real(wp) :: en
    character(len=32) :: eu,fu
    type(extxyz_signatures) :: ext_sigs
    type(extxyz_properties) :: ext_props

    inquire (file=fname,exist=ex)
    if (.not.ex) then
      write (stdout,'(a)') '**ERROR** could not find coord file '//trim(fname)
      call exit(1)
    end if

    call self%deallocate()

    call checkcoordtype(fname,ftype)
    call rdnat(fname,nat,ftype=ftype)

    if (nat > 0) then
      en = 0.0_wp
      allocate (at(nat),xyz(3,nat))
      select case (ftype)
      case (coordtype%PDB)
        call rdPDB(fname,nat,at,xyz,self%pdb) ! ← need to fill self%pdb
        xyz = xyz/bohr

      case (coordtype%extxyz)
        open (newunit=iunit,file=fname)
        call read_extxyz_frame(iunit,ext_sigs,ext_props,nat,en,lat,success, &
        &                      energy_units=eu,forces_units=fu)
        close (iunit)
        if (success) then
          select case (trim(eu))
          case ('hartree','ha','au')
            ! energy already in Hartree, no conversion needed
          case default  !> 'ev' and anything unrecognised
            en = en/autoeV
          end select
          call get_at_from_ext(ext_props,at)
          call get_xyz_from_ext(ext_props,xyz)
          call get_grad_from_ext(ext_props,grad,forces_units=fu)
          if (allocated(lat)) call move_alloc(lat,self%lat)
          if (allocated(grad)) call move_alloc(grad,self%gradient)
        end if

      case default
        call rdcoord(fname,nat,at,xyz,energy=en,ftype=ftype)

      end select
      self%nat = nat
      self%energy = en
      call move_alloc(at,self%at)
      call move_alloc(xyz,self%xyz)
    else
      write (stdout,'(a)') '**ERROR** Format issue while reading coord file '//trim(fname)
      write (stdout,'(a)') '          Number of atoms detected as zero!'
      call exit(1)
    end if

    return
  end subroutine opencoord

! ──────────────────────────────────────────────────────────────────────────────

! subroutine getcoord
! allocate "coord" class and fill with data
  subroutine getcoord(self,convfac,nat,at,xyz)
    implicit none
    class(coord) :: self
    real(wp),intent(in) :: convfac
    integer,intent(in)  :: nat
    integer,intent(in)  :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    call self%deallocate()
    allocate (self%at(nat))
    allocate (self%xyz(3,nat))
    self%nat = nat
    self%at = at
    self%xyz = xyz/convfac
    return
  end subroutine getcoord

! ──────────────────────────────────────────────────────────────────────────────

! function coord_getdistance
! calculate the distance for a given pair of atoms
  function coord_getdistance(self,a1,a2) result(d)
    implicit none
    class(coord) :: self
    integer,intent(in) :: a1,a2
    real(wp) :: d
    d = 0.0_wp
    if (allocated(self%xyz)) then
      d = (self%xyz(1,a1)-self%xyz(1,a2))**2+ &
      &   (self%xyz(2,a1)-self%xyz(2,a2))**2+ &
      &   (self%xyz(3,a1)-self%xyz(3,a2))**2
      d = sqrt(d)
    end if
    return
  end function coord_getdistance

! ──────────────────────────────────────────────────────────────────────────────

! function coord_getangle
! calculate the angle for a given trio of atoms in rad
! A1-A2-A3
  function coord_getangle(self,a1,a2,a3) result(angle)
    implicit none
    class(coord) :: self
    integer,intent(in) :: a1,a2,a3
    real(wp) :: angle,u(3),v(3),o(3)
    real(wp) :: d2ij,d2jk,d2ik,xy,temp
    angle = 0.0_wp
    if (allocated(self%xyz)) then
      u(1:3) = self%xyz(1:3,a1)-self%xyz(1:3,a2)
      v(1:3) = self%xyz(1:3,a3)-self%xyz(1:3,a2)
      angle = tangle(u,v)
    end if
    return
  end function coord_getangle

! ──────────────────────────────────────────────────────────────────────────────

! function coord_getdihedral
! calculate the dihedral angle for a given quartet of atoms in rad
! A1-A2-A3-A4
  function coord_getdihedral(self,a1,a2,a3,a4) result(dihed)
    implicit none
    class(coord) :: self
    integer,intent(in) :: a1,a2,a3,a4
    real(wp) :: dihed
    real(wp) :: u(3),v(3),w(3)
    real(wp) :: n1(3),n2(3)
    real(wp) :: u1(3),u2(3),u3(3)

    dihed = 0.0_wp
    if (allocated(self%xyz)) then

      u(1:3) = self%xyz(1:3,a2)-self%xyz(1:3,a1)
      v(1:3) = self%xyz(1:3,a3)-self%xyz(1:3,a2)
      w(1:3) = self%xyz(1:3,a4)-self%xyz(1:3,a3)
      dihed = dihedral(u,v,w)
    end if
    return
  end function coord_getdihedral

! ──────────────────────────────────────────────────────────────────────────────
! function coord_getgutout
! create a cutout mol object
  function coord_getcutout(self,atlist) result(molout)
    implicit none
    class(coord) :: self
    logical,intent(in) :: atlist(self%nat)
    type(coord) :: molout
    integer :: newnat,i,j,k,l

    newnat = count(atlist,1)
    if (newnat == self%nat) then
      molout = self
    else
      molout%nat = newnat
      allocate (molout%at(newnat),source=0)
      allocate (molout%xyz(3,newnat),source=0.0_wp)
      k = 0
      do i = 1,self%nat
        if (atlist(i)) then
          k = k+1
          molout%at(k) = self%at(i)
          molout%xyz(1:3,k) = self%xyz(1:3,i)
        end if
      end do
    end if
    return
  end function coord_getcutout

! ──────────────────────────────────────────────────────────────────────────────
  subroutine coord_get_CN(self,cn,cn_type,cn_thr,dcndr)
    implicit none
    class(coord) :: self
    real(wp),intent(out),allocatable :: cn(:)
    real(wp),intent(in),optional :: cn_thr
    character(len=*),intent(in),optional :: cn_type
    real(wp),intent(out),optional :: dcndr(3,self%nat,self%nat)
    if (self%nat <= 0) return
    if (.not.allocated(self%xyz).or..not.allocated(self%at)) return
    allocate (cn(self%nat),source=0.0_wp)
    if (allocated(self%lat)) then
      call calculate_CN(self%nat,self%at,self%xyz,cn, &
      & cntype=cn_type,cnthr=cn_thr,dcndr=dcndr,lat=self%lat)
    else
      call calculate_CN(self%nat,self%at,self%xyz,cn, &
      & cntype=cn_type,cnthr=cn_thr,dcndr=dcndr)
    end if
  end subroutine coord_get_CN

! ──────────────────────────────────────────────────────────────────────────────
  function coord_cellvol(self) result(vol)
    !********************************************************************
    !* Cell volume of a periodic structure as the absolute value of the
    !* lattice determinant (triple product of the three lattice vectors).
    !* The lattice columns are the lattice vectors, matching the
    !* convention used elsewhere (e.g. matmul(lat, sfrac) in irmsd_module).
    !* Returns 0.0 in Bohr³ if no lattice is allocated (non-periodic).
    !********************************************************************
    implicit none
    class(coord) :: self
    real(wp) :: vol
    real(wp) :: a(3,3)
    vol = 0.0_wp
    if (.not.allocated(self%lat)) return
    a = self%lat
    vol = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
      & -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
      & +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    vol = abs(vol)
  end function coord_cellvol

! ──────────────────────────────────────────────────────────────────────────────
  subroutine coord_get_z(self,z)
    implicit none
    class(coord) :: self
    real(wp),intent(out),allocatable :: z(:)
    integer :: i,j,k
    if (self%nat <= 0) return
    if (.not.allocated(self%xyz).or..not.allocated(self%at)) return
    allocate (z(self%nat),source=0.0_wp)
    do i = 1,self%nat
      z(i) = real(self%at(i),wp)-real(ncore(self%at(i)))
      if (self%at(i) > 57.and.self%at(i) < 72) z(i) = 3.0_wp
    end do
  end subroutine coord_get_z

! ──────────────────────────────────────────────────────────────────────────────
  subroutine coord_cn_to_bond(self,cn,bond,cn_type,cn_thr)
    implicit none
    class(coord) :: self
    real(wp),intent(out),allocatable :: cn(:)
    real(wp),intent(out),allocatable,optional :: bond(:,:)
    real(wp),intent(in),optional :: cn_thr
    character(len=*),intent(in),optional :: cn_type
    if (self%nat <= 0) return
    if (.not.allocated(self%xyz).or..not.allocated(self%at)) return
    allocate (cn(self%nat),source=0.0_wp)
    if (allocated(self%lat)) then
      call calculate_CN(self%nat,self%at,self%xyz,cn, &
      & cntype=cn_type,cnthr=cn_thr,bond=bond,lat=self%lat)
    else
      call calculate_CN(self%nat,self%at,self%xyz,cn, &
      & cntype=cn_type,cnthr=cn_thr,bond=bond)
    end if
  end subroutine coord_cn_to_bond

! ══════════════════════════════════════════════════════════════════════════════
!  ROUTINES FOR WRITING STRUCTURES AND CONVERTING THEM
! ══════════════════════════════════════════════════════════════════════════════

  subroutine write_extxyz(self,iunit)
!*************************************************************************
!* Write an extended xyz file from the coord object.                     *
!* Output units are controlled by the module variable extxyz_units_global*
!* (default 'Hartree'). Set to 'eV' to use eV/Ang conventions instead.  *
!*************************************************************************
    class(coord) :: self
    integer,intent(in) :: iunit !> assume the unit is open for writing

    character(len=200) :: atmp
    integer :: ii
    logical :: use_hartree

    use_hartree = (trim(extxyz_units_global) .ne. 'ev')

    !> print number of atoms
    write (iunit,'(i10)') self%nat

    !> construct ext comment line bit by bit
    if (use_hartree) then
      write (atmp,'(f20.10)') self%energy
      write (iunit,'(a,a)',advance='no') trim('energy='//adjustl(atmp)),' '
      write (iunit,'(a)',advance='no') 'energy_units=Hartree '
    else
      write (atmp,'(f20.10)') self%energy*autoeV
      write (iunit,'(a,a)',advance='no') trim('energy='//adjustl(atmp)),' '
      write (iunit,'(a)',advance='no') 'energy_units=eV '
    end if
    if (allocated(self%lat)) then
      write (iunit,'(a)',advance='no') 'Lattice="'
      write (iunit,'(9f15.8)',advance='no') reshape(self%lat*autoaa, [9])
      write (iunit,'(a)',advance='no') '"  pbc="T T T"  '
    end if
    if (allocated(self%gradient)) then
      if (use_hartree) then
        write (iunit,'(a)',advance='no') 'forces_units=Ha/Bohr '
      else
        write (iunit,'(a)',advance='no') 'forces_units=eV/Ang '
      end if
    end if
    if (allocated(self%extxyz)) then
      call assemble_properties_tag(self%extxyz,atmp)
    else if (allocated(self%gradient)) then
      write (atmp,'("species:S:1:pos:R:3:forces:R:3")')
    else
      write (atmp,'("species:S:1:pos:R:3")')
    end if
    write (iunit,'(a,a,a)',advance='no') 'Properties=',trim(atmp),' '
    write (iunit,*)

    !> coord block
    if (allocated(self%extxyz)) then
      write (stdout,*) '**ERROR** This extxyz write function is TODO'
      call exit(1)
    else if (allocated(self%gradient)) then
      if (use_hartree) then
        do ii = 1,self%nat
          !> positions in Ang, forces in Ha/Bohr (sign flip: forces = -gradient)
          write (iunit,'(1x,a2,1x,6f20.10)')  &
          &  i2e(self%at(ii)),self%xyz(1:3,ii)*autoaa,self%gradient(1:3,ii)*(-1.0_wp)
        end do
      else
        do ii = 1,self%nat
          !> positions in Ang, forces in eV/Ang
          write (iunit,'(1x,a2,1x,6f20.10)')  &
          &  i2e(self%at(ii)),self%xyz(1:3,ii)*autoaa,self%gradient(1:3,ii)*(-autoeV/autoaa)
        end do
      end if
    else
      do ii = 1,self%nat
        write (iunit,'(1x,a2,1x,3f20.10)') i2e(self%at(ii)),self%xyz(1:3,ii)*autoaa
      end do
    end if
  end subroutine write_extxyz

! ──────────────────────────────────────────────────────────────────────────────

  subroutine xyz2coord(iname,oname)
!***********************************************
!* subroutine xyz2coord                        *
!* simple conversion of a xyz to a coord file. *
!*                                             *
!* On Input: iname  - name of the xyz file     *
!*           oname  - name of the coord file   *
!*                                             *
!* On Output: file written to "oname"          *
!***********************************************
    implicit none
    character(len=*) :: iname
    character(len=*) :: oname
    type(coord) :: struc
    call struc%open(iname)
    call wrc0(oname,struc%nat,struc%at,struc%xyz)
    call struc%deallocate()
    return
  end subroutine xyz2coord

! ──────────────────────────────────────────────────────────────────────────────
  subroutine coord2xyz(iname,oname)
!***********************************************
!* subroutine coord2xyz                        *
!* simple conversion of a coord to a xyz file. *
!*                                             *
!* On Input: iname  - name of the coord file   *
!*           oname  - name of the xyz file     *
!*                                             *
!* On Output: file written to "oname"          *
!***********************************************
    implicit none
    character(len=*) :: iname
    character(len=*) :: oname
    type(coord) :: struc
    call struc%open(trim(iname))
    struc%xyz = struc%xyz*bohr !to Angström
    call wrxyz(oname,struc%nat,struc%at,struc%xyz)
    call struc%deallocate()
    return
  end subroutine coord2xyz

! ──────────────────────────────────────────────────────────────────────────────

  subroutine writecoord(self,fname)
!*************************************************
!* subroutine writecoord                         *
!* is the write procedure for the "coord" class. *
!*************************************************
    implicit none
    class(coord) :: self
    character(len=*),intent(in) :: fname
    character(len=80) :: comment
    integer :: ftype,iunit
    if (.not.allocated(self%xyz)) then
      write (stdout,*) 'Cannot write ',trim(fname),'. No coordinates allocated'
    end if
    call checkcoordtype(fname,ftype)
    open (newunit=iunit,file=trim(fname))
    select case (ftype)
    case (coordtype%xyz)
      if (self%wrextxyz) then
        call self%writeextxyz(iunit)
      else
        call wrxyz(iunit,self%nat,self%at,self%xyz*autoaa,self%energy)
      end if

    case (coordtype%extxyz)
      call self%writeextxyz(iunit)

    case (coordtype%sdf,coordtype%sdfV3000)
      call wrsdfV3000(iunit,self%nat,self%at,self%xyz*autoaa, &
        & self%energy,real(self%chrg,wp),real(self%bond,wp),' written by CREST')

    case (coordtype%sdfV2000)
      call wrsdfV2000(iunit,self%nat,self%at,self%xyz*autoaa, &
        & self%energy,self%chrg,real(self%bond,wp),' written by CREST')

    case (coordtype%PDB)
      write (stdout,'(a)') '**ERROR** PDB file writer not implemented, TODO'
      call exit(1)
    case default
      !> defaults to Turbomole coord type
      call wrc0(iunit,self%nat,self%at,self%xyz)
    end select
    close (iunit)
    return
  end subroutine writecoord

! ──────────────────────────────────────────────────────────────────────────────

  subroutine appendcoord(self,iunit,energy)
!*************************************************
!* subroutine appendcoord                        *
!* is the write procedure for the "coord" class. *
!* coords will be written out in XYZ format!     *
!*************************************************
    implicit none
    class(coord) :: self
    integer,intent(in) :: iunit
    real(wp),intent(in),optional :: energy
    character(len=64) :: atmp
    character(len=32) :: btmp
    real(wp) :: etmp
    if (.not.self%wrextxyz) then !> regular xyz append
      self%xyz = self%xyz*bohr !to Angström
      if (present(energy)) then
        write (btmp,'(f22.10)') energy
      else
        write (btmp,'(f22.10)') self%energy
      end if
      write (atmp,'(a,a)') ' energy= ',adjustl(btmp)
      if (allocated(self%comment)) then
        call wrxyz(iunit,self%nat,self%at,self%xyz, &
        &          trim(atmp)//' '//trim(self%comment))
      else
        call wrxyz(iunit,self%nat,self%at,self%xyz,trim(atmp))
      end if
      self%xyz = self%xyz/bohr !back
    else
      !> extxyz append
      etmp = self%energy
      if (present(energy)) self%energy = energy
      call self%writeextxyz(iunit)
      self%energy = etmp
    end if
    return
  end subroutine appendcoord

! ══════════════════════════════════════════════════════════════════════════════
!  GENERAL UTILITY ROUTINES
! ══════════════════════════════════════════════════════════════════════════════

  subroutine atswp(self,ati,atj)
    !********************************
    !* swap atom ati with atj in mol
    !********************************
    implicit none
    class(coord),intent(inout) :: self
    integer,intent(in) :: ati,atj
    real(wp) :: xyztmp(3)
    integer :: attmp
    xyztmp(1:3) = self%xyz(1:3,ati)
    attmp = self%at(ati)
    self%xyz(1:3,ati) = self%xyz(1:3,atj)
    self%at(ati) = self%at(atj)
    self%xyz(1:3,atj) = xyztmp(1:3)
    self%at(atj) = attmp
  end subroutine atswp

! ──────────────────────────────────────────────────────────────────────────────

  function coord_sumform(self) result(sumformula)
    implicit none
    class(coord) :: self
    character(len=:),allocatable :: sumformula
    sumformula = sumform(self%nat,self%at)
  end function coord_sumform

! ══════════════════════════════════════════════════════════════════════════════
! end of the module
! ══════════════════════════════════════════════════════════════════════════════
end module molecule_type
