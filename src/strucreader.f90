!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020-2023 Philipp Pracht
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

!=========================================================================================!
! STRUCRD is a module for reading and writing molecular structures.
!
! The source is organized as follows:
!   0. Variable declarations
!   1. Routines for reading and writing ensemble files/trajectories in the XYZ format
!   2. Routines for reading single structures in various formats
!   3. Routines for writing structures in various formats
!   4. Utility routines mainly used only within the module
!
! Currently supported formats:
!   .xyz (Xmol) files and trajectories (read and write)
!   coord (turbomole) files (read and write)
!   .sdf/.mol files (V2000, read only)
!   .pdb files (in development)
!
!=========================================================================================!
module strucrd
  use iso_fortran_env,only:wp => real64
  use iso_c_binding
  use geo !> simple geomerty and vector operations
  use miscdata, only: PSE
  implicit none

!=========================================================================================!
!>--- private module variables and parameters
  private

!>--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp
!>-- filetypes as integers
  integer,parameter :: tmcoord = 1
  integer,parameter :: xmol = 2
  integer,parameter :: sdf = 3  !currently unused
  integer,parameter :: sdfV2000 = 31
  integer,parameter :: sdfV3000 = 32
  integer,parameter :: pdbfile = 4  !currently unused
  ! [...]

!>--- private utility subroutines
  private :: upperCase,lowerCase
  private :: convertlable,fextension,sgrep

!!>--- Element symbols
!!&<
!  character(len=2),private,parameter :: PSE(118) = [ &
! & 'H ',                                                                                'He', &
! & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
! & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
! & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
! & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
! & 'Cs','Ba','La',                                                                            &
! &                'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',      &
! &                'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
! & 'Fr','Ra','Ac',                                                                            &
! &                'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',      &
! &                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' ]
!!&>

!=========================================================================================!
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

  public :: xyz2coord
  public :: coord2xyz

  public :: rdensembleparam   !-- read Nat and Nall for a XYZ trajectory
  public :: rdensemble        !-- read a XYZ trajectory
  interface rdensemble
    module procedure rdensemble_conf1
    module procedure rdensemble_conf2
    module procedure rdensemble_conf3

    module procedure rdensemble_mixed2
    
    module procedure rdensemble_coord_type
  end interface rdensemble

  public :: wrensemble
  interface wrensemble
    module procedure wrensemble_conf
    module procedure wrensemble_conf_energy
    module procedure wrensemble_conf_energy_comment
  end interface wrensemble

  public :: pdbdata
  public :: coord
  public :: ensemble
  public :: coordline

!=========================================================================================!
  !coord class. contains a single structure in the PDB format.
  !coordinates by definition are in Angstroem.
  type :: pdbdata

    !--- data
    integer :: nat = 0
    integer :: frag = 0
    !--- arrays
    integer,allocatable  :: athet(:) !ATOM (1) or HETATM (2)
    character(len=4),allocatable :: pdbat(:) !PDB atom specifier
    character(len=3),allocatable :: pdbas(:) !PDB amino acid specifier
    integer,allocatable :: pdbfrag(:) !PDB fragment specifier
    character(len=1),allocatable :: pdbgrp(:)  !PDB group specifier
    real(wp),allocatable :: pdbocc(:) !PDB occupancy
    real(wp),allocatable :: pdbtf(:)  !PDB temperature factor

  contains
    procedure :: deallocate => deallocate_pdb !clear memory space
    procedure :: allocate => allocate_pdb

  end type pdbdata
!=========================================================================================!
  !coord class. contains a single structure
  !by convention coordinates are in atomic units (Bohr) for a single structure!
  type :: coord

    !********************************************!
    !> data that's typically used in coord type <!
    !********************************************!
    !>-- number of atoms
    integer :: nat = 0
    !>-- energy
    real(wp) :: energy = 0.0_wp
    !>-- atom types as integer, dimension will be at(nat)
    integer,allocatable  :: at(:)
    !>-- atomic coordinates, by convention in Bohrs
    real(wp),allocatable :: xyz(:,:)

    !**************************************!
    !> (optional) data, often not present <!
    !**************************************!
    !>-- a comment line
    character(len=:),allocatable :: comment
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

    !--- (optional) PDB data
    type(pdbdata) :: pdb

  contains
    procedure :: deallocate => deallocate_coord !> clear memory space
    procedure :: open => opencoord              !> read an coord file
    procedure :: write => writecoord            !> write
    procedure :: append => appendcoord          !> append
    procedure :: get => getcoord                !> allocate & fill with data
    procedure :: appendlog                      !> append .log file with coordinates and energy
    procedure :: dist => coord_getdistance      !> calculate distance between two atoms
    procedure :: angle => coord_getangle        !> calculate angle between three atoms
    procedure :: dihedral => coord_getdihedral  !> calculate dihedral angle between four atoms
  end type coord
!=========================================================================================!
  !ensemble class. contains all structures of an ensemble
  !by convention coordinates are in Angström for an ensemble!
  type :: ensemble

    !--- data
    integer :: nat = 0             !number of total atoms
    integer :: nall = 0             !number of structures
    integer,allocatable :: vnat(:)     !used instead of nat if not all structures have the same      number of atoms, in which case nat will be  =maxval(vnat,1)

    integer,allocatable  :: at(:)      !atom types as integer, dimension will be at(nat)
    real(wp),allocatable :: xyz(:,:,:) !coordinates, dimension will be xyz(3,nat,nall)
    real(wp),allocatable :: er(:)   !energy of each structure, dimension will be eread(nall)

    real(wp)            :: g         !gibbs free energy
    real(wp)            :: s         !entropy
    real(wp),allocatable :: gt(:)    !gibbs free energy of each member
    real(wp),allocatable :: ht(:)    !enthalpy of each member
    real(wp),allocatable :: svib(:)  !vibrational entropy of each member
    real(wp),allocatable :: srot(:)  !rotational entropy of each member
    real(wp),allocatable :: stra(:)  !translational entropy of each member

  contains
    procedure :: deallocate => deallocate_ensembletype !clear memory space
    procedure :: open => openensemble !read an ensemble file
    procedure :: write => write_ensemble !write to file

  end type ensemble
!=========================================================================================!
!=========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!=========================================================================================!
!=========================================================================================!
!  1. ROUTINES FOR READING ENTIRE ENSEMBLES (OR TRAJECTORIES)
!=========================================================================================!
!=========================================================================================!

!==================================================================!
! subroutine rdensembleparam
! read a ensemble file and get some information from
! it:
! On Input: fname - name of the file, should be in
!                   the Xmol (*.xyz) format.
!
! On Output: nat  - number of atoms
!                   (if different sized structures are present,
!                    nat is the largest)
!            nall - number of structures
!            conform - (optional) do all structures
!                      have the same number of atoms?
!=================================================================!
  subroutine rdensembleparam(fname,nat,nall,conform)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(out) :: nat
    integer,intent(out) :: nall
    logical,optional :: conform
    logical :: conformdum
    integer :: dum,iosum
    integer :: natref
    real(wp) :: x,y,z
    integer :: i,j,k,ich,io
    logical :: ex
    character(len=10) :: str
    conformdum = .true.
    nat = 0
    nall = 0
    natref = 0
    inquire (file=fname,exist=ex)
    if (.not.ex) return
    open (newunit=ich,file=fname)
    do
      read (ich,*,iostat=io) dum
      if (io < 0) exit
      if (io > 0) cycle
      if (nat == 0) natref = dum
      read (ich,*,iostat=io)
      if (io < 0) exit
      iosum = 0
      do i = 1,dum
        read (ich,*,iostat=io) str,x,y,z
        if (io < 0) exit
        iosum = iosum+io
      end do
      if (iosum > 0) cycle
      nat = max(dum,nat)
      if (dum .ne. natref) conformdum = .false.
      nall = nall+1
    end do
    close (ich)
    if (present(conform)) conform = conformdum
    return
  end subroutine rdensembleparam

!==================================================================!
! subroutine rdensemble_conf1
! read a conformer ensemble/a MD trajectory, i.e.,
! all structures have the same number and order of atoms.
! version 1 also reads the energy
!=================================================================!
  subroutine rdensemble_conf1(fname,nat,nall,at,xyz,eread)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(inout) :: nat
    integer,intent(inout) :: nall
    integer,intent(inout),allocatable :: at(:)
    real(wp),intent(inout),allocatable :: xyz(:,:,:)
    real(wp),intent(inout),allocatable :: eread(:)
    integer :: i,j,k,ich,io
    logical :: ex
    integer :: dum
    character(len=512) :: line
    character(len=6) :: sym
    if(.not.allocated(xyz).or. .not.allocated(at))then
      call rdensembleparam(fname,nat,nall)
    endif
    if(.not.allocated(xyz)) allocate(xyz(3,nat,nall))
    if(.not.allocated(at)) allocate(at(nat))
    if(.not.allocated(eread)) allocate(eread(nall))

    eread = 0.0_wp
    xyz = 0.0_wp
    open (newunit=ich,file=fname)
    do i = 1,nall
      read (ich,*,iostat=io) dum
      if (io < 0) exit
      if (io > 0) cycle
      if (dum .ne. nat) then
        call ensemble_strucskip(ich,nat,io)
        if (io < 0) exit
      end if
      read (ich,'(a)',iostat=io) line
      if (io < 0) exit
      eread(i) = grepenergy(line)
      do j = 1,dum
        read (ich,'(a)',iostat=io) line
        if (io < 0) exit
        call coordline(line,sym,xyz(1:3,j,i),io)
        if (io .ne. 0) then
          backspace (ich)
          exit
        end if
        at(j) = e2i(sym)
      end do
    end do
    close (ich)

    if (io < 0) then
      error stop 'error while reading ensemble file.'
    end if

    return
  end subroutine rdensemble_conf1

!==================================================================!
! subroutine rdensemble_conf2
! read a conformer ensemble/a MD trajectory, i.e.,
! all structures have the same number and order of atoms.
! version 2 does not read the energy
!=================================================================!
  subroutine rdensemble_conf2(fname,nat,nall,at,xyz)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(inout) :: nat
    integer,intent(inout) :: nall
    integer,intent(inout),allocatable :: at(:)
    real(wp),intent(inout),allocatable :: xyz(:,:,:)
    integer :: i,j,k,ich,io
    logical :: ex
    integer :: dum,nallnew
    character(len=512) :: line
    character(len=6) :: sym
    if(.not.allocated(xyz).or. .not.allocated(at))then
      call rdensembleparam(fname,nat,nall)
    endif
    if(.not.allocated(xyz)) allocate(xyz(3,nat,nall))
    if(.not.allocated(at)) allocate(at(nat))
    io = 0
    xyz = 0.0_wp
    open (newunit=ich,file=fname)
    do i = 1,nall
      read (ich,*,iostat=io) dum
      if (io < 0) exit
      if (io > 0) cycle
      if (dum .ne. nat) then
        call ensemble_strucskip(ich,nat,io)
        if (io < 0) exit
      end if
      read (ich,'(a)',iostat=io) line
      if (io < 0) exit
      do j = 1,dum
        read (ich,'(a)',iostat=io) line
        if (io < 0) exit
        call coordline(line,sym,xyz(1:3,j,i),io)
        if (io .ne. 0) then
          backspace (ich)
          exit
        end if
        at(j) = e2i(sym)
      end do
    end do
    close (ich)

    if (io < 0) then
      error stop 'error while reading ensemble file.'
    end if

    return
  end subroutine rdensemble_conf2

!==================================================================!
! subroutine rdensemble_conf3
! read a conformer ensemble/a MD trajectory, i.e.,
! all structures have the same number and order of atoms.
! version 3 saves the comment line for each structure
!=================================================================!
  subroutine rdensemble_conf3(fname,nat,nall,at,xyz,comments)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(inout) :: nat
    integer,intent(inout) :: nall
    integer :: at(nat)
    integer,allocatable :: atdum(:)
    real(wp) :: xyz(3,nat,nall)
    character(len=*) :: comments(nall)
    integer :: i,j,k,ich,io
    logical :: ex
    integer :: dum,nallnew
    character(len=512) :: line
    character(len=6) :: sym
    io = 0
    xyz = 0.0_wp
    k = 0
    open (newunit=ich,file=fname)
    do i = 1,nall
      read (ich,*,iostat=io) dum
      if (io < 0) exit
      if (io > 0) cycle
      if (dum .ne. nat) then
        call ensemble_strucskip(ich,nat,io)
        if (io < 0) exit
      end if
      read (ich,'(a)',iostat=io) line
      if (io < 0) exit
      comments(i) = trim(line)
      do j = 1,dum
        k = k+1
        read (ich,'(a)',iostat=io) line
        if (io < 0) exit
        call coordline(line,sym,xyz(1:3,j,i),io)
        if (io .ne. 0) then
          backspace (ich)
          exit
        end if
        at(j) = e2i(sym)
      end do
    end do
    close (ich)

    if (io < 0) then
      error stop 'error while reading ensemble file.'
    end if

    return
  end subroutine rdensemble_conf3

  subroutine ensemble_strucskip(ich,nat,io)
    implicit none
    integer,intent(in) :: ich
    integer,intent(in) :: nat
    integer,intent(out) :: io
    integer :: io2,dum,k
    io = 0
    dum = 0
    k = 0
    do while (dum .ne. nat)
      read (ich,*,iostat=io) dum
      if (io < 0) exit
      k = k+1
      if (io > 0) cycle
    end do
  end subroutine ensemble_strucskip

!==================================================================!
! subroutine rdensemble_mixed2
! read an ensemble of mixed strcutres, i.e., all stuctures
! can have a diferent number and order of atoms.
! version 2 does not read energies
!=================================================================!
  subroutine rdensemble_mixed2(fname,natmax,nall,nats,ats,xyz)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: natmax
    integer,intent(in) :: nall
    integer  :: nats(nall)
    integer  :: ats(natmax,nall)
    real(wp) :: xyz(3,natmax,nall)
    integer :: i,j,k,ich,io
    logical :: ex
    integer :: dum
    character(len=512) :: line
    character(len=6) :: sym
    open (newunit=ich,file=fname)
    do i = 1,nall
      read (ich,*,iostat=io) dum
      if (io < 0) exit
      if (io > 0) cycle
      nats(i) = dum
      read (ich,'(a)',iostat=io) line
      if (io < 0) exit
      do j = 1,dum
        read (ich,'(a)',iostat=io) line
        if (io < 0) exit
        call coordline(line,sym,xyz(1:3,j,i),io)
        if (io < 0) exit
        ats(j,i) = e2i(sym)
      end do
    end do
    close (ich)

    if (io < 0) then
      error stop 'error while reading ensemble file.'
    end if

    return
  end subroutine rdensemble_mixed2

!========================================================================================!
  subroutine rdensemble_coord_type(fname,nall,ensemble)
!*********************************************************
!* subroutine rdensemble_coord_type
!* A variant of the rdensemble routine that automatically
!* produces an array of coord containers
!*********************************************************
    implicit none
    character(len=*),intent(in) :: fname !> name of the ensemble file
    integer,intent(out) :: nall  !> number of structures in ensemble
    type(coord),intent(out),allocatable :: ensemble(:)

    real(wp),allocatable :: xyz(:,:,:)
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: eread(:)
    integer :: i,j,k,ich,io
    logical :: ex

    call rdensembleparam(fname,nat,nall)
    allocate (xyz(3,nat,nall),at(nat),eread(nall))
    call rdensemble(fname,nat,nall,at,xyz,eread)
    !>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<!
    !>--- Important: coord types must be in Bohrs
    xyz = xyz/bohr
    !>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<!
 
    allocate (ensemble(nall)) 
    do i=1,nall
      ensemble(i)%nat = nat 
      allocate(ensemble(i)%at(nat))
      ensemble(i)%at(:) = at(:)
      allocate(ensemble(i)%xyz(3,nat))
      ensemble(i)%xyz(:,:) = xyz(:,:,i)
      ensemble(i)%energy = eread(i)
    enddo

    deallocate(eread,at,xyz)
  end subroutine rdensemble_coord_type

!=================================================================!
! subroutine wrensemble_conf
! write a ensemble file/a trajectory from memory.
!=================================================================!
  subroutine wrensemble_conf(fname,nat,nall,at,xyz)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(in) :: nall
    integer :: at(nat)
    real(wp) :: xyz(3,nat,nall)
    integer :: i,j,k,ich,io
    logical :: ex

    open (newunit=ich,file=fname,status='replace')
    do i = 1,nall
      call wrxyz(ich,nat,at,xyz(:,:,i))
    end do
    close (ich)

    return
  end subroutine wrensemble_conf

!=================================================================!
! subroutine wrensemble_conf_energy
! write a ensemble file/a trajectory from memory.
!=================================================================!
  subroutine wrensemble_conf_energy(fname,nat,nall,at,xyz,er)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(in) :: nall
    integer :: at(nat)
    real(wp) :: xyz(3,nat,nall)
    real(wp) :: er(nall)
    integer :: i,j,k,ich,io
    logical :: ex

    open (newunit=ich,file=fname,status='replace')
    do i = 1,nall
      call wrxyz(ich,nat,at,xyz(:,:,i),er(i))
    end do
    close (ich)

    return
  end subroutine wrensemble_conf_energy

!=================================================================!
! subroutine wrensemble_conf_energy_comment
! write a ensemble file/a trajectory from memory.
!=================================================================!
  subroutine wrensemble_conf_energy_comment(fname,nat,nall,at,xyz,er,comments)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(in) :: nall
    integer :: at(nat)
    real(wp) :: xyz(3,nat,nall)
    real(wp) :: er(nall)
    character(len=*) :: comments(nall)
    character(len=512) :: line
    integer :: i,j,k,ich,io
    logical :: ex

    open (newunit=ich,file=fname,status='replace')
    do i = 1,nall
      write (line,'(2x,f18.8,2x,a)') er(i),trim(comments(i))
      call wrxyz(ich,nat,at,xyz(:,:,i),trim(line))
    end do
    close (ich)

    return
  end subroutine wrensemble_conf_energy_comment

!==================================================================!
! subroutine write_ensemble
! wrapper to write an ensemble from the "ensemble" class
!==================================================================!
  subroutine write_ensemble(self,fname)
    implicit none
    class(ensemble) :: self
    character(len=*),intent(in) :: fname
    call wrensemble_conf_energy(fname,self%nat,self%nall,self%at,self%xyz,self%er)
    return
  end subroutine write_ensemble

!==================================================================!
! subroutine deallocate_ensembletype
! is used to clear memory for the ensemble type
!==================================================================!
  subroutine deallocate_ensembletype(self)
    implicit none
    class(ensemble) :: self
    self%nat = 0
    self%nall = 0
    if (allocated(self%vnat)) deallocate (self%vnat)
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%er)) deallocate (self%er)
    if (allocated(self%gt)) deallocate (self%gt)
    if (allocated(self%ht)) deallocate (self%ht)
    if (allocated(self%svib)) deallocate (self%svib)
    if (allocated(self%srot)) deallocate (self%srot)
    if (allocated(self%stra)) deallocate (self%stra)

    return
  end subroutine deallocate_ensembletype

!==================================================================!
! subroutine openensemble
! is the open procedure for the "ensemble" class.
! a ensemble (trajectory) fname is read into a new ensemble object
!==================================================================!
  subroutine openensemble(self,fname)
    implicit none
    class(ensemble) :: self
    character(len=*),intent(in) :: fname
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:,:)
    real(wp),allocatable :: eread(:)
    integer :: nall
    integer :: i,j,k,ich,io
    logical :: ex

    inquire (file=fname,exist=ex)
    if (.not.ex) then
      error stop 'ensemble file does not exist.'
    end if

    call rdensembleparam(fname,nat,nall)

    if (nat > 0.and.nall > 0) then
      call self%deallocate()
      allocate (at(nat),xyz(3,nat,nall),eread(nall))
      call rdensemble(fname,nat,nall,at,xyz,eread)

      self%nat = nat
      self%nall = nall
      call move_alloc(at,self%at)
      call move_alloc(xyz,self%xyz)
      call move_alloc(eread,self%er)
    else
      error stop 'format error while reading ensemble file.'
    end if

    return
  end subroutine openensemble

!=========================================================================================!
!=========================================================================================!
!  2. ROUTINES FOR READING SINGLE STRUCTURES (COORDS)
!=========================================================================================!
!=========================================================================================!

!============================================================!
! subroutine checkcoordtype
! try to identify the filetype of the coord type.
! first based on file extension, if that fails by
! a keyword within the file.
!============================================================!
  subroutine checkcoordtype(fname,typint)
    implicit none
    character(len=*) :: fname
    integer,intent(out) :: typint
    typint = 0
    !-- check file extension first
    select case (fextension(fname))
    case ('.coord','.COORD')
      typint = tmcoord
    case ('.xyz','.XYZ','.trj','.TRJ','.sorted')
      typint = xmol
    case ('.sd','.sdf','.SDF','.mol','.MOL')
      typint = sdf
      if (sgrep(fname,'V2000')) then
        typint = sdfV2000
      end if
      if (sgrep(fname,'V3000')) then
        typint = sdfV3000
      end if
    case ('.pdb','.PDB')
      typint = pdbfile
    case default
      typint = 0
    end select
    if (typint .ne. 0) return !-- file extension was recognized
    !-- grep for keywords otherwise
    if (sgrep(fname,'$coord')) then
      typint = tmcoord
    else !--no match found
      typint = 0
    end if
    return
  end subroutine checkcoordtype

!============================================================!
! subroutine rdnat
! read number of atoms "nat" form file
!
! On Input: fname  - name of the coord file
!           ftype  - (OPTIONAL) format of the input coord file
!                    if ftype is not present, it is determined
! On Output: nat   - number of atoms
!============================================================!
  subroutine rdnat(fname,nat,ftype)
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
      error stop 'file does not exist.'
    end if
    if (present(ftype)) then
      ftypedum = ftype
    else
      call checkcoordtype(fname,ftypedum)
    end if
    open (newunit=ich,file=fname)
    select case (ftypedum)
      !--- *.xyz files
    case (xmol)
      read (ich,*,iostat=io) nat
      !--- TM coord file
    case (tmcoord)
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
      !--- sdf V2000 (or *.mol) file
    case (sdfV2000)
      do i = 1,3 !-- first three comment lines
        read (ich,'(a)',iostat=io) atmp
        if (io < 0) exit
      end do
      read (ich,'(a)',iostat=io) atmp
      if (index(atmp,'V2000') .ne. 0) then
        read (atmp,'(i3)') nat !- first argument is nat
      end if
      !--- sdf V3000 file
    case (sdfV3000)
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
      !--- pdb file
    case (pdbfile)
      !write(*,*) 'PDB file format not supported yet.'
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

!============================================================!
! subroutine rdcoord
! read in a structure. The format is determined automatically
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (always in Bohr)
!            energy - (OPTIONAL) if present, try to get energy
!                      mainly from xyz files
!============================================================!
  subroutine rdcoord(fname,nat,at,xyz,energy)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    real(wp),optional :: energy
    character(len=256) :: atmp
    integer :: ftype
    type(pdbdata) :: pdbdummy

    !--- determine the file type
    call checkcoordtype(fname,ftype)

    select case (ftype)
    case (tmcoord)  !-- TM coord file, is already in Bohr
      call rdtmcoord(fname,nat,at,xyz)
    case (xmol)     !-- XYZ file, is Angström, needs conversion
      if (present(energy)) then
        call rdxmol(fname,nat,at,xyz,atmp)
        energy = grepenergy(atmp)
      else
        call rdxmol(fname,nat,at,xyz)
      end if
      xyz = xyz/bohr
    case (sdfV2000)      !-- SDF/MOL V2000 file, also Angström
      call rdsdf(fname,nat,at,xyz)
      xyz = xyz/bohr
    case (sdfV3000)     !-- SDF V3000 file, Angström
      call rdsdfV3000(fname,nat,at,xyz)
      xyz = xyz/bohr
    case (pdbfile)  !-- PDB file, Angström
      !error stop 'PDB file format not supported yet.'
      call rdPDB(fname,nat,at,xyz,pdbdummy)
      xyz = xyz/bohr
      call pdbdummy%deallocate()
    case default
      continue
    end select

    return
  end subroutine rdcoord

!============================================================!
! subroutine rdtmcoord
! read a struncture in the TM coord style.
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (always in Bohr)
!============================================================!
  subroutine rdtmcoord(fname,nat,at,xyz)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nat
    integer,intent(inout)  :: at(nat)
    real(wp),intent(inout) :: xyz(3,nat)
    character(len=6) :: sym
    integer :: ich,io,i
    character(len=256) :: atmp
    open (newunit=ich,file=fname)
    do
      read (ich,'(a)',iostat=io) atmp
      if (io < 0) exit
      atmp = adjustl(atmp)
      if (index(atmp,"$coord") .eq. 1) exit
    end do
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
    return
  end subroutine rdtmcoord

!============================================================!
! subroutine rdxmol
! read a struncture in the *.xyz (Xmol) style.
! The commentary (second) line is ignored
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            comment - (OPTIONAL) commentary line of the file
!============================================================!
  subroutine rdxmol(fname,nat,at,xyz,comment)
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
      error stop 'error while reading input coordinates'
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
    close (ich)
    return
  end subroutine rdxmol

!============================================================!
! subroutine rdsdf
! read a struncture in the .sdf/.mol V2000 style.
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            comment - (OPTIONAL) commentary line of the file
!============================================================!
  subroutine rdsdf(fname,nat,at,xyz,comment)
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
      error stop 'error while reading input coordinates'
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

!============================================================!
! subroutine rdsdfV3000
! read a struncture in the .sdf/.mol V3000 style.
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            comment - (OPTIONAL) commentary line of the file
!============================================================!
  subroutine rdsdfV3000(fname,nat,at,xyz,comment)
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
      error stop 'error while reading input coordinates'
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

!============================================================!
! subroutine rdPDB
! read a struncture in the .PDB style.
!
! On Input: fname  - name of the coord file
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!            pdb  - pdbdata object
!============================================================!
  subroutine rdPDB(fname,nat,at,xyz,pdb)
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

!============================================================!
! subroutine rdxmolselec
! Read a file with multiple structures in the *.xyz (Xmol) style.
! Picks one structure.
! The commentary (second) line is ignored
!
! On Input: fname  - name of the coord file
!           m      - position of the desired structure
!           nat    - number of atoms
!
! On Output: at   - atom number as integer
!            xyz  - coordinates (in Angström)
!============================================================!

  subroutine rdxmolselec(fname,m,nat,at,xyz,comment)
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
        error stop 'error while reading input coordinates'
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

!==================================================================!
! subroutine deallocate_coord
! is used to clear memory for the coord type
!==================================================================!
  subroutine deallocate_coord(self)
    implicit none
    class(coord) :: self
    self%nat = 0
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    call self%pdb%deallocate()
    return
  end subroutine deallocate_coord

!==================================================================!
! subroutine deallocate_pdb
! is used to clear memory for the pdbdata type
!==================================================================!
  subroutine deallocate_pdb(self)
    implicit none
    class(pdbdata) :: self
    self%nat = 0
    self%frag = 0
    if (allocated(self%athet)) deallocate (self%athet)
    if (allocated(self%pdbat)) deallocate (self%pdbat)
    if (allocated(self%pdbas)) deallocate (self%pdbas)
    if (allocated(self%pdbfrag)) deallocate (self%pdbfrag)
    if (allocated(self%pdbgrp)) deallocate (self%pdbgrp)
    if (allocated(self%pdbocc)) deallocate (self%pdbocc)
    if (allocated(self%pdbtf)) deallocate (self%pdbtf)
    return
  end subroutine deallocate_pdb

!==================================================================!
! subroutine allocate_pdb
! is used to clear memory for the pdbdata type
!==================================================================!
  subroutine allocate_pdb(self,nat)
    implicit none
    class(pdbdata) :: self
    integer :: nat
    call deallocate_pdb(self)
    self%nat = nat
    allocate (self%athet(nat))
    allocate (self%pdbat(nat))
    allocate (self%pdbas(nat))
    allocate (self%pdbfrag(nat))
    allocate (self%pdbgrp(nat))
    allocate (self%pdbocc(nat))
    allocate (self%pdbtf(nat))
    return
  end subroutine allocate_pdb

!==================================================================!
! subroutine opencoord
! is the open procedure for the "coord" class.
!==================================================================!
  subroutine opencoord(self,fname)
    implicit none
    class(coord) :: self
    character(len=*),intent(in) :: fname
    integer :: nat
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    integer :: ftype
    integer :: i,j,k,ich,io
    logical :: ex

    inquire (file=fname,exist=ex)
    if (.not.ex) then
      error stop 'coord file does not exist.'
    end if

    call self%deallocate()

    call checkcoordtype(fname,ftype)
    call rdnat(fname,nat)

    if (nat > 0) then
      allocate (at(nat),xyz(3,nat))
      if (ftype == pdbfile) then
        call rdPDB(fname,nat,at,xyz,self%pdb)
        xyz = xyz/bohr
      else
        call rdcoord(fname,nat,at,xyz)
      end if

      self%nat = nat
      call move_alloc(at,self%at)
      call move_alloc(xyz,self%xyz)
    else
      error stop 'format error while reading coord file.'
    end if

    return
  end subroutine opencoord
!==================================================================!
! subroutine getcoord
! allocate "coord" class and fill with data
!==================================================================!
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

!==================================================================!
! function coord_getdistance
! calculate the distance for a given pair of atoms
!==================================================================!
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

!==================================================================!
! function coord_getangle
! calculate the angle for a given trio of atoms
! A1-A2-A3
!==================================================================!
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

!==================================================================!
! function coord_getdihedral
! calculate the dihedral angle for a given quartet of atoms
! A1-A2-A3-A4
!==================================================================!
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
    write (ch,'(2x,i0)') nat
    write (ch,'(2x,f18.8)') er
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

!============================================================!
! subroutine xyz2coord
! simple conversion of a xyz to a coord file.
!
! On Input: iname  - name of the xyz file
!           oname  - name of the coord file
!
! On Output: file written to "oname"
!============================================================!
  subroutine xyz2coord(iname,oname)
    implicit none
    character(len=*) :: iname
    character(len=*) :: oname
    type(coord) :: struc
    call struc%open(iname)
    call wrc0(oname,struc%nat,struc%at,struc%xyz)
    call struc%deallocate()
    return
  end subroutine xyz2coord

!============================================================!
! subroutine coord2xyz
! simple conversion of a coord to a xyz file.
!
! On Input: iname  - name of the coord file
!           oname  - name of the xyz file
!
! On Output: file written to "oname"
!============================================================!
  subroutine coord2xyz(iname,oname)
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

!==================================================================!
! subroutine writecoord
! is the write procedure for the "coord" class.
!==================================================================!
  subroutine writecoord(self,fname)
    implicit none
    class(coord) :: self
    character(len=*),intent(in) :: fname
    if (.not.allocated(self%xyz)) then
      write (*,*) 'Cannot write ',trim(fname),'. not allocated'
    end if
    if (index(fname,'.xyz') .ne. 0) then
      self%xyz = self%xyz*bohr !to Angström
      call wrxyz(fname,self%nat,self%at,self%xyz)
      self%xyz = self%xyz/bohr !back
    else
      call wrc0(fname,self%nat,self%at,self%xyz)
    end if
    return
  end subroutine writecoord

!==================================================================!
! subroutine appendcoord
! is the write procedure for the "coord" class.
! coords will be written out in XYZ format!
!==================================================================!
  subroutine appendcoord(self,io)
    implicit none
    class(coord) :: self
    integer :: io
    character(len=64) :: atmp
    self%xyz = self%xyz*bohr !to Angström
    if (allocated(self%comment)) then
      call wrxyz(io,self%nat,self%at,self%xyz,trim(self%comment))
    else if (self%energy .ne. 0.0_wp) then
      write (atmp,'(a,f22.10)') ' Etot= ',self%energy
      call wrxyz(io,self%nat,self%at,self%xyz,trim(atmp))
    else
      call wrxyz(io,self%nat,self%at,self%xyz)
    end if
    self%xyz = self%xyz/bohr !back
    return
  end subroutine appendcoord

  subroutine appendlog(self,io,energy,gnorm)
    implicit none
    class(coord) :: self
    integer :: io
    real(wp),optional :: energy
    real(wp),optional :: gnorm
    character(len=64) :: atmp
    self%xyz = self%xyz*bohr !to Angström
    if (present(gnorm).and.present(energy)) then
      write (atmp,'(a,f22.10,a,f16.8)') ' Etot= ',energy,' grad.norm.= ',gnorm
    else if (present(energy)) then
      write (atmp,'(a,f22.10)') ' Etot= ',energy
    else
      atmp = ''
    end if
    call wrxyz(io,self%nat,self%at,self%xyz,trim(atmp))
    self%xyz = self%xyz/bohr !back
    return
  end subroutine appendlog

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
      !if(io.ne.0)then
      !  error stop 'error while reading coord line'
      !endif
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
    integer :: i,j,k,ich,io
    logical :: ex
    c = trim(convertlable(cin))
    if (any(PSE(:) .eq. c)) then
      do i = 1,118
        if (trim(PSE(i)) .eq. c) then
          iout = i
          exit
        end if
      end do
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
  function sgrep(fname,key)
    implicit none
    character(len=*),intent(in) :: fname
    character(len=*),intent(in) :: key
    logical :: sgrep
    character(len=256) :: atmp
    integer :: ic,io
    sgrep = .false.
    open (newunit=ic,file=fname)
    do
      read (ic,'(a)',iostat=io) atmp
      if (io < 0) exit !EOF
      if (index(atmp,key) .ne. 0) then
        sgrep = .true.
        exit
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
    integer :: i,io
    atmp = trim(line)
    energy = 0.0_wp
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
    grepenergy = energy
    return
  end function grepenergy

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
!=========================================================================================!
! end of the module
!=========================================================================================!
!=========================================================================================!
end module strucrd
