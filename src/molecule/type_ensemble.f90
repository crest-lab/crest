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



! ══════════════════════════════════════════════════════════════════════════════
!
! NOTE: While there are the two types "ensemble" and "mollist", the best way 
!       to handle a list of structures is as "type(coord),allocatable :: structures(:)"
!       which is mapped in the rdensemble_coord_type routine to rdensemble.
!       Both of the other types wrap that, but it would be better to not use them.
!
! ──────────────────────────────────────────────────────────────────────────────

module molecule_type_ensemble
  use iso_c_binding
  use molecule_parameters
  use molecule_io
  use molecule_type
  use molecule_type_components
  implicit none

! ══════════════════════════════════════════════════════════════════════════════
!>--- private module variables and parameters
  private

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

    module procedure wrensemble_coord_name
    module procedure wrensemble_coord_channel
  end interface wrensemble

  public :: ensemble
  public :: mollist

! ──────────────────────────────────────────────────────────────────────────────
  !> ensemble class. contains all structures of an ensemble
  !> by convention coordinates are in Angström for an ensemble!
  type :: ensemble

    logical :: mixed = .false.   !> if all molecules were the same == .false.

    !> data
    integer :: nat = 0              !> (max) number of total atoms
    integer :: nall = 0             !> number of structures

    !> if all structures were the same molecule these are filled
    !> mixed==.false.
    integer,allocatable  :: at(:)      !> atom types as integer, dimension will be at(nat)
    real(wp),allocatable :: xyz(:,:,:) !> coordinates, dimension will be xyz(3,nat,nall)
    real(wp),allocatable :: er(:)      !> energy of each structure, dimension will be eread(nall)

    !> otherwise this is filled
    !> mixed == .true.
    type(coord),allocatable :: structures(:)

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
    procedure :: get_mol => ensemble_get_mol !extract the i-th mol from ensemble type
  end type ensemble

!==========================================================================================!
  type :: mollist
    integer :: nall = 0
    type(coord),allocatable :: structure(:)
  end type mollist

! ══════════════════════════════════════════════════════════════════════════════
contains  !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════
!  1. ROUTINES FOR READING ENTIRE ENSEMBLES (OR TRAJECTORIES)
! ──────────────────────────────────────────────────────────────────────────────

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
    logical,intent(out),optional :: conform
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
    if (.not.allocated(xyz).or..not.allocated(at)) then
      call rdensembleparam(fname,nat,nall)
    end if
    if (.not.allocated(xyz)) allocate (xyz(3,nat,nall))
    if (.not.allocated(at)) allocate (at(nat))
    if (.not.allocated(eread)) allocate (eread(nall))

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
    if (.not.allocated(xyz).or..not.allocated(at)) then
      call rdensembleparam(fname,nat,nall)
    end if
    if (.not.allocated(xyz)) allocate (xyz(3,nat,nall))
    if (.not.allocated(at)) allocate (at(nat))
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
  subroutine rdensemble_mixed2(fname,natmax,nall,nats,ats,xyz,comments)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: natmax
    integer,intent(in) :: nall
    integer  :: nats(nall)
    integer  :: ats(natmax,nall)
    real(wp) :: xyz(3,natmax,nall)
    character(len=*) :: comments(nall)
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
      comments(i) = trim(line)
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

! ══════════════════════════════════════════════════════════════════════════════
!> most important routine in this module: read into a list of coord objects!
! ══════════════════════════════════════════════════════════════════════════════
  subroutine rdensemble_coord_type(fname,nall,structures)
!**********************************************************
!* subroutine rdensemble_coord_type                       *
!* A variant of the rdensemble routine that automatically *
!* produces an array of coord containers                  *
!**********************************************************
    implicit none
    character(len=*),intent(in) :: fname !> name of the ensemble file
    integer,intent(out) :: nall  !> number of structures in ensemble
    type(coord),intent(out),allocatable :: structures(:)

    real(wp),allocatable :: xyz(:,:,:)
    integer :: nat
    integer,allocatable :: nats(:)
    integer,allocatable :: at(:)
    integer,allocatable :: ats(:,:)
    real(wp),allocatable :: eread(:)
    character(len=512),allocatable :: comments(:)
    integer :: i,j,k,ich,io,nat_i,iunit,ii
    logical :: ex,multiple_sizes,is_extxyz,success

    type(extxyz_signatures) :: ext_sigs
    type(extxyz_properties) :: ext_props
    real(wp),allocatable :: exyz(:,:),egrd(:,:),lat(:,:)
    integer,allocatable :: eat(:)
    real(wp) :: energy
    character(len=32) :: eu,fu

    is_extxyz = sgrep(fname,'Properties=',casesensitive=.false.)

    call rdensembleparam(fname,nat,nall,multiple_sizes)

    !>--- multiple sizes possible
    allocate (structures(nall))

    if (is_extxyz) then
      !>-- extended xyz case
      open (newunit=iunit,file=trim(fname))
      do ii = 1,nall
        call read_extxyz_frame(iunit,ext_sigs,ext_props,nat,energy,lat,success, &
        &                      energy_units=eu,forces_units=fu)
        if (success) then
          select case (trim(eu))
          case ('hartree','ha','au')
            ! energy already in Hartree, no conversion needed
          case default  !> 'ev' and anything unrecognised
            energy = energy/autoeV
          end select
          call get_at_from_ext(ext_props,eat)
          call get_xyz_from_ext(ext_props,exyz)
          call get_grad_from_ext(ext_props,egrd,forces_units=fu)
          if (allocated(eat)) call move_alloc(eat,structures(ii)%at)
          if(allocated(exyz)) call move_alloc(exyz,structures(ii)%xyz)
          if (allocated(lat)) call move_alloc(lat,structures(ii)%lat)
          if (allocated(egrd)) call move_alloc(egrd,structures(ii)%gradient)
          structures(ii)%energy = energy
          structures(ii)%wrextxyz = .true.
          structures(ii)%nat = nat
        end if
      end do
      close (iunit)
    else
      !>-- regular xyz case
      allocate (xyz(3,nat,nall),ats(nat,nall),nats(nall),eread(nall))
      allocate (comments(nall))
      call rdensemble_mixed2(fname,nat,nall,nats,ats,xyz,comments)
      !>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<!
      !>--- Important: coord types must be in Bohrs
      xyz = xyz/bohr
      !>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<!
      do i = 1,nall
        nat_i = nats(i)
        structures(i)%nat = nats(i)
        allocate (structures(i)%at(nat_i))
        structures(i)%at(:) = ats(1:nat_i,i)
        allocate (structures(i)%xyz(3,nat_i))
        structures(i)%xyz(:,:) = xyz(1:3,1:nat_i,i)
        eread(i) = grepenergy(comments(i))
        structures(i)%energy = eread(i)
        structures(i)%comment = trim(comments(i))
      end do
      deallocate (comments,eread,nats,ats,xyz)
    end if
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
    if (.not.self%mixed) then
      call wrensemble_conf_energy(fname,self%nat,self%nall,self%at,self%xyz,self%er)
    else
      self%structures(:)%energy = self%er(:)
      call wrensemble_coord_name(fname,self%nall,self%structures)
    end if
    return
  end subroutine write_ensemble

  subroutine wrensemble_coord_name(fname,nall,structures)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nall
    type(coord) :: structures(nall)
    integer :: ich,i
    open (newunit=ich,file=fname,status='replace')
    do i = 1,nall
      call structures(i)%append(ich)
    end do
    close (ich)
    return
  end subroutine wrensemble_coord_name

  subroutine wrensemble_coord_channel(ich,nall,structures)
    implicit none
    integer,intent(in) :: ich
    integer,intent(in) :: nall
    type(coord) :: structures(nall)
    integer :: i
    do i = 1,nall
      call structures(i)%append(ich)
    end do
    return
  end subroutine wrensemble_coord_channel

!==================================================================!
! subroutine deallocate_ensembletype
! is used to clear memory for the ensemble type
!==================================================================!
  subroutine deallocate_ensembletype(self)
    implicit none
    class(ensemble) :: self

    self%mixed = .false.
    self%nat = 0
    self%nall = 0
    if (allocated(self%at)) deallocate (self%at)
    if (allocated(self%xyz)) deallocate (self%xyz)
    if (allocated(self%er)) deallocate (self%er)

    if (allocated(self%structures)) deallocate (self%structures)

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
    logical :: ex,conform
    type(coord),allocatable :: structures(:)

    inquire (file=fname,exist=ex)
    if (.not.ex) then
      error stop 'ensemble file does not exist.'
    end if

    !> we check if all the structures in the file
    !> are actually the same length (nat), if not we need to
    !> take care of this and read into self%structures instead
    call rdensembleparam(fname,nat,nall,conform)
    self%mixed = .not.conform

    if (conform) then
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
    else
      call rdensemble_coord_type(fname,self%nall,self%structures)
      allocate (self%er(nall),source=0.0_wp)
      self%er(:) = self%structures(:)%energy
    end if

    return
  end subroutine openensemble

  subroutine ensemble_get_mol(self,i,mol)
    class(ensemble) :: self
    integer,intent(in) :: i
    class(coord),intent(inout) :: mol
    integer :: n
    logical :: reinitialize
    if (i > self%nall) error stop 'can´t get molecule from ensemble. i>nall'
    if (i < 1) error stop 'can´t get molecule from ensemble. i<1'
    if (.not.self%mixed) then
      n = self%nat
      reinitialize = .not. (mol%nat == n)
      if (reinitialize) then
        mol%nat = n
        if (allocated(mol%at)) deallocate (mol%at)
        allocate (mol%at(n),source=0)
        if (allocated(mol%xyz)) deallocate (mol%xyz)
        allocate (mol%xyz(3,n),source=0.0_wp)
      end if
      mol%energy = self%er(i)
      mol%at(:) = self%at(:)
      !> Important, ens is in Angström, mol is in Bohrs
      mol%xyz(1:3,1:n) = self%xyz(1:3,1:n,i)*aatoau
    else !> self%mixed == .true.
      n = self%structures(i)%nat
      reinitialize = .not. (mol%nat == n)
      if (reinitialize) then
        if (allocated(mol%at)) deallocate (mol%at)
        allocate (mol%at(n),source=0)
        if (allocated(mol%xyz)) deallocate (mol%xyz)
        allocate (mol%xyz(3,n),source=0.0_wp)
      end if
      mol%nat = self%structures(i)%nat
      mol%at(:) = self%structures(i)%at(:)
      mol%xyz(:,:) = self%structures(i)%xyz(:,:)
      mol%energy = self%structures(i)%energy
    end if
  end subroutine ensemble_get_mol

! ══════════════════════════════════════════════════════════════════════════════
! ══════════════════════════════════════════════════════════════════════════════
end module molecule_type_ensemble
