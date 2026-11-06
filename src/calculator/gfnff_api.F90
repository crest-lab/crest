!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023-2026 Philipp Pracht
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

!====================================================!
! module gfnff_api
! An interface to GFN-FF standalone calculations
!====================================================!

module gfnff_api
  use iso_fortran_env,only:wp => real64,stdout => output_unit,stderr => error_unit
  use strucrd
#ifdef WITH_GFNFF
  use gfnff_interface
#endif
  implicit none
  private

  !> link to subproject routines/datatypes
  public :: gfnff_data
  !> routines from this file
  public :: gfnff_api_setup
  public :: gfnff_sp
  public :: gfnff_printout
  public :: gfnff_getwbos
  public :: gfnff_dump_sasa
  public :: gfnff_set_fragments
  public :: gfnff_set_refcharges
  public :: gfnff_set_bondmat

#ifndef WITH_GFNFF
  !> these are placeholders if no gfnff module is used!
  type :: gfnff_data
    integer :: id = 0
    character(len=:),allocatable :: parametrisation
    logical :: restart = .false.
    character(len=:),allocatable :: restartfile
    character(len=:),allocatable :: refgeo
  end type gfnff_data
#endif

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine gfnff_api_setup(mol,chrg,ff_dat,io,pr,iunit)
!*************************************************************
!* Set up (initialize) a GFN-FF calculator from a coord mol. *
!* Lattice vectors are read from mol%lat when present, so    *
!* PBC calculations are automatically activated.             *
!*                                                           *
!* INPUT:                                                    *
!*   mol    - molecule (coords + optional lattice)           *
!*   chrg   - total molecular charge                         *
!*   pr     - optional verbosity flag (logical)              *
!*   iunit  - optional output unit                           *
!* OUTPUT:                                                   *
!*   ff_dat - initialized GFN-FF data object                 *
!*   io     - error status (0 = success)                     *
!*************************************************************
    implicit none
    type(coord),intent(in)      :: mol
    integer,intent(in)          :: chrg
    integer,intent(out)         :: io
    logical,intent(in),optional :: pr
    integer,intent(in),optional :: iunit
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
    type(coord) :: refmol
    !> LOCAL
    integer :: mylevel,myunit
    io = 0

    ! ── map legacy pr/iunit to integer printlevel/printunit ──────────────────
    mylevel = 0
    if (present(pr)) then
      if (pr) mylevel = 2
    end if
    if (present(iunit)) then
      myunit = iunit
    else
      myunit = stdout
    end if

#ifdef WITH_GFNFF
    if (allocated(ff_dat%refgeo)) then
      ! ── initialize from a separate reference structure ──────────────────────
      call refmol%open(ff_dat%refgeo)
      if (allocated(refmol%lat)) then
        call gfnff_initialize(refmol%nat,refmol%at,refmol%xyz,ff_dat, &
        &   ichrg=chrg,printlevel=mylevel,printunit=myunit,iostat=io, &
        &   lattice=refmol%lat,npbc=3)
      else
        call gfnff_initialize(refmol%nat,refmol%at,refmol%xyz,ff_dat, &
        &   ichrg=chrg,printlevel=mylevel,printunit=myunit,iostat=io)
      end if
      call refmol%deallocate()
    else
      ! ── initialize from mol directly ────────────────────────────────────────
      if (allocated(mol%lat)) then
        call gfnff_initialize(mol%nat,mol%at,mol%xyz,ff_dat, &
        &   ichrg=chrg,printlevel=mylevel,printunit=myunit,iostat=io, &
        &   lattice=mol%lat,npbc=3)
      else
        call gfnff_initialize(mol%nat,mol%at,mol%xyz,ff_dat, &
        &   ichrg=chrg,printlevel=mylevel,printunit=myunit,iostat=io)
      end if
    end if

#else /* WITH_GFNFF */
    write (stdout,*) 'Error: Compiled without GFN-FF support!'
    write (stdout,*) 'Use -DWITH_GFNFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfnff_api_setup

!========================================================================================!

  subroutine gfnff_sp(mol,ff_dat,energy,gradient,iostatus,sigma)
!*************************************************************
!* GFN-FF single-point energy + gradient for mol.           *
!* When mol%lat is allocated the periodic singlepoint call  *
!* is used automatically (lattice passed to gfnff).         *
!*                                                           *
!* INPUT:                                                    *
!*   mol      - molecule (coords + optional lattice)        *
!*   ff_dat   - initialized GFN-FF data object              *
!* OUTPUT:                                                   *
!*   energy   - total energy (Hartree)                      *
!*   gradient - gradient (Eh/Bohr)                          *
!*   iostatus - error status (0 = success)                  *
!*   sigma    - optional stress tensor (Eh); zero non-PBC   *
!*************************************************************
    implicit none
    !> INPUT
    type(coord),intent(in) :: mol
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
    !> OUTPUT
    real(wp),intent(out) :: energy
    real(wp),intent(out) :: gradient(3,mol%nat)
    integer,intent(out)  :: iostatus
    real(wp),intent(out),optional :: sigma(3,3)
    !> LOCAL
    real(wp) :: sigma_loc(3,3)
    energy = 0.0_wp
    gradient = 0.0_wp
    iostatus = 0
    sigma_loc = 0.0_wp
#ifdef WITH_GFNFF
    if (allocated(mol%lat)) then
      call gfnff_singlepoint(mol%nat,mol%at,mol%xyz,ff_dat, &
      &   energy,gradient,lattice=mol%lat,sigma=sigma_loc,iostat=iostatus)
    else
      call gfnff_singlepoint(mol%nat,mol%at,mol%xyz,ff_dat, &
      &   energy,gradient,iostat=iostatus)
    end if
    if (present(sigma)) sigma = sigma_loc
#else
    write (stdout,*) 'Error: Compiled without GFN-FF support!'
    write (stdout,*) 'Use -DWITH_GFNFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfnff_sp

!========================================================================================!

  subroutine gfnff_printout(iunit,ff_dat)
    implicit none
    !> INPUT
    integer,intent(in)  :: iunit
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
#ifdef WITH_GFNFF
    call ff_dat%resultprint(printunit=iunit)
#else
    write (stdout,*) 'Error: Compiled without GFN-FF support!'
    write (stdout,*) 'Use -DWITH_GFNFF=true in the setup to enable this function'
    error stop
#endif
  end subroutine gfnff_printout

!========================================================================================!
  subroutine gfnff_getwbos(ff_dat,nat,wbo)
!********************************************************
!* obtain connectivity information from GFN-FF topology
!* This is obviously not a true WBO
!********************************************************
    implicit none
    type(gfnff_data),intent(in) :: ff_dat
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)

    wbo = 0.0_wp
#ifdef WITH_GFNFF
    call gfnff_get_fake_wbo(ff_dat,nat,wbo)
#endif
  end subroutine gfnff_getwbos

!========================================================================================!

  subroutine gfnff_set_fragments(ff_dat,nat,fragments)
!***********************************************************************
!* Resolve a list of atom-range strings into a per-atom fragment index
!* and hand it to the GFN-FF data object (ff_dat%userinput%fraglist).
!*
!* Each entry of fragments(:) is an atom-range string (e.g. "1-432" or
!* "1,5-9,20") defining one fragment. Fragments are numbered 1..nfrag in
!* the given order. Any atom not covered by any string is left in
!* fragment 0 (the implicit "rest" group). GFN-FF will subsequently not
!* create bonds between atoms of different fragments.
!*
!* INPUT:
!*   nat       - number of atoms
!*   fragments - array of atom-range strings, one per fragment
!* IN/OUT:
!*   ff_dat    - GFN-FF data object; userinput%fraglist is (re)built
!***********************************************************************
    implicit none
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
    integer,intent(in) :: nat
    character(len=*),intent(in) :: fragments(:)
    !> LOCAL
    integer :: f,nsel
    integer,allocatable :: fraglist(:),mark(:)

    if (.not.allocated(ff_dat)) return

    allocate (fraglist(nat),source=0)
    allocate (mark(nat),source=0)
    do f = 1,size(fragments)
      call parse_atlist(fragments(f),nsel,nat,mark)
      where (mark == 1) fraglist = f
    end do
    deallocate (mark)

#ifdef WITH_GFNFF
    if (.not.allocated(ff_dat%userinput)) allocate (ff_dat%userinput)
    ff_dat%userinput%fraglist = fraglist
#endif
    deallocate (fraglist)
  end subroutine gfnff_set_fragments

!========================================================================================!

  subroutine gfnff_set_refcharges(ff_dat,q)
!***********************************************************************
!* Hand a set of atomic reference charges to the GFN-FF data object via
!* the in-memory bundle (ff_dat%userinput%refq). GFN-FF sums these over
!* each fragment to set the integer per-fragment net-charge constraint
!* used by its EEQ model. This is how the CEH guess assigns the correct
!* charge to each fragment (in-memory equivalent of the older, file-based
!* refcharges route).
!*
!* INPUT:
!*   q      - atomic reference charges (size = nat)
!* IN/OUT:
!*   ff_dat - GFN-FF data object; userinput%refq is (re)set
!***********************************************************************
    implicit none
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
    real(wp),intent(in) :: q(:)
    if (.not.allocated(ff_dat)) return
#ifdef WITH_GFNFF
    if (.not.allocated(ff_dat%userinput)) allocate (ff_dat%userinput)
    ff_dat%userinput%refq = q
#endif
  end subroutine gfnff_set_refcharges

!========================================================================================!

  subroutine gfnff_set_bondmat(ff_dat,nat,fname,io)
!***********************************************************************
!* Read a user-defined molecular graph from file and hand it to the
!* GFN-FF data object (ff_dat%userinput%bondmat). GFN-FF then skips its
!* distance-based bond perception and uses this connectivity instead.
!*
!* The file lists one bond per line as "i j [order]" with 1-based atom
!* indices; the bond order is an optional positive integer (default 1).
!* Blank lines and lines starting with '#' or '!' are skipped.
!* GFN-FF currently only uses the zero/nonzero pattern of the matrix.
!*
!* INPUT:
!*   nat    - number of atoms
!*   fname  - bond list file
!* IN/OUT:
!*   ff_dat - GFN-FF data object; userinput%bondmat is (re)built
!* OUTPUT:
!*   io     - error status (0 = success)
!***********************************************************************
    implicit none
    type(gfnff_data),allocatable,intent(inout) :: ff_dat
    integer,intent(in) :: nat
    character(len=*),intent(in) :: fname
    integer,intent(out) :: io
    !> LOCAL
    integer,allocatable :: bondmat(:,:)
    character(len=256) :: line
    integer :: ich,i,j,k,iline

    io = 0
    if (.not.allocated(ff_dat)) return

    open (newunit=ich,file=fname,status='old',action='read',iostat=io)
    if (io /= 0) then
      write (stderr,'(a,a)') 'could not open bond matrix file ',trim(fname)
      return
    end if

    allocate (bondmat(nat,nat),source=0)
    iline = 0
    do
      read (ich,'(a)',iostat=io) line
      if (io /= 0) exit
      iline = iline+1
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#'.or.line(1:1) == '!') cycle
      read (line,*,iostat=io) i,j,k
      if (io /= 0) then
        read (line,*,iostat=io) i,j
        k = 1
      end if
      if (io /= 0.or.i < 1.or.i > nat.or.j < 1.or.j > nat.or.i == j.or.k < 1) then
        write (stderr,'(a,i0,a,a,a)') 'invalid bond in line ',iline,' of ',trim(fname), &
        &                             ': '//trim(line)
        io = 1
        exit
      end if
      bondmat(i,j) = k
      bondmat(j,i) = k
    end do
    close (ich)
    if (io > 0) return
    io = 0

#ifdef WITH_GFNFF
    if (.not.allocated(ff_dat%userinput)) allocate (ff_dat%userinput)
    ff_dat%userinput%bondmat = bondmat
#endif
    deallocate (bondmat)
  end subroutine gfnff_set_bondmat

!========================================================================================!

  subroutine gfnff_dump_sasa(ff_dat,nat,atlist)
!***********************************************************
!* Dumps cumulative SASA for all the atoms .true. in atlist
!* into an file called fort.5454
!***********************************************************
    implicit none
    type(gfnff_data),intent(in) :: ff_dat
    logical,intent(in) :: atlist(nat)
    integer,intent(in) :: nat
    integer :: i
    real(wp) :: sumsasa
    if (allocated(ff_dat%solvation)) then
      if (allocated(ff_dat%solvation%sasa)) then
        sumsasa = 0.0_wp
        do i = 1,nat
          if (atlist(i)) sumsasa = sumsasa+ff_dat%solvation%sasa(i)
        end do
        write (5454,*) sumsasa
      end if
    end if
  end subroutine gfnff_dump_sasa
!========================================================================================!
!========================================================================================!
end module gfnff_api
