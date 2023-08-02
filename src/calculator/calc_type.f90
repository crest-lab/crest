!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2023 Philipp Pracht
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

module calc_type
  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use constraints
!>--- api types
  use tblite_api
  use gfn0_api
  use gfnff_api,only:gfnff_data
  use xhcff_api,only:xhcff_calculator
  implicit none

  character(len=1),public,parameter :: sep = '/'
  character(len=12),public,parameter :: dev0 = ' 2>/dev/null'

!&<
  !> job type enumerator
  type ,private:: enum_jobtype
    integer :: unknown   = 0
    integer :: xtbsys    = 1
    integer :: generic   = 2
    integer :: turbomole = 3
    integer :: orca      = 4
    integer :: terachem  = 5
    integer :: tblite    = 6
    integer :: gfn0      = 7
    integer :: gfn0occ   = 8
    integer :: gfnff     = 9
    integer :: xhcff     = 10
  end type enum_jobtype
  type(enum_jobtype), parameter,public :: jobtype = enum_jobtype()

  character(len=45),parameter,private :: jobdescription(11) = [ &
     & 'Unknown calculation type                    ', &
     & 'xTB calculation via external binary         ', &
     & 'Generic script execution                    ', &
     & 'Systemcall to the Turbomole program package ', &
     & 'Systemcall to the ORCA program package      ', &
     & 'Systemcall to the TeraChem program package  ', &
     & 'xTB calculation via tblite lib              ', &
     & 'GFN0-xTB calculation via GFN0 lib           ', &
     & 'GFN0*-xTB calculation via GFN0 lib          ', &
     & 'GFN-FF calculation via GFNFF lib            ', &
     & 'XHCFF calculation via XHCFF-lib             ' ]
!&>

!=========================================================================================!
!>--- data object that contains the data for a *SINGLE* calculation
  public :: calculation_settings
  type :: calculation_settings

    integer :: id  = 0  !> calculation type (see "jobtype" parameter above)
    integer :: prch = stdout

    integer :: chrg = 0
    integer :: uhf = 0
    real(wp) :: weight = 1.0_wp

    character(len=:),allocatable :: calcspace  !> subdirectory to perform the calculation in
    character(len=:),allocatable :: calcfile
    character(len=:),allocatable :: gradfile
    character(len=:),allocatable :: path
    character(len=:),allocatable :: other
    character(len=:),allocatable :: binary     !> binary or generic script
    character(len=:),allocatable :: systemcall !> systemcall for running generic scripts
    character(len=:),allocatable :: description

!>--- gradient format specifications
    integer :: gradtype = 0
    integer :: gradfmt = 0
    character(len=:),allocatable :: gradkey
    character(len=:),allocatable :: efile

!>--- results/property requests
    real(wp) :: epot = 0.0_wp
    real(wp) :: efix = 0.0_wp
    real(wp) :: etot = 0.0_wp

    !> bond orders
    logical :: rdwbo = .false.
    real(wp),allocatable :: wbo(:,:)

    !> dipole and dipole gradient
    logical :: rddip = .false.
    real(wp) :: dip(3) = 0
    logical :: rddipgrad = .false.
    real(wp),allocatable :: dipgrad(:,:,:)

!>--- API constructs
    integer  :: tblitelvl = 2
    real(wp) :: etemp = 300.0_wp
    real(wp) :: accuracy = 1.0_wp
    logical  :: apiclean = .true.
    integer  :: maxscc = 500
    logical  :: saveint = .false.
    character(len=:),allocatable :: solvmodel
    character(len=:),allocatable :: solvent

    !> tblite data
    type(wavefunction_type),allocatable  :: wfn
    type(tblite_calculator),allocatable  :: tbcalc
    type(tblite_ctx),allocatable         :: ctx
    type(tblite_resultstype),allocatable :: tbres
    type(wavefunction_type),allocatable  :: wfn_backup

    !> GFN0-xTB data
    type(gfn0_data),allocatable          :: g0calc
    integer :: nconfig = 0
    integer,allocatable :: config(:)
    real(wp),allocatable :: occ(:)

    !> GFN-FF data
    type(gfnff_data),allocatable :: ff_dat
    !> XHCFF data
    integer :: ngrid = 230             !>  lebedev grid points per atom
    real(wp) :: extpressure = 0.0_wp   !>  hydorstatic pressure in Gpa
    real(wp) :: proberad = 1.5_wp       !>  proberadius in a.u.
    integer :: vdwset = 0              !>  Set of VDW radii to use in sas calculation -> default D3, 1 -> Bondi
    type(xhcff_calculator),allocatable :: xhcff

!>--- Type procedures
  contains
    procedure :: deallocate => calculation_settings_deallocate
    procedure :: addconfig => calculation_settings_addconfig
    procedure :: autocomplete => calculation_settings_autocomplete
    procedure :: info => calculation_settings_info
  end type calculation_settings
!=========================================================================================!

!=========================================================================================!
!> data object that collects settings for *ALL* calculations and constraints.
  public :: calcdata
  type :: calcdata
    integer :: id = 0  !> this parameter will decide how to return or add up energies and gradients

!>--- calculations
    integer :: ncalculations = 0
    type(calculation_settings),allocatable :: calcs(:)
    real(wp),allocatable :: etmp(:)
    real(wp),allocatable :: grdtmp(:,:,:)
    real(wp),allocatable :: eweight(:)
    real(wp),allocatable :: etmp2(:)
    real(wp),allocatable :: grdtmp2(:,:,:)
    real(wp),allocatable :: eweight2(:)

!>--- constraints
    integer :: nconstraints = 0
    type(constraint),allocatable :: cons(:)

!>--- scans
    integer :: nscans = 0
    logical :: relaxscan = .true.
    real(wp) :: scansforce = 0.5_wp
    type(scantype),allocatable :: scans(:)

!>--- results/property requests
    real(wp) :: epot
    real(wp) :: efix
    real(wp) :: etot

!>--- optimization settings
    integer  :: optlev = 0
    integer  :: micro_opt = 20
    integer  :: maxcycle = 0
    real(wp) :: maxdispl_opt = 1.0_wp
    real(wp) :: hlow_opt = 0.010_wp
    real(wp) :: hmax_opt = 5.0_wp
    real(wp) :: acc_opt = 1.0_wp
    logical  :: exact_rf = .false.
    logical  :: average_conv = .false.
    logical  :: tsopt = .false.
    integer  :: iupdat = 0  !> 0=BFGS, 1=Powell, 2=SR1, 3=Bofill, 4=Schlegel

!>--- GFN0* data, needed for special MECP application
    type(gfn0_data),allocatable  :: g0calc

!>--- printouts and io
    logical :: pr_energies = .false.
    integer :: eout_unit = stdout
    character(len=:),allocatable :: elog

!>--- Type procedures
  contains
    procedure :: reset => calculation_reset
    procedure :: init => calculation_init
    generic,public :: add => calculation_add_constraint,calculation_add_settings,calculation_add_scan
    procedure,private :: calculation_add_constraint,calculation_add_settings,calculation_add_scan
    procedure :: copy => calculation_copy
    procedure :: printconstraints => calculation_print_constraints
    procedure :: removeconstraint => calculation_remove_constraint
    procedure :: info => calculation_info
  end type calcdata

!========================================================================================!
!========================================================================================!
contains  !>--- Module routines start here
!========================================================================================!
!========================================================================================!

  subroutine calculation_reset(self)
    implicit none
    class(calcdata) :: self

    self%id = 0

    if (allocated(self%calcs)) deallocate (self%calcs)
    self%ncalculations = 0

    if (allocated(self%cons)) deallocate (self%cons)
    self%nconstraints = 0

    self%optlev = 0
    self%micro_opt = 20
    self%maxcycle = 0
    self%maxdispl_opt = 1.0_wp
    self%hlow_opt = 0.010_wp
    self%hmax_opt = 5.0_wp
    self%acc_opt = 1.0_wp
    self%exact_rf = .false.
    self%average_conv = .false.
    self%tsopt = .false.
    self%iupdat = 0

    self%pr_energies = .false.
    self%eout_unit = stdout
    if (allocated(self%elog)) deallocate (self%elog)

    if (allocated(self%etmp)) deallocate (self%etmp)
    if (allocated(self%grdtmp)) deallocate (self%grdtmp)
    if (allocated(self%etmp2)) deallocate (self%etmp2)
    if (allocated(self%grdtmp2)) deallocate (self%grdtmp2)

    if (allocated(self%g0calc)) deallocate (self%g0calc)

    return
  end subroutine calculation_reset

  subroutine calculation_settings_deallocate(self)
    implicit none
    class(calculation_settings) :: self

    if (allocated(self%calcspace)) deallocate (self%calcspace)
    if (allocated(self%calcfile)) deallocate (self%calcfile)
    if (allocated(self%gradfile)) deallocate (self%gradfile)
    if (allocated(self%path)) deallocate (self%path)
    if (allocated(self%other)) deallocate (self%other)
    if (allocated(self%binary)) deallocate (self%binary)
    if (allocated(self%systemcall)) deallocate (self%systemcall)
    if (allocated(self%description)) deallocate (self%description)
    if (allocated(self%wbo)) deallocate (self%wbo)
    if (allocated(self%dipgrad)) deallocate (self%dipgrad)
    if (allocated(self%gradkey)) deallocate (self%gradkey)
    if (allocated(self%efile)) deallocate (self%efile)
    if (allocated(self%solvmodel)) deallocate (self%solvmodel)
    if (allocated(self%solvent)) deallocate (self%solvent)
    if (allocated(self%wfn)) deallocate (self%wfn)
    if (allocated(self%tbcalc)) deallocate (self%tbcalc)
    if (allocated(self%ctx)) deallocate (self%ctx)
    if (allocated(self%tbres)) deallocate (self%tbres)
    if (allocated(self%wfn_backup)) deallocate (self%wfn_backup)
    if (allocated(self%g0calc)) deallocate (self%g0calc)
    if (allocated(self%ff_dat)) deallocate (self%ff_dat)
    if (allocated(self%xhcff)) deallocate(self%xhcff)

    self%id = 0
    self%prch = stdout
    self%chrg = 0
    self%uhf = 0
    self%epot = 0.0_wp
    self%efix = 0.0_wp
    self%etot = 0.0_wp

    self%rdwbo = .false.
    self%rddip = .false.
    self%dip = 0.0_wp
    self%rddipgrad = .false.
    self%gradtype = 0
    self%gradfmt = 0

    self%tblitelvl = 2
    self%etemp = 300.0_wp
    self%accuracy = 1.0_wp
    self%apiclean = .false.
    self%maxscc = 500
    self%saveint = .false.

    self%ngrid = 230 
    self%extpressure = 0.0_wp
    self%proberad = 1.5_wp


    return
  end subroutine calculation_settings_deallocate

!=========================================================================================!

  subroutine calculation_add_settings(self,cal)
    implicit none
    class(calcdata) :: self
    type(calculation_settings) :: cal
    type(calculation_settings),allocatable :: callist(:)
    integer :: i,j

    if (self%ncalculations < 1) then
      allocate (self%calcs(1))
      self%ncalculations = 1
      self%calcs(1) = cal
    else
      i = self%ncalculations+1
      j = self%ncalculations
      allocate (callist(i))
      callist(1:j) = self%calcs(1:j)
      callist(i) = cal
      call move_alloc(callist,self%calcs)
      self%ncalculations = i
    end if

    return
  end subroutine calculation_add_settings

!=========================================================================================!

  subroutine calculation_add_constraint(self,constr)
    implicit none
    class(calcdata) :: self
    type(constraint) :: constr
    type(constraint),allocatable :: conslist(:)
    integer :: i,j

    if (self%nconstraints < 1) then
      allocate (self%cons(1))
      self%nconstraints = 1
      self%cons(1) = constr
    else
      i = self%nconstraints+1
      j = self%nconstraints
      allocate (conslist(i))
      conslist(1:j) = self%cons(1:j)
      conslist(i) = constr
      call move_alloc(conslist,self%cons)
      self%nconstraints = i
    end if

    return
  end subroutine calculation_add_constraint

!=========================================================================================!

  subroutine calculation_add_scan(self,scn)
    implicit none
    class(calcdata) :: self
    type(scantype) :: scn
    type(scantype),allocatable :: scnlist(:)
    integer :: i,j

    if (self%nscans < 1) then
      allocate (self%scans(1))
      self%nscans = 1
      self%scans(1) = scn
    else
      i = self%nscans+1
      j = self%nscans
      allocate (scnlist(i))
      scnlist(1:j) = self%scans(1:j)
      scnlist(i) = scn
      call move_alloc(scnlist,self%scans)
      self%nscans = i
    end if

    return
  end subroutine calculation_add_scan

!=========================================================================================!

  subroutine calculation_remove_constraint(self,d)
    implicit none
    class(calcdata) :: self
    type(constraint) :: constr
    type(constraint),allocatable :: conslist(:)
    integer :: i,j,d,d1,d2

    if (self%nconstraints < d) return

    i = self%nconstraints-1
    j = self%nconstraints
    allocate (conslist(i))
    if (d == 1) then
      conslist(1:i) = self%cons(2:j)
    else if (d == j) then
      conslist(1:i) = self%cons(1:i)
    else
      d1 = d-1
      d2 = d+1
      conslist(1:d1) = self%cons(1:d1)
      conslist(d:i) = self%cons(d2:j)
    end if
    call move_alloc(conslist,self%cons)
    self%nconstraints = i

    return
  end subroutine calculation_remove_constraint

!=========================================================================================!
  subroutine calculation_print_constraints(self,chnl)
    implicit none
    class(calcdata) :: self
    integer :: i,j
    integer,optional :: chnl

    if (self%nconstraints < 1) then
      return
    else
      if (present(chnl)) then
        i = chnl
      else
        i = stdout
      end if
      !write(i,*)
      do j = 1,self%nconstraints
        call self%cons(j)%print(i)
      end do
    end if

    return
  end subroutine calculation_print_constraints

!=========================================================================================!

  subroutine calculation_init(self)
    class(calcdata) :: self

    if (allocated(self%elog)) then
      self%pr_energies = .true.
      open (newunit=self%eout_unit,file=self%elog)
    end if

  end subroutine calculation_init

!=========================================================================================!
!> copy a calcdata object from src to self
  subroutine calculation_copy(self,src)
    class(calcdata) :: self
    type(calcdata) :: src
    integer :: i
    !call self%reset()

    self%id = src%id

    self%ncalculations = src%ncalculations
    if (allocated(self%calcs)) deallocate (self%calcs)
    !self%calcs = src%calcs
    do i = 1,self%ncalculations
      call self%add(src%calcs(i))
    end do

    self%nconstraints = src%nconstraints
    if (allocated(self%cons)) deallocate (self%cons)
    !self%cons = src%cons
    do i = 1,self%nconstraints
      call self%add(src%cons(i))
    end do

    self%optlev = src%optlev
    self%micro_opt = src%micro_opt
    self%maxcycle = src%maxcycle
    self%maxdispl_opt = src%maxdispl_opt
    self%hlow_opt = src%hlow_opt
    self%hmax_opt = src%hmax_opt
    self%acc_opt = src%acc_opt
    self%exact_rf = src%exact_rf
    self%average_conv = src%average_conv
    self%tsopt = src%tsopt
    self%iupdat = src%iupdat

    self%pr_energies = src%pr_energies
    self%eout_unit = src%eout_unit
    self%elog = src%elog

    return
  end subroutine calculation_copy

!=========================================================================================!

  subroutine calculation_settings_addconfig(self,config)
    implicit none
    class(calculation_settings) :: self
    integer,intent(in)  :: config(:)
    integer :: i,j
    integer :: l,lold,lnew,n,nnew
    integer,allocatable :: configtmp(:,:)

    l = size(config,1)
    if (allocated(self%config)) deallocate (self%config)
    allocate (self%config(l))
    self%config = config
  end subroutine calculation_settings_addconfig

!=========================================================================================!

!>-- check for missing settings in a calculation_settings object
  subroutine calculation_settings_autocomplete(self,id)
    implicit none
    class(calculation_settings) :: self
    integer,intent(in)  :: id
    integer :: i,j
    character(len=50) :: nmbr

    !> add a short description
    self%description = trim(jobdescription(self%id+1))

    if (.not.allocated(self%calcspace)) then
      !> I've decided to perform all calculations in a separate directory to
      !> avoid accumulation of files in the main workspace
      write (nmbr,'(i0)') id
      self%calcspace = 'calculation.level.'//trim(nmbr)
    end if
  end subroutine calculation_settings_autocomplete

!=========================================================================================!
  subroutine calculation_settings_info(self,iunit)
    implicit none
    class(calculation_settings) :: self
    integer,intent(in) :: iunit
    integer :: i,j
    character(len=*),parameter :: fmt1 = '(1x,a20," : ",i5)'
    character(len=*),parameter :: fmt2 = '(1x,a20," : ",f12.5)'
    character(len=*),parameter :: fmt3 = '(1x,a20," : ",a)'

    if (allocated(self%description)) then
      write (iunit,'(2x,a)') trim(self%description)
    else
      write (iunit,fmt1) 'Job type',self%id
    end if

    write (iunit,fmt1) 'Molecular charge',self%chrg
    if (self%uhf /= 0) then
      write (iunit,fmt1) 'UHF parameter',self%uhf
    end if

    if (allocated(self%solvmodel)) then
      write (iunit,fmt3) 'Solvation model',trim(self%solvmodel)
    end if
    if (allocated(self%solvent)) then
      write (iunit,fmt3) 'Solvent',trim(self%solvent)
    end if

    !> xTB specific parameters
    if (any((/jobtype%tblite,jobtype%xtbsys,jobtype%gfn0,jobtype%gfn0occ/) == self%id)) then
      write (iunit,fmt2) 'Fermi temperature',self%etemp
      write (iunit,fmt2) 'Accuracy',self%accuracy
      write (iunit,fmt1) 'max SCC cycles',self%maxscc
    end if

    if (self%apiclean) write (iunit,fmt3) 'Reset data?','  yes'
    if (self%rdwbo) write (iunit,fmt3) 'Read WBOs?','  yes'
    if (self%rddip) write (iunit,fmt3) 'Read dipoles?','  yes'

  end subroutine calculation_settings_info
!========================================================================================!
  subroutine calculation_info(self,iunit)
    implicit none
    class(calcdata) :: self
    integer,intent(in) :: iunit
    integer :: i
    character(len=*),parameter :: fmt1 = '(1x,a20," : ",i5)'
    character(len=*),parameter :: fmt2 = '(1x,a20," : ",f12.5)'

    write (iunit,'(1x,a)') '----------------'
    write (iunit,'(1x,a)') 'Calculation info'
    write (iunit,'(1x,a)') '----------------'
    if (self%ncalculations <= 0) then
      write (iunit,'("> ",a)') 'No calculation levels set up!'
    else if (self%ncalculations > 1) then
      do i = 1,self%ncalculations
        write (iunit,'("> ",a,i0)') 'Calculation level ',i
        call self%calcs(i)%info(iunit)
        write (iunit,fmt2) 'weight',self%calcs(i)%weight
      end do
    else
      write (iunit,'("> ",a)') 'Calculation level'
      call self%calcs(1)%info(iunit)
    end if
    write (iunit,*)

    if (self%nconstraints > 0) then
      write (iunit,'("> ",a)') 'User-defined constraints:'
      do i = 1,self%nconstraints
        call self%cons(i)%print(iunit)
      end do
      write (iunit,*)
    end if

    if (self%ncalculations > 1) then
      write (iunit,'("> ",a)') 'Energy and gradient processing:'
      select case (self%id)
      case (1:)
        write (iunit,'(1x,a,i0)') 'Energies and gradients of calculation level ',self%id, &
        &  ' will be used'
      case (-1)
        write (iunit,'(1x,a)') 'Special MECP energy and gradients will be constructed'
        write (iunit,'(1x,a)') 'See https://doi.org/10.1021/acs.jpclett.3c00494'
      case default
        write (iunit,'(1x,a)') 'Energies and gradients of all calculation levels will be'// &
        & ' added according to their weights'
      end select
      write (iunit,*)
    end if

    return
  end subroutine calculation_info

!=========================================================================================!
end module calc_type
