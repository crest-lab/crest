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
  use crest_external_engrad,only:engrad_interface
  use constraints
  use strucrd,only:coord
!>--- api types
  use tblite_api
  use gfn0_api
  use crest_solvation,only:solvation_data
  use crest_electrostatic,only:electrostatic_data
  use gfnff_api,only:gfnff_data
  use libpvol_api,only:libpvol_calculator
!>--- other types
  use orca_type
  use lwoniom_module
  use hessian_reconstruct
  use approxg_module, only: approxg_params
  use penalty_module, only: penalty_params
  use mlip_sc, only: mlip_params
  implicit none

  character(len=*),public,parameter :: sep = '/'
  character(len=*),public,parameter :: dev0 = ' 2>/dev/null'

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
    integer :: libpvol   = 10
    integer :: lj        = 11
    integer :: approxg   = 12
    integer :: penalty   = 13
    integer :: mlip      = 14
    integer :: solvation = 15
    integer :: electrostatic = 16
    integer :: external  = 17
  end type enum_jobtype
  type(enum_jobtype), parameter,public :: jobtype = enum_jobtype()

  character(len=45),parameter,private :: jobdescription(18) = [character(len=45) :: &
     & 'Unknown calculation type                    ', &
     & 'xTB calculation via external binary         ', &
     & 'Generic script execution                    ', &
     & 'Systemcall with Turbomole-style in/output   ', &
     & 'Systemcall to the ORCA program package      ', &
     & 'Systemcall to the TeraChem program package  ', &
     & 'xTB calculation via tblite lib              ', &
     & 'GFN0-xTB calculation via GFN0 lib           ', &
     & 'GFN0*-xTB calculation via GFN0 lib          ', &
     & 'GFN-FF calculation via GFNFF lib            ', &
     & 'external pressure calculation via libpvol   ', &
     & 'Lennard-Jones potential calculation         ', &
     & 'Approximate free energy computation         ', &
     & 'Empirical penalty function                  ', &
     & 'MLIP via persistent python socket           ', &
     & 'Standalone implicit solvation contribution  ', &
     & 'Charge-equilibration electrostatics         ', &
     & 'Externally supplied potential (host callback)']
!&>

!=========================================================================================!
!>--- RE-EXPORT of the externally-supplied-potential interface
!>    (defined in module crest_external_engrad)
  public :: engrad_interface

!=========================================================================================!

  public :: calculation_settings
  type :: calculation_settings
!**********************************************************************
!* data object that contains the data for a *SINGLE* calculation level
!**********************************************************************
    integer :: id = 0         !> calculation type (see "jobtype" parameter above)
    integer :: prch = stdout  !> printout channel
    logical :: pr = .false.   !> allow the calculation to produce printout? Results in a lot I/O
    logical :: prappend = .false. !> append printout
    logical :: prstdout = .false. !> special case, fwd some printout to stdout
    integer :: refine_lvl = 0 !> to allow defining different refinement levels

    integer :: chrg = 0          !> molecular charge
    integer :: uhf = 0           !> uhf = Nα-Nβ = 2S (xtb/tblite/gfn0/gfnff/turbomole)
    integer :: multiplicity = 1  !> spin multiplicity = 2S+1 = uhf+1 (ORCA, fmlip-relay)
                                 !> kept aligned with uhf via sync_multiplicity/set_multiplicity
    logical :: active = .true.   !> active setting to disable the calculation (this is different from weight=0)
    real(wp) :: weight = 1.0_wp  !> calculation weight (when adding them up)
    integer :: threads = 0       !> max number of cores this level may use per engrad call
                                 !> (e.g. ORCA %pal nprocs, MLIP server threads).
                                 !> 0 (or <1) means "unset/inherit" and is ignored.

    character(len=:),allocatable :: calcspace  !> subdirectory to perform the calculation in
    character(len=:),allocatable :: calcfile
    character(len=:),allocatable :: gradfile
    character(len=:),allocatable :: path
    character(len=:),allocatable :: other
    character(len=:),allocatable :: binary      !> binary or generic script
    character(len=:),allocatable :: systemcall  !> systemcall for running generic scripts
    character(len=:),allocatable :: description !> see above jobdescription parameter
    character(len=:),allocatable :: shortflag   !> shorter job description

!>--- gradient format specifications
    logical :: numgrad = .false.      !> run numerical gradient (expensive!)
    real(wp) :: gradstep = 0.0005_wp  !> displacement for numerical gradient
    logical :: rdgrad = .true.
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

    !> atomic charges
    logical :: rdqat = .false.
    logical :: dumpq = .false.
    real(wp),allocatable :: qat(:)

    !> dipole and dipole gradient
    logical :: rddip = .false.
    real(wp) :: dipole(3) = 0.0_wp
    logical :: rddipgrad = .false.
    real(wp),allocatable :: dipgrad(:,:)

    !> other properties
    logical,allocatable :: getsasa(:)
    logical :: getlmocent = .false.
    integer :: nprot = 0
    real(wp),allocatable :: protxyz(:,:)
    real(wp),allocatable :: efield(:)  !> in V/Å

!>--- API constructs
    integer  :: tblitelvl = 2
    real(wp) :: etemp = 300.0_wp
    logical  :: etemp_user_set = .false.
    real(wp) :: accuracy = 1.0_wp
    logical  :: apiclean = .true.
    integer  :: maxscc = 500
    logical  :: saveint = .false.
    character(len=:),allocatable :: solvmodel
    character(len=:),allocatable :: solvent
!>--- Some optional file name storages
    character(len=:),allocatable :: parametrisation
    logical  :: restart = .false.  !> restart option (some potentials can do this)
    character(len=:),allocatable :: restartfile
    character(len=:),allocatable :: refgeo
    character(len=:),allocatable :: refcharges

!>--- tblite data
    type(tblite_data),allocatable :: tblite
    character(len=:),allocatable :: tbliteparam
    logical :: ceh_guess = .false.
    logical :: spin_polarized = .false.

!>--- GFN0-xTB data
    type(gfn0_data),allocatable :: g0calc
    integer :: nconfig = 0
    integer,allocatable :: config(:)
    real(wp),allocatable :: occ(:)

!>--- GFN-FF data
    type(gfnff_data),allocatable :: ff_dat
    !> optional user-defined fragmentation (one atom-range string per fragment,
    !> e.g. ["1-432","433-458"]). Atoms not covered by any string form one
    !> additional "rest" fragment. Passed to GFN-FF to prevent artificial bonds
    !> between fragments. Resolved to a per-atom list once nat is known.
    character(len=:),allocatable :: gff_fragments(:)

!>--- libpvol data
    integer  :: pvmodel = 1            !> libpvol model type (0=XHCFF, 1=PV)
    integer  :: ngrid = 1202           !> lebedev grid points per atom
    real(wp) :: extpressure = 0.0_wp   !> hydorstatic pressure in Gpa
    real(wp) :: proberad = 1.5_wp      !> proberadius in Angstroem
    integer  :: vdwset = 0             !> Type of VDW radii -> 0 (default) D3, 1 -> Bondi
    real(wp) :: pvradscal = 1.0_wp     !> Scaling factor for SAS radii
    type(libpvol_calculator),allocatable :: libpvol

!>--- implicit solvation composite data
    type(solvation_data),allocatable :: solv

!>--- charge-equilibration electrostatics data
    type(electrostatic_data),allocatable :: eeq

!>--- approxg data
    type(approxg_params) :: ag

!>--- penalty params
    type(penalty_params) :: penalty

    !> ONIOM fragment IDs
    integer :: ONIOM_highlowroot = 0
    integer :: ONIOM_id = 0

!>--- ORCA job template
    type(orca_input) :: ORCA

!>--- MLIP settings
    type(mlip_params) :: MPAR

!>--- externally supplied potential (host-program callback)
!>    native Fortran procedure pointer + unlimited polymorphic context
    procedure(engrad_interface),pointer,nopass :: ext_engrad => null()
    class(*),pointer :: ext_userdata => null()  !> opaque host context

!>--- Type procedures
  contains
    procedure :: deallocate => calculation_settings_deallocate
    procedure :: set_external => calculation_settings_set_external
    procedure :: addconfig => calculation_settings_addconfig
    procedure :: autocomplete => calculation_settings_autocomplete
    procedure :: printid => calculation_settings_printid
    procedure :: info => calculation_settings_info
    procedure :: create => create_calclevel_shortcut
    procedure :: norestarts => calculation_settings_norestarts
    procedure :: dumpdipgrad => calculation_dump_dipgrad
    procedure :: copy => calculation_settings_copy
    procedure :: sync_multiplicity => calc_sync_multiplicity
    procedure :: set_multiplicity => calc_set_multiplicity
  end type calculation_settings

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!

  public :: calcdata
  type :: calcdata
!*****************************************************************************
!* data object that collects settings for *ALL* calculations and constraints.
!*****************************************************************************
    integer :: id = 0  !> this parameter will decide how to return or add up energies and gradients
    integer :: refine_stage = 0 !> to allow iterating different refinement stages

!>--- calculations
    integer :: ncalculations = 0
    type(calculation_settings),allocatable :: calcs(:)
    real(wp),allocatable :: etmp(:)
    real(wp),allocatable :: grdtmp(:,:,:)
    real(wp),allocatable :: eweight(:)
    real(wp),allocatable :: weightbackup(:)
    logical,allocatable  :: activebackup(:)
    real(wp),allocatable :: etmp2(:)
    real(wp),allocatable :: grdtmp2(:,:,:)
    real(wp),allocatable :: eweight2(:)
    real(wp),allocatable :: grdfix(:,:)

!>--- constraints
    integer :: nconstraints = 0
    type(constraint),allocatable :: cons(:)

!>--- scans
    integer :: nscans = 0
    logical :: relaxscan = .true.
    real(wp) :: scansforce = 0.5_wp
    type(scantype),allocatable :: scans(:)

!>--- frozen atoms
    integer :: nfreeze = 0
    logical,allocatable :: freezelist(:)

!>--- results/property requests
    real(wp) :: epot
    real(wp) :: efix
    real(wp) :: etot

!>--- optimization settings
    logical  :: optnewinit = .false.  !> ensure fresh calc/param setup at beginning of opt
    logical  :: anopt = .false.  !> allow collecting structures that ran into maxcycle?
    integer  :: optlev = 0
    integer  :: micro_opt = 20
    integer  :: maxcycle = 0
    real(wp) :: maxdispl_opt = 1.0_wp
    real(wp) :: ethr_opt = -1.0_wp  !> ΔE convergence
    real(wp) :: gthr_opt = -1.0_wp  !> RMS force convergence
    real(wp) :: hlow_opt = 0.010_wp
    real(wp) :: hmax_opt = 5.0_wp
    real(wp) :: acc_opt = 1.0_wp
    real(wp) :: maxerise = 1.0e-5_wp
    real(wp) :: hguess = 0.02_wp !> Hessian guess (for some optimizers)
    logical  :: exact_rf = .false.
    logical  :: average_conv = .false.
    logical  :: tsopt = .false.
    integer  :: iupdat = 0  !> 0=BFGS, 1=Powell, 2=SR1, 3=Bofill, 4=Schlegel
    integer  :: opt_engine = 0 !> default: ANCOPT
    integer  :: lbfgs_histsize = 20  !> L-BFGS history size
    integer  :: hess_init = 5 !> Initialization of the hessian, standard modhess lindh95
    logical  :: logextxyz = .true.  !> write extended xyz files from optimization trajectories

!>--- GFN0* data, needed for special MECP application
    type(gfn0_data),allocatable  :: g0calc

!>--- printouts and IO
    logical :: pr_energies = .false.
    integer :: eout_unit = stdout
    character(len=:),allocatable :: elog

!>--- ONIOM calculator data
    type(lwoniom_data),allocatable :: ONIOM
    type(coord),allocatable :: ONIOMmols(:)
    integer,allocatable :: ONIOMmap(:) !> map lwONIOM jobs to calculation_settings

!>--- Hessian Reconstructor and Thermo data
    type(cashed_hessian),allocatable :: chess
    logical :: do_HR = .false.
    logical :: full_HR = .false. !> Keyword for HR with all opt steps
    integer :: hu_steps = 10 !> default number of update steps
    integer :: nt !> following all required for thermochemistry
    real(wp),allocatable :: temperatures(:)
    real(wp),allocatable :: et(:),ht(:),gt(:),stot(:)
    real(wp) :: ithr,fscal,sthr
    character(len=:),allocatable :: emodel
    integer :: initialize_hr_type !> case defining initialization
    integer :: mh_type = 0
    integer :: hr_hu_type = 0
    logical :: deform_opt_hess = .false.
    real(wp) :: doh_stepsize = 0.10_wp !>stepsize for the deformation/reoptimization hessian generation
    real(wp) :: chess_id_guess = 0.1_wp
    logical :: g_sampling = .false. !>Do sampling on free energy surface as approximated using the lindh95 hessian
    integer :: gs_hess_type = 5

!>--- Parameters for smooth function within optimizer
    real(wp) :: L = 1.50_wp
    real(wp) :: k = 5000.0_wp
    real(wp) :: shift = 0.0006_wp
    real(wp) :: scaling = 0.1_wp

!>--- Type procedures
  contains
    procedure :: reset => calculation_reset
    procedure :: init => calculation_init
    procedure :: create => calcdata_create_shortcut
    generic,public :: add => calculation_add_constraint,calculation_add_settings, &
    & calculation_add_scan,calculation_add_constraintlist
    procedure,private :: calculation_add_constraint,calculation_add_settings, &
    & calculation_add_scan,calculation_add_constraintlist
    procedure :: copy => calculation_copy
    procedure :: printconstraints => calculation_print_constraints
    procedure :: removeconstraint => calculation_remove_constraint
    procedure :: info => calculation_info
    procedure :: ONIOMexpand => calculation_ONIOMexpand
    procedure :: active => calc_set_active
    procedure :: active_restore => calc_set_active_restore
    generic,public :: set_freeze => calculation_set_freeze_range,calculation_set_freeze_bools
    procedure,private :: calculation_set_freeze_range,calculation_set_freeze_bools
    procedure :: freezegrad => calculation_freezegrad
    procedure :: set_charge => calculation_set_charge
    procedure :: increase_charge => calculation_increase_charge
    procedure :: decrease_charge => calculation_decrease_charge
    procedure :: dealloc_params => calculation_deallocate_params
    procedure :: maxthreads => calculation_maxthreads
    procedure :: maxthreads_capped => calculation_maxthreads_capped
  end type calcdata

  public :: get_dipoles

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
    self%ethr_opt = -1.0_wp  !> ΔE convergence
    self%gthr_opt = -1.0_wp  !> RMS force convergence
    self%hlow_opt = 0.010_wp
    self%hmax_opt = 5.0_wp
    self%acc_opt = 1.0_wp
    self%maxerise = 1.0e-5_wp
    self%exact_rf = .false.
    self%average_conv = .false.
    self%tsopt = .false.
    self%iupdat = 0  !> 0=BFGS, 1=Powell, 2=SR1, 3=Bofill, 4=Schlegel
    self%opt_engine = 0 !> default: ANCOPT

    self%pr_energies = .false.
    self%eout_unit = stdout
    if (allocated(self%elog)) deallocate (self%elog)

    if (allocated(self%etmp)) deallocate (self%etmp)
    if (allocated(self%grdtmp)) deallocate (self%grdtmp)
    if (allocated(self%etmp2)) deallocate (self%etmp2)
    if (allocated(self%grdtmp2)) deallocate (self%grdtmp2)

    if (allocated(self%g0calc)) deallocate (self%g0calc)

    if (allocated(self%ONIOM)) deallocate (self%ONIOM)
    if (allocated(self%ONIOMmols)) deallocate (self%ONIOMmols)
    if (allocated(self%ONIOMmap)) deallocate (self%ONIOMmap)

    return
  end subroutine calculation_reset

!=========================================================================================!

  subroutine calculation_deallocate_params(self)
!**********************************************
!* Deallocate parametrization/model setup
!* This only applies (or is useful) for a few
!* selected model like GFN-FF, gfn0 and tblite
!**********************************************
    class(calcdata) :: self
    integer :: i,j,k
    if (self%ncalculations > 0) then
      do i = 1,self%ncalculations
        if (allocated(self%calcs(i)%tblite)) deallocate (self%calcs(i)%tblite)
        if (allocated(self%calcs(i)%g0calc)) deallocate (self%calcs(i)%g0calc)
        if (allocated(self%calcs(i)%ff_dat)) deallocate (self%calcs(i)%ff_dat)
        if (allocated(self%calcs(i)%libpvol)) deallocate (self%calcs(i)%libpvol)
      end do
    end if
  end subroutine calculation_deallocate_params

!=========================================================================================!

  function calculation_maxthreads(self) result(maxthreads)
!************************************************************
!* Return the largest per-level thread count requested by
!* any active calculation level (default 1). Intended to
!* inform the outer OMP thread-budgeting (cores per job).
!************************************************************
    class(calcdata),intent(in) :: self
    integer :: maxthreads
    integer :: i
    maxthreads = 1
    if (self%ncalculations > 0) then
      do i = 1,self%ncalculations
        if (.not.self%calcs(i)%active) cycle
        maxthreads = max(maxthreads,self%calcs(i)%threads)
      end do
    end if
  end function calculation_maxthreads

!=========================================================================================!

  function calculation_maxthreads_capped(self) result(maxcap)
!************************************************************
!* Largest thread reservation among active levels whose
!* backend is HARD-CAPPED at that count, i.e. uses exactly
!* threads cores and does not grow into cores_per_job via the
!* OMP_NUM_THREADS env var / nested OpenMP. Hard-capped are the
!* ORCA subprocess (its %pal is written from threads) and the
!* fmlip-relay MLIP servers (spawned with --max-threads from
!* threads); internal API calculators and generic subprocesses
!* pick up cores_per_job instead and thus soak leftover cores.
!*
!* Returns 0 if no capped level is active. Used to report
!* genuinely idle cores (see new_ompautoset).
!************************************************************
    class(calcdata),intent(in) :: self
    integer :: maxcap
    integer :: i
    maxcap = 0
    if (self%ncalculations > 0) then
      do i = 1,self%ncalculations
        if (.not.self%calcs(i)%active) cycle
        if (self%calcs(i)%id == jobtype%orca .or. &
        &   self%calcs(i)%id == jobtype%mlip) then
          maxcap = max(maxcap,self%calcs(i)%threads)
        end if
      end do
    end if
  end function calculation_maxthreads_capped

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

  subroutine calcdata_create_shortcut(self,levelstring, &
      & chrg,uhf,solvmodel,solvent)
!*********************************************************************
!* subroutine calcdata_create_shortcut called with %create(...)
!* Quick clean setup of a *fresh* calcdata object holding exactly
!* one calculation level. The calcdata is reset first, then a single
!* calculation_settings object is built (via its own %create) for the
!* requested level of theory and registered in self%calcs(1).
!*
!* This is the calcdata-level counterpart to the %create shortcut of
!* the calculation_settings type and exists as an internal code
!* shortcut only (it does not sanity-check the chosen settings).
!*
!* Required argument:
!*  - levelstring : level of theory, passed on to settings%create
!* Optional arguments (forwarded to settings%create):
!*  - chrg        : molecular charge  (integer)
!*  - uhf         : uhf parameter     (integer)
!*  - solvmodel   : solvation model   (only together with solvent)
!*  - solvent     : implicit solvent  (only together with solvmodel)
!*********************************************************************
    implicit none
    class(calcdata),intent(inout) :: self
    character(len=*),intent(in) :: levelstring
    integer,intent(in),optional :: chrg
    integer,intent(in),optional :: uhf
    character(len=*),intent(in),optional :: solvmodel
    character(len=*),intent(in),optional :: solvent
    type(calculation_settings) :: job

    ! ── start from a clean calcdata object ──────────────────────────
    call self%reset()

    ! ── build the single calculation level (absent optionals are ────
    ! ── forwarded as absent to settings%create) ─────────────────────
    call job%create(levelstring,chrg=chrg,uhf=uhf, &
    &               solvmodel=solvmodel,solvent=solvent)

    ! ── register it as the only level in this calcdata ──────────────
    call self%add(job)

    return
  end subroutine calcdata_create_shortcut

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

  subroutine calculation_add_constraintlist(self,k,constr)
    implicit none
    class(calcdata) :: self
    integer :: k
    type(constraint) :: constr(k)
    type(constraint),allocatable :: conslist(:)
    integer :: i,j

    if (self%nconstraints < 1) then
      allocate (self%cons(k))
      self%nconstraints = k
      self%cons(1:k) = constr(1:k)
    else
      i = self%nconstraints+k
      j = self%nconstraints
      allocate (conslist(i))
      conslist(1:j) = self%cons(1:j)
      conslist(j+1:i) = constr(1:k)
      call move_alloc(conslist,self%cons)
      self%nconstraints = i
    end if

    return
  end subroutine calculation_add_constraintlist

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
  subroutine calculation_copy(self,src,ignore_constraints)
    class(calcdata) :: self
    type(calcdata),intent(in) :: src
    logical,intent(in),optional :: ignore_constraints
    type(calculation_settings) :: newset
    type(constraint) :: newcons
    integer :: i
    logical :: igno

    call self%reset()

! ── identity ─────────────────────────────────────────────────────────────────
    self%id           = src%id
    self%refine_stage = src%refine_stage

! ── calculation levels ───────────────────────────────────────────────────────
    if (allocated(self%calcs)) deallocate(self%calcs)
    self%ncalculations = 0
    do i = 1,src%ncalculations
      call newset%copy(src%calcs(i))
      call self%add(newset)
    end do

! ── constraints ──────────────────────────────────────────────────────────────
    igno = .false.
    if (present(ignore_constraints)) igno = ignore_constraints
    if (allocated(self%cons)) deallocate(self%cons)
    self%nconstraints = 0
    if (.not.igno) then
      do i = 1,src%nconstraints
        call newcons%copy(src%cons(i))
        call self%add(newcons)
      end do
    end if

! ── scans ────────────────────────────────────────────────────────────────────
    self%nscans    = src%nscans
    self%relaxscan = src%relaxscan
    self%scansforce = src%scansforce
    if (allocated(src%scans)) then
      allocate(self%scans(src%nscans))
      do i = 1,src%nscans
        self%scans(i)%type        = src%scans(i)%type
        self%scans(i)%n           = src%scans(i)%n
        self%scans(i)%steps       = src%scans(i)%steps
        self%scans(i)%minval      = src%scans(i)%minval
        self%scans(i)%maxval      = src%scans(i)%maxval
        self%scans(i)%constrnmbr  = src%scans(i)%constrnmbr
        self%scans(i)%restore     = src%scans(i)%restore
        self%scans(i)%currentstep = src%scans(i)%currentstep
        if (allocated(src%scans(i)%atms))   self%scans(i)%atms   = src%scans(i)%atms
        if (allocated(src%scans(i)%points)) self%scans(i)%points = src%scans(i)%points
      end do
    end if

! ── frozen atoms ─────────────────────────────────────────────────────────────
    self%nfreeze = src%nfreeze
    if (allocated(src%freezelist)) self%freezelist = src%freezelist

! ── optimization settings ────────────────────────────────────────────────────
    self%optnewinit     = src%optnewinit
    self%anopt          = src%anopt
    self%optlev         = src%optlev
    self%micro_opt      = src%micro_opt
    self%maxcycle       = src%maxcycle
    self%maxdispl_opt   = src%maxdispl_opt
    self%ethr_opt       = src%ethr_opt
    self%gthr_opt       = src%gthr_opt
    self%hlow_opt       = src%hlow_opt
    self%hmax_opt       = src%hmax_opt
    self%acc_opt        = src%acc_opt
    self%maxerise       = src%maxerise
    self%hguess         = src%hguess
    self%exact_rf       = src%exact_rf
    self%average_conv   = src%average_conv
    self%tsopt          = src%tsopt
    self%iupdat         = src%iupdat
    self%opt_engine     = src%opt_engine
    self%lbfgs_histsize = src%lbfgs_histsize
    self%hess_init      = src%hess_init
    self%logextxyz      = src%logextxyz

! ── smooth-function parameters ───────────────────────────────────────────────
    self%L       = src%L
    self%k       = src%k
    self%shift   = src%shift
    self%scaling = src%scaling

! ── printout and I/O ─────────────────────────────────────────────────────────
    self%pr_energies = src%pr_energies
    self%eout_unit   = src%eout_unit
    if (allocated(src%elog)) self%elog = src%elog

! ── ONIOM integer maps ───────────────────────────────────────────────────────
    if (allocated(src%ONIOMmap))    self%ONIOMmap    = src%ONIOMmap

! ── thermochemistry ──────────────────────────────────────────────────────────
    self%do_HR              = src%do_HR
    self%full_HR            = src%full_HR
    self%hu_steps           = src%hu_steps
    self%nt                 = src%nt
    self%ithr               = src%ithr
    self%fscal              = src%fscal
    self%sthr               = src%sthr
    self%initialize_hr_type = src%initialize_hr_type
    self%mh_type            = src%mh_type
    self%hr_hu_type         = src%hr_hu_type
    self%deform_opt_hess    = src%deform_opt_hess
    self%doh_stepsize       = src%doh_stepsize
    self%chess_id_guess     = src%chess_id_guess
    self%g_sampling         = src%g_sampling
    self%gs_hess_type       = src%gs_hess_type
    if (allocated(src%emodel))       self%emodel       = src%emodel
    if (allocated(src%temperatures)) self%temperatures = src%temperatures
    if (allocated(src%et))   self%et   = src%et
    if (allocated(src%ht))   self%ht   = src%ht
    if (allocated(src%gt))   self%gt   = src%gt
    if (allocated(src%stot)) self%stot = src%stot

! ── ONIOM model data ─────────────────────────────────────────────────────────
    !> lwoniom_data is plain (allocatable-only) Fortran data, so a deep copy
    !> is safe and REQUIRED: engrad consults calc%ONIOM of the object it is
    !> called with, and algos like crest_singlepoint work on a copy.
    !> (The ONIOMmols scratch buffers are recreated on demand in engrad.)
    if (allocated(src%ONIOM)) self%ONIOM = src%ONIOM

!>  NOTE: API handle objects (g0calc, ONIOMmols, chess) are NOT copied;
!>        they hold C-level or heavy reconstructed state and are re-initialized.
    return
  end subroutine calculation_copy

!=========================================================================================!
  subroutine calculation_set_freeze_range(self,nat,start,finish)
    class(calcdata) :: self
    integer,intent(in) :: nat,start,finish
    integer :: i,k
    if (allocated(self%freezelist)) deallocate (self%freezelist)
    allocate (self%freezelist(nat),source=.false.)
    k = 0
    do i = 1,nat

      if (i >= start.and.i <= finish) then
        k = k+1
        self%freezelist(i) = .true.
      end if
    end do
    self%nfreeze = k
  end subroutine calculation_set_freeze_range

  subroutine calculation_set_freeze_bools(self,freezetmp)
    class(calcdata) :: self
    logical,intent(in) :: freezetmp(:)
    integer :: nat
    if (allocated(self%freezelist)) deallocate (self%freezelist)
    nat = size(freezetmp,1)
    allocate (self%freezelist(nat),source=.false.)
    self%nfreeze = count(freezetmp)
    self%freezelist(:) = freezetmp(:)
  end subroutine calculation_set_freeze_bools

  subroutine calculation_freezegrad(self,grad)
    class(calcdata) :: self
    real(wp),intent(inout) :: grad(:,:)
    integer :: nat,i
    if (self%nfreeze > 0) then
      nat = size(grad,2)
      if (nat == self%nfreeze) then
        error stop '*** ERROR *** Must not freeze all atoms!'
      end if
      do i = 1,nat
        if (self%freezelist(i)) grad(:,i) = 0.0_wp
      end do
    end if
  end subroutine calculation_freezegrad

!=========================================================================================!

  subroutine calculation_set_charge(self,dchrg)
!***********************************************************
!* set the charge of all calculation_settings objects to 
!* the specified dchrg
!***********************************************************
    implicit none
    class(calcdata) :: self
    integer,intent(in) :: dchrg
    integer :: i,j
    if (self%ncalculations > 0) then
       j = dchrg
      do i = 1,self%ncalculations
        self%calcs(i)%chrg = j
      end do
    end if
    return
  end subroutine calculation_set_charge

  subroutine calculation_increase_charge(self,dchrg)
!******************************************************************
!* increase the charge of all calculation_settings objects by one
!* or the specified dchrg
!******************************************************************
    implicit none
    class(calcdata) :: self
    integer,intent(in),optional :: dchrg
    integer :: i,j
    if (self%ncalculations > 0) then
      if (present(dchrg)) then
        j = dchrg
      else
        j = 1
      end if
      do i = 1,self%ncalculations
        self%calcs(i)%chrg = self%calcs(i)%chrg+j
      end do
    end if
    return
  end subroutine calculation_increase_charge

!=========================================================================================!

  subroutine calculation_decrease_charge(self,dchrg)
!******************************************************************
!* decrease the charge of all calculation_settings objects by one
!* or the specified dchrg
!******************************************************************
    implicit none
    class(calcdata) :: self
    integer,intent(in),optional :: dchrg
    integer :: i,j
    if (self%ncalculations > 0) then
      if (present(dchrg)) then
        j = dchrg
      else
        j = 1
      end if
      do i = 1,self%ncalculations
        self%calcs(i)%chrg = self%calcs(i)%chrg-j
      end do
    end if
    return
  end subroutine calculation_decrease_charge

!=========================================================================================!

  subroutine calc_set_active(self,ids)
    implicit none
    class(calcdata) :: self
    integer,intent(in) :: ids(:)
    integer :: i,j,k,l
    if (allocated(self%weightbackup)) deallocate (self%weightbackup)
    if (allocated(self%activebackup)) deallocate (self%activebackup)
!>--- on-the-fly multiscale definition
    allocate (self%weightbackup(self%ncalculations),source=1.0_wp)
    allocate (self%activebackup(self%ncalculations),source=.true.)
    do i = 1,self%ncalculations
!>--- save backup weights
      self%weightbackup(i) = self%calcs(i)%weight
      self%activebackup(i) = self%calcs(i)%active
!>--- set the weight of all unwanted calculations to 0
      if (.not.any(ids(:) .eq. i)) then
        self%calcs(i)%weight = 0.0_wp
        self%calcs(i)%active = .false.
      else
!>--- and all other to active
        if (self%calcs(i)%weight == 0.0_wp) then
          self%calcs(i)%weight = 1.0_wp
        end if
        self%calcs(i)%active = .true.
      end if
    end do
  end subroutine calc_set_active

  subroutine calc_set_active_restore(self)
    implicit none
    class(calcdata) :: self
    integer :: i,j,k,l
    if (.not.allocated(self%weightbackup)) return
    if (.not.allocated(self%activebackup)) return
    do i = 1,self%ncalculations
!>--- set all to active and restore saved weights
      self%calcs(i)%weight = self%weightbackup(i)
      self%calcs(i)%active = self%activebackup(i)
      self%eweight(i) = self%weightbackup(i)
    end do
    deallocate (self%weightbackup)
  end subroutine calc_set_active_restore

!==========================================================================================!

  subroutine get_dipoles(calc,diptmp)
!*********************************************
!* Collect the x y and z dipole moments for
!* each calculation level in one array diptmp
!*********************************************
    implicit none
    type(calcdata),intent(inout) :: calc
    real(wp),allocatable,intent(out) :: diptmp(:,:)
    integer :: i,j,k,l,m

    m = calc%ncalculations
    allocate (diptmp(3,m),source=0.0_wp)
    do i = 1,m
      if (calc%calcs(i)%rddip) then
        diptmp(1:3,i) = calc%calcs(i)%dipole(1:3)
      end if
    end do
  end subroutine get_dipoles

!=========================================================================================!

  subroutine calculation_ONIOMexpand(self)
!*******************************************************
!* for an ONIOM calculation some of the calculators
!* have to be duplicated, which is done by this routine.
!* The canonical list of required calculations is taken
!* from lwONIOM's job table (lwoniom_get_jobs); crest
!* only clones its calculation_settings templates
!* accordingly and records the job -> calculator index
!* mapping in self%ONIOMmap.
!*******************************************************
    class(calcdata) :: self
    type(lwoniom_job),allocatable :: jobs(:)
    type(calculation_settings) :: calculator
    integer :: j,ref,newid,maxid
    character(len=40) :: atmp
    if (.not.allocated(self%ONIOM)) return

! ── canonical job list from lwONIOM (frag-major, high before low, root once) ──
    call lwoniom_get_jobs(self%ONIOM,jobs)

    maxid = maxval(jobs(:)%theoryid,1)
    if (maxid > self%ncalculations) then
      write (stdout,'(a)') '**ERROR** in ONIOM setup: not enough calculators defined!'
      error stop
    end if

    write (stdout,'(a)',advance='no') 'Assigning and duplicating calculators for ONIOM setup ...'
    flush (stdout)

! ── one crest calculator per job: reuse each template once, then clone ────────
    if (allocated(self%ONIOMmap)) deallocate (self%ONIOMmap)
    allocate (self%ONIOMmap(size(jobs)),source=0)
    do j = 1,size(jobs)
      !> the reference calculation_settings object (from layerlevel input)
      ref = jobs(j)%theoryid
      if (any(self%ONIOMmap(1:j-1) .eq. ref)) then
        !> template already assigned to another job: duplicate the calculator.
        !> However, we MUST not use restart I/O options for duplicates!
        call calculator%deallocate()
        calculator = self%calcs(ref)
        call calculator%norestarts()
        call self%add(calculator)
        newid = self%ncalculations
      else
        !> otherwise (i.e. it's not yet present), we can simply use it in place
        newid = ref
      end if
      self%ONIOMmap(j) = newid

      self%calcs(newid)%ONIOM_id = jobs(j)%fragid
      self%calcs(newid)%ONIOM_highlowroot = jobs(j)%level
      select case (jobs(j)%level)
      case (oniom_high)
        write (atmp,'(a,i0,a)') 'ONIOM.',jobs(j)%fragid,'.high'
      case (oniom_low)
        write (atmp,'(a,i0,a)') 'ONIOM.',jobs(j)%fragid,'.low'
      case (oniom_root)
        write (atmp,'(a,i0,a)') 'ONIOM.',jobs(j)%fragid,'.root'
      end select
      self%calcs(newid)%calcspace = trim(atmp)

      !> ALWAYS set the weight of ONIOM calcs to 0!
      self%calcs(newid)%weight = 0.0_wp

      !> transfer fragment charge/multiplicity, if attached
      if (allocated(jobs(j)%chrg)) self%calcs(newid)%chrg = jobs(j)%chrg
      if (allocated(jobs(j)%uhf)) self%calcs(newid)%uhf = jobs(j)%uhf
    end do
    !> NOTE: self%ONIOM%calcids remains untouched (it keeps the theory level ids)

    write (stdout,*) 'done.'
  end subroutine calculation_ONIOMexpand

!========================================================================================!

  subroutine calculation_info(self,iunit,printhdr)
    implicit none
    class(calcdata) :: self
    integer,intent(in) :: iunit
    logical,intent(in),optional :: printhdr
    integer :: i,j
    character(len=*),parameter :: fmt1 = '(1x,a20," : ",i5)'
    character(len=*),parameter :: fmt2 = '(1x,a20," : ",f12.5)'
    character(len=20) :: atmp
    integer :: constraintype(8)
    logical :: prhdr

    prhdr=.true.; if(present(printhdr)) prhdr = printhdr
    if(prhdr)then
      write (iunit,'(1x,a)') '----------------'
      write (iunit,'(1x,a)') 'Calculation info'
      write (iunit,'(1x,a)') '----------------'
    endif
    if (self%ncalculations <= 0) then
      write (iunit,'("> ",a)') 'No calculation levels set up!'
    else if (self%ncalculations > 1) then
      do i = 1,self%ncalculations
        if (self%calcs(i)%ONIOM_id > 0) then
          write (iunit,'("> ",a,i0,a)') 'Automatic ONIOM setup calculation level ',i,':'
        else
          write (iunit,'("> ",a,i0,a)') 'User-defined calculation level ',i,':'
        end if
        call self%calcs(i)%info(iunit)
      end do
    else
      write (iunit,'("> ",a)') 'User-defined calculation level:'
      call self%calcs(1)%info(iunit)
    end if
    write (iunit,*)

    if (self%nconstraints > 0) then
      write (iunit,'("> ",a)') 'User-defined constraints:'
      if (self%nconstraints <= 20) then
        do i = 1,self%nconstraints
          if (.not.self%cons(i)%active) cycle
          call self%cons(i)%print(iunit)
        end do
      else
        constraintype(:) = 0
        do i = 1,self%nconstraints
          if (.not.self%cons(i)%active) cycle
          j = self%cons(i)%type
          if (j > 0.and.j < 9) then
            constraintype(j) = constraintype(j)+1
          end if
        end do
        if (constraintype(1) > 0) write (iunit,'(2x,a,i0)') '# bond constraints    : ',constraintype(1)
        if (constraintype(2) > 0) write (iunit,'(2x,a,i0)') '# angle constraints   : ',constraintype(2)
        if (constraintype(3) > 0) write (iunit,'(2x,a,i0)') '# dihedral constraints: ',constraintype(3)
        if (constraintype(4) > 0) write (iunit,'(2x,a,i0)') '# wall potential      : ',constraintype(4)
        if (constraintype(5) > 0) write (iunit,'(2x,a,i0)') '# wall(logfermi) potential  : ',constraintype(5)
        if (constraintype(8) > 0) write (iunit,'(2x,a,i0)') '# bondrange constraints    : ',constraintype(8)
      end if
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
      if (allocated(self%ONIOM)) then
        write (iunit,'(1x,a)') 'ONIOM energy and gradient will be constructed from calculations:'
        do i = 1,self%ncalculations
          if (any(self%ONIOMmap(:) .eq. i)) then
            write (stdout,'(1x,i0)',advance='no') i
          end if
        end do
        write (iunit,*)
      end if
      write (iunit,*)
    end if

    return
  end subroutine calculation_info

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!> CALCULATION_SETTINGS associated routines
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!

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
    if (allocated(self%tblite)) deallocate (self%tblite)
    if (allocated(self%g0calc)) deallocate (self%g0calc)
    if (allocated(self%solv)) deallocate (self%solv)
    if (allocated(self%eeq)) deallocate (self%eeq)
    if (allocated(self%ff_dat)) deallocate (self%ff_dat)
    if (allocated(self%gff_fragments)) deallocate (self%gff_fragments)
    if (allocated(self%libpvol)) deallocate (self%libpvol)

    !> external callback: drop the references (we do not own the targets)
    self%ext_engrad => null()
    self%ext_userdata => null()

    self%id = 0
    self%prch = stdout
    self%chrg = 0
    self%uhf = 0
    self%epot = 0.0_wp
    self%efix = 0.0_wp
    self%etot = 0.0_wp

    self%rdwbo = .false.
    self%rddip = .false.
    self%dipole = 0.0_wp
    self%rddipgrad = .false.
    self%gradtype = 0
    self%gradfmt = 0

    self%tblitelvl = 2
    self%etemp = 300.0_wp
    self%etemp_user_set = .false.
    self%accuracy = 1.0_wp
    self%apiclean = .false.
    self%maxscc = 500
    self%saveint = .false.

    self%ngrid = 230
    self%extpressure = 0.0_wp
    self%proberad = 1.5_wp

    self%ONIOM_highlowroot = 0
    self%ONIOM_id = 0

    return
  end subroutine calculation_settings_deallocate

  subroutine calculation_settings_copy(self,src)
    implicit none
    class(calculation_settings),intent(out) :: self
    type(calculation_settings) :: src

! ── identity and printout ────────────────────────────────────────────────────
    self%id         = src%id
    self%prch       = src%prch
    self%pr         = src%pr
    self%prappend   = src%prappend
    self%prstdout   = src%prstdout
    self%refine_lvl = src%refine_lvl

! ── system ───────────────────────────────────────────────────────────────────
    self%chrg       = src%chrg
    self%uhf        = src%uhf
    self%multiplicity = src%multiplicity
    self%spin_polarized = src%spin_polarized
    self%active     = src%active
    self%weight     = src%weight
    self%threads    = src%threads

! ── allocatable strings ──────────────────────────────────────────────────────
    if (allocated(src%calcspace))      self%calcspace      = src%calcspace
    if (allocated(src%calcfile))       self%calcfile       = src%calcfile
    if (allocated(src%gradfile))       self%gradfile       = src%gradfile
    if (allocated(src%path))           self%path           = src%path
    if (allocated(src%other))          self%other          = src%other
    if (allocated(src%binary))         self%binary         = src%binary
    if (allocated(src%systemcall))     self%systemcall     = src%systemcall
    if (allocated(src%description))    self%description    = src%description
    if (allocated(src%shortflag))      self%shortflag      = src%shortflag
    if (allocated(src%gradkey))        self%gradkey        = src%gradkey
    if (allocated(src%efile))          self%efile          = src%efile
    if (allocated(src%solvmodel))      self%solvmodel      = src%solvmodel
    if (allocated(src%solvent))        self%solvent        = src%solvent
    if (allocated(src%parametrisation))self%parametrisation= src%parametrisation
    if (allocated(src%restartfile))    self%restartfile    = src%restartfile
    if (allocated(src%refgeo))         self%refgeo         = src%refgeo
    if (allocated(src%refcharges))     self%refcharges     = src%refcharges
    if (allocated(src%gff_fragments))  self%gff_fragments  = src%gff_fragments
    if (allocated(src%tbliteparam))    self%tbliteparam    = src%tbliteparam
    if (allocated(src%solv))           self%solv           = src%solv
    if (allocated(src%eeq))            self%eeq            = src%eeq

! ── gradient settings ────────────────────────────────────────────────────────
    self%numgrad    = src%numgrad
    self%gradstep   = src%gradstep
    self%rdgrad     = src%rdgrad
    self%gradtype   = src%gradtype
    self%gradfmt    = src%gradfmt

! ── property requests ────────────────────────────────────────────────────────
    self%rdwbo      = src%rdwbo
    self%rdqat      = src%rdqat
    self%dumpq      = src%dumpq
    self%rddip      = src%rddip
    self%dipole     = src%dipole
    self%rddipgrad  = src%rddipgrad
    self%getlmocent = src%getlmocent
    self%nprot      = src%nprot
    if (allocated(src%getsasa))   self%getsasa   = src%getsasa
    if (allocated(src%efield))    self%efield    = src%efield
    if (allocated(src%wbo))       self%wbo       = src%wbo
    if (allocated(src%qat))       self%qat       = src%qat
    if (allocated(src%dipgrad))   self%dipgrad   = src%dipgrad
    if (allocated(src%protxyz))   self%protxyz   = src%protxyz

! ── API / backend settings ───────────────────────────────────────────────────
    self%tblitelvl  = src%tblitelvl
    self%etemp          = src%etemp
    self%etemp_user_set = src%etemp_user_set
    self%accuracy   = src%accuracy
    self%apiclean   = src%apiclean
    self%maxscc     = src%maxscc
    self%saveint    = src%saveint
    self%ceh_guess  = src%ceh_guess
    self%restart    = src%restart

    self%ngrid       = src%ngrid
    self%extpressure = src%extpressure
    self%proberad    = src%proberad
    self%pvmodel     = src%pvmodel
    self%vdwset      = src%vdwset
    self%pvradscal   = src%pvradscal

    self%nconfig     = src%nconfig
    if (allocated(src%config))    self%config    = src%config
    if (allocated(src%occ))       self%occ       = src%occ

! ── ONIOM identifiers ────────────────────────────────────────────────────────
    self%ONIOM_highlowroot = src%ONIOM_highlowroot
    self%ONIOM_id          = src%ONIOM_id

! ── ORCA input block ─────────────────────────────────────────────────────────
    self%ORCA%mpi     = src%ORCA%mpi
    self%ORCA%nlines  = src%ORCA%nlines
    self%ORCA%maxcore = src%ORCA%maxcore
    self%ORCA%srckind = src%ORCA%srckind
    if (allocated(src%ORCA%cmd))   self%ORCA%cmd   = src%ORCA%cmd
    if (allocated(src%ORCA%input)) self%ORCA%input = src%ORCA%input

! ── inline potentials ────────────────────────────────────────────────────────
    self%ag      = src%ag
    self%penalty = src%penalty

! ── MLIP settings ────────────────────────────────────────────────────────────
    self%MPAR     = src%MPAR
    self%MPAR%iid = 0  !> reset instance ID for parallelization

! ── external callback ─────────────────────────────────────────────────────────
!>  Pointer-copy the host callback and its context: both copies refer to the
!>  same host-side routine/data (the targets are owned by the host program).
    self%ext_engrad    => src%ext_engrad
    self%ext_userdata  => src%ext_userdata

!>  NOTE: API handle objects (tblite, g0calc, ff_dat, libpvol) are NOT copied;
!>        they hold C-level state and are re-initialized on first use.
    return
  end subroutine calculation_settings_copy

!=========================================================================================!

  subroutine calculation_settings_set_external(self,fptr,userdata)
!*********************************************************************
!* Register an externally supplied energy+gradient routine on this
!* calculation level. This is the public entry point for a host
!* program that links CREST as a library: provide a routine matching
!* the engrad_interface and (optionally) an opaque context object,
!* and the level is switched to jobtype%external.
!*
!*  fptr     : procedure matching engrad_interface (required)
!*  userdata : optional opaque host context. Must have the TARGET or
!*             POINTER attribute on the caller side and must outlive
!*             all engrad calls; CREST only stores a reference to it
!*             and passes it back to fptr verbatim.
!*********************************************************************
    implicit none
    class(calculation_settings) :: self
    procedure(engrad_interface)              :: fptr
    class(*),pointer,intent(in),optional     :: userdata

    self%id = jobtype%external
    self%ext_engrad => fptr
    if (present(userdata)) then
      self%ext_userdata => userdata
    else
      self%ext_userdata => null()
    end if

    call self%autocomplete(self%id)
  end subroutine calculation_settings_set_external

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

  subroutine calculation_settings_norestarts(self)
!*************************************************
!* remove restart options from this calculation
!*************************************************
    implicit none
    class(calculation_settings) :: self
    self%restart = .false.
    if (allocated(self%refgeo)) deallocate (self%refgeo)
    if (allocated(self%restartfile)) deallocate (self%restartfile)
  end subroutine calculation_settings_norestarts

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
    call calculation_settings_shortflag(self)

    !> ── keep the spin multiplicity aligned with the uhf carrier ─────────────────
    call self%sync_multiplicity()

    if (.not.allocated(self%calcspace)) then
      !> I've decided to perform all calculations in a separate directory to
      !> avoid accumulation of files in the main workspace
      write (nmbr,'(i0)') id
      self%calcspace = 'calculation.level.'//trim(nmbr)
    end if

    if (self%pr.and.self%prch .ne. stdout) then
      self%prch = self%prch+id
    end if
  end subroutine calculation_settings_autocomplete

!========================================================================================!

  subroutine calc_sync_multiplicity(self)
    !***************************************************************
    !* Align the spin multiplicity with the uhf carrier:          *
    !* multiplicity = uhf + 1 = (Nα-Nβ) + 1 = 2S+1.               *
    !* uhf is the single source of truth that survives copies and *
    !* direct assignments; call this before any interface that    *
    !* consumes a multiplicity (ORCA, fmlip-relay).               *
    !***************************************************************
    implicit none
    class(calculation_settings) :: self
    self%multiplicity = self%uhf+1
  end subroutine calc_sync_multiplicity

!========================================================================================!

  subroutine calc_set_multiplicity(self,mult)
    !***************************************************************
    !* Set the spin state from a multiplicity (2S+1) value while  *
    !* keeping the uhf carrier (=2S=Nα-Nβ) aligned: uhf = mult-1. *
    !***************************************************************
    implicit none
    class(calculation_settings) :: self
    integer,intent(in) :: mult
    self%multiplicity = mult
    self%uhf = mult-1
  end subroutine calc_set_multiplicity

!>--- create a short calculation info flag
  subroutine calculation_settings_shortflag(self)
    implicit none
    class(calculation_settings) :: self
    integer :: i,j

    select case (self%id)
    case (jobtype%xtbsys)
      self%shortflag = 'xtb subprocess'
    case (jobtype%generic)
      self%shortflag = 'generic subprocess'
    case (jobtype%turbomole)
      self%shortflag = 'TURBOMOLE subprocess'
    case (jobtype%orca)
      self%shortflag = 'ORCA subprocess'
    case (jobtype%terachem)
      self%shortflag = 'TeraChem subprocess'
    case (jobtype%tblite)
      select case (self%tblitelvl)
      case (xtblvl%gfn2)
        self%shortflag = 'GFN2-xTB'
      case (xtblvl%gfn1)
        self%shortflag = 'GFN1-xTB'
      case (xtblvl%ipea1)
        self%shortflag = 'IPEA1-xTB'
      case (xtblvl%ceh)
        self%shortflag = 'CEH'
      case (xtblvl%eeq)
        self%shortflag = 'EEQ(D4)'
      case (xtblvl%param)
        self%shortflag = 'parameter file: '//trim(self%tbliteparam)
      end select
    case (jobtype%gfn0)
      self%shortflag = 'GFN0-xTB'
    case (jobtype%gfn0occ)
      self%shortflag = 'GFN0-xTB*'
    case (jobtype%gfnff)
      self%shortflag = 'GFN-FF'
    case (jobtype%libpvol)
      self%shortflag = 'libpvol'
    case (jobtype%lj)
      self%shortflag = 'LJ'
    case (jobtype%approxg)
      self%shortflag = 'model Hessian gradient'
    case (jobtype%penalty)
      self%shortflag = 'penalty potential'
    case (jobtype%mlip)
      self%shortflag = 'MLIP'
    case (jobtype%solvation)
      self%shortflag = 'solvation'
      if (allocated(self%solv)) then
        if (allocated(self%solv%smodel)) &
        & self%shortflag = self%shortflag//'/'//trim(self%solv%smodel)
        if (allocated(self%solv%solvent)) &
        & self%shortflag = self%shortflag//'('//trim(self%solv%solvent)//')'
      end if
    case (jobtype%electrostatic)
      self%shortflag = 'EEQ'
      if (allocated(self%eeq)) then
        if (allocated(self%eeq%charge_model)) &
        & self%shortflag = trim(self%eeq%charge_model)
      end if
    case (jobtype%external)
      self%shortflag = 'external (host callback)'
    case default
      self%shortflag = 'undefined'
    end select
    if (allocated(self%solvmodel).and.allocated(self%solvent)) then
      self%shortflag = self%shortflag//'/'//trim(self%solvmodel)
      self%shortflag = self%shortflag//'('//trim(self%solvent)//')'
    end if
  end subroutine calculation_settings_shortflag

!>-- generate a unique print id for the calculation
  subroutine calculation_settings_printid(self,thread,id)
    implicit none
    class(calculation_settings) :: self
    integer,intent(in)  :: thread,id
    integer :: i,j,dum
    character(len=50) :: nmbr
    dum = 100*(thread+1)
    dum = dum+id
    self%prch = dum
  end subroutine calculation_settings_printid

  subroutine calculation_dump_dipgrad(self,filename)
    implicit none
    class(calculation_settings) :: self
    character(len=*),intent(in) :: filename
    integer :: i,j,ich
    if (.not.allocated(self%dipgrad)) return
    open (newunit=ich,file=filename)
    do i = 1,size(self%dipgrad,2)
      write (ich,'(3f20.10)') self%dipgrad(1:3,i)
    end do
    close (ich)
  end subroutine calculation_dump_dipgrad

!=========================================================================================!

  subroutine calculation_settings_info(self,iunit)
    implicit none
    class(calculation_settings) :: self
    integer,intent(in) :: iunit
    integer :: i,j
    character(len=*),parameter :: fmt1 = '(" :",2x,a20," : ",i0)'
    character(len=*),parameter :: fmt2 = '(" :",2x,a20," : ",f0.5)'
    character(len=*),parameter :: fmt3 = '(" :",2x,a20," : ",a)'
    character(len=*),parameter :: fmt4 = '(" :",1x,a)'
    character(len=20) :: atmp
    logical :: gxtbwarn

    gxtbwarn = .false.

    if (allocated(self%description)) then
      write (iunit,'(" :",1x,a)') trim(self%description)
    else
      write (atmp,*) 'Job type'
      write (iunit,fmt1) atmp,self%id
    end if
    !> more info
    if (self%id == jobtype%tblite) then
      select case (self%tblitelvl)
      case (xtblvl%gfn2)
        write (iunit,fmt4) 'GFN2-xTB level'
      case (xtblvl%gfn1)
        write (iunit,fmt4) 'GFN1-xTB level'
      case (xtblvl%ceh)
        write (iunit,fmt4) 'Charge Extended Hückel (CEH) model'
      end select
    end if

    !> MLIP (fmlip-relay) backend details — print only what is meaningful
    !> for the selected backend
    if (self%id == jobtype%mlip) then
      block
        character(len=:),allocatable :: bk
        bk = 'unknown'
        if (allocated(self%MPAR%backend)) bk = trim(self%MPAR%backend)
        ! ── friendly backend headline ───────────────────────────────────────
        select case (bk)
        case ('uma')
          write (iunit,fmt4) 'FairChem UMA foundation model (fairchem)'
        case ('mace_off')
          write (iunit,fmt4) 'MACE-OFF23 organic force field'
        case ('mace_mp')
          write (iunit,fmt4) 'MACE-MP foundation model (Materials Project)'
        case ('mace')
          write (iunit,fmt4) 'Custom MACE model'
        case ('lj')
          write (iunit,fmt4) 'Lennard-Jones potential'
        case ('dummy')
          write (iunit,fmt4) 'Dummy (zero) potential'
        case default
          write (iunit,fmt4) 'fmlip-relay backend: '//bk
        end select
        write (atmp,*) 'MLIP backend'
        write (iunit,fmt3) atmp,bk
        ! ── backend-specific model / task selection ─────────────────────────
        select case (bk)
        case ('uma')
          write (atmp,*) 'UMA model'
          if (allocated(self%MPAR%umamodel)) then
            write (iunit,fmt3) atmp,trim(self%MPAR%umamodel)
          else
            write (iunit,fmt3) atmp,'uma-s-1p2 (default)'
          end if
          write (atmp,*) 'UMA task'
          if (allocated(self%MPAR%umatask)) then
            write (iunit,fmt3) atmp,trim(self%MPAR%umatask)
          else
            write (iunit,fmt3) atmp,'omol (default)'
          end if
        case ('mace_off','mace_mp')
          if (allocated(self%MPAR%modelpath)) then
            write (atmp,*) 'MACE model file'
            write (iunit,fmt3) atmp,trim(self%MPAR%modelpath)
          else
            write (atmp,*) 'MACE model size'
            if (allocated(self%MPAR%modelsize)) then
              write (iunit,fmt3) atmp,trim(self%MPAR%modelsize)
            else
              write (iunit,fmt3) atmp,'medium (default)'
            end if
          end if
        case default
          if (allocated(self%MPAR%modelpath)) then
            write (atmp,*) 'Model file'
            write (iunit,fmt3) atmp,trim(self%MPAR%modelpath)
          end if
        end select
        ! ── compute device: always shown for the torch-based backends ───────
        select case (bk)
        case ('uma','mace_off','mace_mp','mace')
          write (atmp,*) 'Compute device'
          if (allocated(self%MPAR%device)) then
            write (iunit,fmt3) atmp,trim(self%MPAR%device)
          else
            write (iunit,fmt3) atmp,'cpu (default)'
          end if
        case default
          !> lj / dummy etc. have no torch device; show only if explicitly set
          if (allocated(self%MPAR%device)) then
            write (atmp,*) 'Compute device'
            write (iunit,fmt3) atmp,trim(self%MPAR%device)
          end if
        end select
        if (self%threads > 0) then
          write (atmp,*) 'Server thread cap'
          write (iunit,fmt1) atmp,self%threads
        end if
        if (self%MPAR%BASE_PORT /= 54320) then
          write (atmp,*) 'Socket base port'
          write (iunit,fmt1) atmp,self%MPAR%BASE_PORT
        end if
        if (self%MPAR%TIMEOUT_SEC /= 120) then
          write (atmp,*) 'Server timeout [s]'
          write (iunit,fmt1) atmp,self%MPAR%TIMEOUT_SEC
        end if
      end block
    end if

    !> composite implicit-solvation calculator details
    if (self%id == jobtype%solvation .and. allocated(self%solv)) then
      write (iunit,fmt4) 'Composite implicit solvation (ddX-based)'
      if (allocated(self%solv%solvent)) then
        write (atmp,*) 'solvent'
        write (iunit,fmt3) atmp,trim(self%solv%solvent)
      end if
      if (allocated(self%solv%smodel)) then
        write (atmp,*) 'continuum model'
        write (iunit,fmt3) atmp,trim(self%solv%smodel)
      end if
      if (allocated(self%solv%charge_model)) then
        write (atmp,*) 'charge model'
        write (iunit,fmt3) atmp,trim(self%solv%charge_model)
      end if
      write (atmp,*) 'nonpolar params'
      write (iunit,fmt3) atmp,'GFN2/ALPB'
      write (atmp,*) 'H-bond term'
      if (self%solv%do_hbond) then
        write (iunit,fmt3) atmp,'on'
      else
        write (iunit,fmt3) atmp,'off'
      end if
    end if

    !> charge-equilibration electrostatics details
    if (self%id == jobtype%electrostatic .and. allocated(self%eeq)) then
      write (iunit,fmt4) 'Charge-equilibration electrostatics (multicharge)'
      if (allocated(self%eeq%charge_model)) then
        write (atmp,*) 'charge model'
        write (iunit,fmt3) atmp,trim(self%eeq%charge_model)
      end if
    end if

    !> externally supplied (host-program) potential details
    if (self%id == jobtype%external) then
      write (iunit,fmt4) 'Externally supplied potential (host callback)'
      write (atmp,*) 'callback'
      if (associated(self%ext_engrad)) then
        write (iunit,fmt3) atmp,'Fortran procedure'
      else
        write (iunit,fmt3) atmp,'NOT associated (!)'
      end if
      write (atmp,*) 'user context'
      if (associated(self%ext_userdata)) then
        write (iunit,fmt3) atmp,'provided'
      else
        write (iunit,fmt3) atmp,'none'
      end if
    end if

    if (any((/jobtype%orca,jobtype%xtbsys,jobtype%turbomole, &
    &  jobtype%generic,jobtype%terachem/) == self%id)) then
      if (index(self%binary,'gxtb') .ne. 0) then
        write (iunit,fmt4) 'g-xTB (development version)'
        gxtbwarn = .true.
      else
        write (iunit,'(" :",3x,a,a)') 'selected binary : ',trim(self%binary)
      end if
    end if
    if (self%refine_lvl > 0) then
      write (atmp,*) 'refinement stage'
      block
        character(len=20) :: rtmp
        select case (self%refine_lvl)
        case (1);  rtmp = 'singlepoint'
        case (2);  rtmp = 'correction'
        case (3);  rtmp = 'geoopt'
        case (5);  rtmp = 'deltaG'
        case (10); rtmp = 'post_opt'
        case (11); rtmp = 'post_sp'
        case (12); rtmp = 'post_reopt'
        case default
          write (rtmp,'(i0)') self%refine_lvl
        end select
        write (iunit,fmt3) atmp,trim(rtmp)
      end block
    end if

    !> system data
    write (atmp,*) 'Molecular charge'
    write (iunit,fmt1) atmp,self%chrg
    !> spin: ORCA and the UMA (omol) backend consume a multiplicity (2S+1)
    !> from the multiplicity store, the xtb-type methods consume uhf = Nα-Nβ
    block
      logical :: use_mult
      use_mult = (self%id == jobtype%orca)
      if (self%id == jobtype%mlip .and. allocated(self%MPAR%backend)) then
        if (trim(self%MPAR%backend) == 'uma') use_mult = .true.
      end if
      if (use_mult) then
        write (atmp,*) 'Multiplicity'
        write (iunit,fmt1) atmp,self%multiplicity
      else if (self%uhf /= 0) then
        write (atmp,*) 'UHF parameter'
        write (iunit,fmt1) atmp,self%uhf
      end if
    end block
    if (self%id == jobtype%tblite .and. self%spin_polarized) then
      write (atmp,*) 'Spin-polarization'
      write (iunit,fmt3) atmp,'yes'
    end if

    if (allocated(self%solvmodel)) then
      write (atmp,*) 'Solvation model'
      write (iunit,fmt3) atmp,trim(self%solvmodel)
    end if
    if (allocated(self%solvent)) then
      write (atmp,*) 'Solvent'
      write (iunit,fmt3) atmp,trim(self%solvent)
    end if

    !> xTB specific parameters
    if (any((/jobtype%tblite,jobtype%xtbsys,jobtype%gfn0,jobtype%gfn0occ/) == self%id)) then
      write (atmp,*) 'Fermi temperature'
      write (iunit,fmt2) atmp,self%etemp
      write (atmp,*) 'Accuracy'
      write (iunit,fmt2) atmp,self%accuracy
      write (atmp,*) 'max SCC cycles'
      write (iunit,fmt1) atmp,self%maxscc
    end if

    write (atmp,*) 'Reset data?'
    if (self%apiclean) write (iunit,fmt3) atmp,'yes'
    write (atmp,*) 'Read WBOs?'
    if (self%rdwbo) write (iunit,fmt3) atmp,'yes'
    write (atmp,*) 'Read dipoles?'
    if (self%rddip) write (iunit,fmt3) atmp,'yes'

    if (self%ONIOM_highlowroot /= 0) then
      select case (self%ONIOM_highlowroot)
      case (1)
        write (atmp,*) 'ONIOM frag ("high")'
      case (2)
        write (atmp,*) 'ONIOM frag ("low")'
      case (3)
        write (atmp,*) 'ONIOM frag ("root")'
      end select
      write (iunit,fmt1) trim(atmp),self%ONIOM_id
    else
      if (self%weight .ne. 1.0_wp) then
        write (atmp,*) 'Weight'
        write (iunit,fmt2) atmp,self%weight
      end if
    end if

    if (gxtbwarn) then
      write (iunit,fmt4) 'WARNING: This currently is the development version of g-xTB.'
      write (iunit,fmt4) 'WARNING: Gradients are NUMERICAL (i.e., expensive and noisy!)'
    end if

  end subroutine calculation_settings_info

!=========================================================================================!

  subroutine create_calclevel_shortcut(self,levelstring, &
      & chrg,uhf,solvmodel,solvent)
!*********************************************************************
!* subroutine create_calclevel_shortcut called with %create(...)
!* Set up a calculation_settings object for a given level of theory
!* More shortcuts can be added as required.
!*
!* Optional settings are for:
!*  - molecular charge  (integer)
!*  - uhf parameter (integer)
!*  - solvent/solventmodel (either none or BOTH must be present to work)
!*
!* Be careful about the intent(out) setting!
!* Also, the routine is "dumb" and does not check if the user-provided
!* settings actually make sense for a create_calclevel_shortcutation. It very much
!* exists as an internal code shortcut only.
!*********************************************************************
    implicit none
    class(calculation_settings),intent(out) :: self
    character(len=*),intent(in) :: levelstring
    integer,intent(in),optional :: chrg
    integer,intent(in),optional :: uhf
    character(len=*),intent(in),optional :: solvmodel
    character(len=*),intent(in),optional :: solvent
    call self%deallocate()
    select case (trim(levelstring))
    case ('gfnff','--gff','--gfnff')
      self%id = jobtype%gfnff
    case ('gfn0','--gfn0')
      self%id = jobtype%gfn0
    case ('gfn2','--gfn2')
      self%id = jobtype%tblite
      self%tblitelvl = 2
    case ('gfn1','--gfn1')
      self%id = jobtype%tblite
      self%tblitelvl = 1
    case ('gp3')
      self%id = jobtype%turbomole
      self%rdgrad = .false.
      self%binary = 'gp3'
    case ('gxtb','--gxtb')
      if (have_gxtb) then
        self%id = jobtype%tblite
        self%tblitelvl = xtblvl%gxtb
        self%etemp = 0.0_wp  ! g-xTB uses integer occupations (T=0)
      else
        self%id = jobtype%xtbsys
        self%other = '--gxtb'
      end if
    case ('orca','--orca')
      self%id = jobtype%orca

    case ('uma','--uma')
      !> FairChem UMA via fmlip-relay; default to the molecular (omol) task head
      self%id = jobtype%mlip
      self%MPAR%backend = 'uma'
      self%MPAR%umatask = 'omol'

    case ('maceoff','mace-off','mace_off','--maceoff')
      !> MACE-OFF23 organic force field via fmlip-relay; default to medium size
      self%id = jobtype%mlip
      self%MPAR%backend = 'mace_off'
      self%MPAR%modelsize = 'medium'

    case ('generic')
      self%id = jobtype%generic

    case ('external','--external')
      !> the host program must still register its callback via %set_external
      self%id = jobtype%external

    end select

    if (present(chrg)) then
      self%chrg = chrg
    end if

    if (present(uhf)) then
      self%uhf = uhf
    end if

    !> both must be present to work
    if (present(solvmodel).and.present(solvent)) then
      !> the first two if-cases exist to convert cli args
      !> into sensible keywords (required for legacy compatibility)
      if (index(solvmodel,'gbsa') .ne. 0) then
        self%solvmodel = 'gbsa'
      else if (index(solvmodel,'alpb') .ne. 0) then
        self%solvmodel = 'alpb'
      else
        self%solvmodel = trim(solvmodel)
      end if
      self%solvent = trim(solvent)
    end if

    call self%autocomplete(self%id)
  end subroutine create_calclevel_shortcut

!=========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!=========================================================================================!
end module calc_type
