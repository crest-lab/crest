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
  use strucrd,only:coord
!>--- api types
  use tblite_api
  use gfn0_api
  use gfnff_api,only:gfnff_data
  use xhcff_api,only:xhcff_calculator
!>--- other types
  use orca_type
  use lwoniom_module
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
    integer :: xhcff     = 10
    integer :: lj        = 11
  end type enum_jobtype
  type(enum_jobtype), parameter,public :: jobtype = enum_jobtype()

  character(len=45),parameter,private :: jobdescription(12) = [ &
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
     & 'XHCFF calculation via XHCFF-lib             ', &
     & 'Lennard-Jones potential calculation         ' ]
!&>

!=========================================================================================!
!>--- data object that contains the data for a *SINGLE* calculation
  public :: calculation_settings
  type :: calculation_settings

    integer :: id = 0         !> calculation type (see "jobtype" parameter above)
    integer :: prch = stdout  !> printout channel
    logical :: pr = .false.   !> allow the calculation to produce printout? Results in a lot I/O
    logical :: prappend = .false. !> append printout 
    integer :: refine_lvl = 0 !> to allow defining different refinement levels

    integer :: chrg = 0          !> molecular charge
    integer :: uhf = 0           !> uhf parameter (xtb) or multiplicity (other)
    logical :: active = .true.   !> active setting to disable the calculation (this is different from weight=0)
    real(wp) :: weight = 1.0_wp  !> calculation weight (when adding them up)

    character(len=:),allocatable :: calcspace  !> subdirectory to perform the calculation in
    character(len=:),allocatable :: calcfile
    character(len=:),allocatable :: gradfile
    character(len=:),allocatable :: path
    character(len=:),allocatable :: other
    character(len=:),allocatable :: binary      !> binary or generic script
    character(len=:),allocatable :: systemcall  !> systemcall for running generic scripts
    character(len=:),allocatable :: description !> see above jobdescription parameter

!>--- gradient format specifications
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
    type(tblite_data),allocatable :: tblite

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
    real(wp) :: proberad = 1.5_wp      !>  proberadius in Angstroem
    integer :: vdwset = 0              !>  Set of VDW radii to use in sas calculation -> default D3, 1 -> Bondi
    type(xhcff_calculator),allocatable :: xhcff

    !> ONIOM fragment IDs
    integer :: ONIOM_highlowroot = 0 
    integer :: ONIOM_id = 0

    !> ORCA job template
    type(orca_input) :: ORCA

!>--- Type procedures
  contains
    procedure :: deallocate => calculation_settings_deallocate
    procedure :: addconfig => calculation_settings_addconfig
    procedure :: autocomplete => calculation_settings_autocomplete
    procedure :: printid => calculation_settings_printid
    procedure :: info => calculation_settings_info
    procedure :: create => create_calclevel_shortcut
  end type calculation_settings
!=========================================================================================!

!=========================================================================================!
!> data object that collects settings for *ALL* calculations and constraints.
  public :: calcdata
  type :: calcdata

    integer :: id = 0  !> this parameter will decide how to return or add up energies and gradients
    integer :: refine_stage = 0 !> to allow iterating different refinement stages

!>--- calculations
    integer :: ncalculations = 0
    type(calculation_settings),allocatable :: calcs(:)
    real(wp),allocatable :: etmp(:)
    real(wp),allocatable :: grdtmp(:,:,:)
    real(wp),allocatable :: eweight(:)
    real(wp),allocatable :: weightbackup(:)
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
    integer  :: optlev = 0
    integer  :: micro_opt = 20
    integer  :: maxcycle = 0
    real(wp) :: maxdispl_opt = -1.0_wp
    real(wp) :: ethr_opt = -1.0_wp  !> ΔE convergence
    real(wp) :: gthr_opt = -1.0_wp  !> RMS force convergence
    real(wp) :: hlow_opt = 0.010_wp
    real(wp) :: hmax_opt = 5.0_wp
    real(wp) :: acc_opt = 1.0_wp
    real(wp) :: maxerise = 1.0e-5_wp
    logical  :: exact_rf = .false.
    logical  :: average_conv = .false.
    logical  :: tsopt = .false.
    integer  :: iupdat = 0  !> 0=BFGS, 1=Powell, 2=SR1, 3=Bofill, 4=Schlegel
    integer  :: opt_engine = 0 !> default: ANCOPT

!>--- GFN0* data, needed for special MECP application
    type(gfn0_data),allocatable  :: g0calc

!>--- printouts and IO
    logical :: pr_energies = .false.
    integer :: eout_unit = stdout
    character(len=:),allocatable :: elog

!>--- ONIOM calculator data
    type(lwoniom_data),allocatable :: ONIOM
    type(coord),allocatable :: ONIOMmols(:)
    integer,allocatable :: ONIOMmap(:) !> map ONIOM fragments to calculation_settings
    integer,allocatable :: ONIOMrevmap(:) !> map calculation settings to ONIOM frags (or zero)

!>--- Type procedures
  contains
    procedure :: reset => calculation_reset
    procedure :: init => calculation_init
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
    procedure :: freezegrad => calculation_freezegrad
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
    self%maxdispl_opt = -1.0_wp
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
    !if (allocated(self%wfn)) deallocate (self%wfn)
    !if (allocated(self%tbcalc)) deallocate (self%tbcalc)
    !if (allocated(self%ctx)) deallocate (self%ctx)
    !if (allocated(self%tbres)) deallocate (self%tbres)
    if (allocated(self%tblite)) deallocate(self%tblite)
    if (allocated(self%g0calc)) deallocate (self%g0calc)
    if (allocated(self%ff_dat)) deallocate (self%ff_dat)
    if (allocated(self%xhcff)) deallocate (self%xhcff)

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

    self%ONIOM_highlowroot = 0
    self%ONIOM_id = 0
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
  subroutine calculation_copy(self,src)
    class(calcdata) :: self
    type(calcdata) :: src
    integer :: i

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

  subroutine calculation_freezegrad(self,grad)
    class(calcdata) :: self
    real(wp),intent(inout) :: grad(:,:)
    integer :: nat,i
    if (self%nfreeze > 0) then
       nat = size(grad,2)
       if(nat == self%nfreeze)then
         error stop '**ERROR** Must not freeze all atoms!'
       endif 
       do i =1,nat
         if(self%freezelist(i)) grad(:,i) = 0.0_wp
       enddo
    end if
  end subroutine calculation_freezegrad


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

    if (self%pr) then
      self%prch = self%prch+id
    end if
  end subroutine calculation_settings_autocomplete

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

!=========================================================================================!

  subroutine calculation_ONIOMexpand(self)
!*******************************************************
!* for an ONIOM calculations some of the calculators
!* have to be duplikated, which is done by this routine
!*******************************************************
    class(calcdata) :: self
    integer :: ncalcs
    integer :: maxid
    integer :: i,j,k,l,newid,j2
    type(calculation_settings) :: calculator
    integer,allocatable :: newids(:,:)
    character(len=40) :: atmp
    if (.not.allocated(self%ONIOM)) return
    ncalcs = self%ONIOM%ncalcs
    maxid = maxval(self%ONIOM%calcids(1,:),1)
    if (maxid > self%ncalculations) then
      write (stdout,'(a)') '**ERROR** in ONIOM setup: not enough calculators defined!'
      error stop
    end if

    write(stdout,'(a)',advance='no') 'Assigning and duplicating calculators for ONIOM setup ...'
    flush(stdout)

    allocate (self%ONIOMmap(ncalcs),source=0)
    allocate (newids(2,self%ONIOM%nfrag), source = 0)
    k = 0
    do i = 1,self%ONIOM%nfrag

      do l = 1,2
        !> j is now the ID of the reference calculation_settings object
        j = self%ONIOM%calcids(l,i)
        if (l == 2) then
          !> to exlcude the highest ONIOM layer, we need to cycle
          j2 = self%ONIOM%calcids(1,i)
          if (j == j2)then
             newid = newids(1,i)
             newids(2,i) = newid
             self%calcs(newid)%ONIOM_highlowroot = 3
             self%calcs(newid)%ONIOM_id = i
             cycle
          endif
        end if

        if (any(self%ONIOMmap(:) .eq. j)) then
          !> If one of this type is already in the mapping, duplicate the calculator and add it
          call calculator%deallocate()
          calculator = self%calcs(j)
          call self%add(calculator)
          newid = self%ncalculations
          k = k+1
          self%ONIOMmap(k) = newid
        else
          !> otherwise (i.e. it's not yet present), we can simply add it
          k = k+1
          newid = j
          self%ONIOMmap(k) = newid
        end if
        newids(l,i) = newid

        self%calcs(newid)%ONIOM_highlowroot = l
        self%calcs(newid)%ONIOM_id = i
        select case(l)
        case( 1 ) 
        write(atmp,'(a,i0,a)') 'ONIOM.',i,'.high'
        case( 2 )
        write(atmp,'(a,i0,a)') 'ONIOM.',i,'.low'
        case( 3 )
        write(atmp,'(a,i0,a)') 'ONIOM.',i,'.root'
        end select
        self%calcs(newid)%calcspace = trim(atmp)

        !> ALWAYS set the weight of ONIOM calcs to 0!
        self%calcs(newid)%weight = 0.0_wp
      end do
    end do
    self%ONIOM%calcids = newids
    deallocate(newids)

    if(allocated(self%ONIOMrevmap)) deallocate(self%ONIOMrevmap)
    allocate(self%ONIOMrevmap( self%ncalculations ), source=0)
    do i =1,self%ncalculations
       do j =1,self%ONIOM%ncalcs
         if( self%ONIOMmap(j) == i)then
            self%ONIOMrevmap(i) = j
         endif
       enddo
    enddo

    write(stdout,*) 'done.'
  end subroutine calculation_ONIOMexpand

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

    if (allocated(self%description)) then
      write (iunit,'(" :",1x,a)') trim(self%description)
    else
      write (atmp,*) 'Job type'
      write (iunit,fmt1) atmp,self%id
    end if
    !> more info
    if (self%id == jobtype%tblite) then
      select case (self%tblitelvl)
      case (2)
        write (iunit,fmt4) 'GFN2-xTB level'
      case (1)
        write (iunit,fmt4) 'GFN1-xTB level'
      end select
    end if
    if (any((/jobtype%orca,jobtype%xtbsys,jobtype%turbomole, &
    &  jobtype%generic,jobtype%terachem/) == self%id)) then
      write (iunit,'(" :",3x,a,a)') 'selected binary : ',trim(self%binary)
    end if
    if (self%refine_lvl > 0) then
      write (atmp,*) 'refinement stage'
      write (iunit,fmt1) atmp,self%refine_lvl
    end if

    !> system data
    write (atmp,*) 'Molecular charge'
    write (iunit,fmt1) atmp,self%chrg
    if (self%uhf /= 0) then
      write (atmp,*) 'UHF parameter'
      write (iunit,fmt1) atmp,self%uhf
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

    write (atmp,*) 'Weight'
    write (iunit,fmt2) atmp,self%weight

    if(self%ONIOM_highlowroot /= 0)then
      select case(self%ONIOM_highlowroot)
      case( 1 )
      write (atmp,*) 'ONIOM frag ("high")'
      case(2)
      write (atmp,*) 'ONIOM frag ("low")'
      case(3)
      write (atmp,*) 'ONIOM frag ("root")'
      end select
      write (iunit,fmt1) trim(atmp),self%ONIOM_id
    endif

  end subroutine calculation_settings_info

!========================================================================================!

  subroutine calculation_info(self,iunit)
    implicit none
    class(calcdata) :: self
    integer,intent(in) :: iunit
    integer :: i,j
    character(len=*),parameter :: fmt1 = '(1x,a20," : ",i5)'
    character(len=*),parameter :: fmt2 = '(1x,a20," : ",f12.5)'
    character(len=20) :: atmp
    integer :: constraintype(8)

    write (iunit,'(1x,a)') '----------------'
    write (iunit,'(1x,a)') 'Calculation info'
    write (iunit,'(1x,a)') '----------------'
    if (self%ncalculations <= 0) then
      write (iunit,'("> ",a)') 'No calculation levels set up!'
    else if (self%ncalculations > 1) then
      do i = 1,self%ncalculations
        if(self%calcs(i)%ONIOM_id > 0 )then
        write (iunit,'("> ",a,i0,a)') 'Automatic ONIOM setup calculation level ',i,':'
        else
        write (iunit,'("> ",a,i0,a)') 'User-defined calculation level ',i,':'
        endif
        call self%calcs(i)%info(iunit)
      end do
    else
      write (iunit,'("> ",a)') 'User-defined calculation level:'
      call self%calcs(1)%info(iunit)
    end if
    write (iunit,*)

    if (self%nconstraints > 0) then
      write (iunit,'("> ",a)') 'User-defined constraints:'
      if(self%nconstraints <= 20)then
        do i = 1,self%nconstraints
          call self%cons(i)%print(iunit)
        end do
      else
         constraintype(:) = 0
         do i = 1,self%nconstraints
          j = self%cons(i)%type
          if(j > 0 .and. j < 9)then
            constraintype(j) = constraintype(j) + 1
          endif 
         end do
         if(constraintype(1) > 0) write (iunit,'(2x,a,i0)') '# bond constraints    : ',constraintype(1)
         if(constraintype(2) > 0) write (iunit,'(2x,a,i0)') '# angle constraints   : ',constraintype(2)
         if(constraintype(3) > 0) write (iunit,'(2x,a,i0)') '# dihedral constraints: ',constraintype(3)
         if(constraintype(4) > 0) write (iunit,'(2x,a,i0)') '# wall potential      : ',constraintype(4)
         if(constraintype(5) > 0) write (iunit,'(2x,a,i0)') '# wall(logfermi) potential  : ',constraintype(5)
         if(constraintype(8) > 0) write (iunit,'(2x,a,i0)') '# bondrange constraints    : ',constraintype(8)
      endif 
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
      if(allocated(self%ONIOM))then
        write (iunit,'(1x,a)') 'ONIOM energy and gradient will be constructed from calculations:'
        do i=1,self%ncalculations
          if(any(self%ONIOMmap(:).eq.i))then
             write(stdout,'(1x,i0)',advance='no') i
          endif
        enddo
        write (iunit,*)
      endif
      write (iunit,*)
    end if

    return
  end subroutine calculation_info

!=========================================================================================!

  subroutine create_calclevel_shortcut(self,levelstring)
!*********************************************************************
!* subroutine create_calclevel_shortcut called with %create(...)
!* Set up a calculation_settings object for a given level of theory
!* More shortcuts can be added as required.
!* Be careful about the intent(out) setting!
!*********************************************************************
    implicit none
    class(calculation_settings),intent(out) :: self
    character(len=*) :: levelstring
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
    case ('orca')
      self%id = jobtype%orca

    case ('generic')
      self%id = jobtype%generic

    end select
  end subroutine create_calclevel_shortcut

!=========================================================================================!

  subroutine calc_set_active(self,ids)
     implicit none
     class(calcdata) :: self
     integer,intent(in) :: ids(:) 
     integer :: i,j,k,l
     if(allocated(self%weightbackup)) deallocate(self%weightbackup)
!>--- on-the-fly multiscale definition
      allocate(self%weightbackup(self%ncalculations), source = 1.0_wp)
      do i=1,self%ncalculations
!>--- save backup weights
        self%weightbackup(i) =  self%calcs(i)%weight
!>--- set the weight of all unwanted calculations to 0
        if(.not.any(ids(:).eq.i))then
           self%calcs(i)%weight = 0.0_wp
           self%calcs(i)%active = .false.
        else
!>--- and all other to 1
           self%calcs(i)%weight = 1.0_wp
           self%calcs(i)%active = .true.
        endif
      enddo
  end subroutine calc_set_active
  
  subroutine calc_set_active_restore(self)
     implicit none
     class(calcdata) :: self
     integer :: i,j,k,l
     if(.not.allocated(self%weightbackup)) return
     do i=1,self%ncalculations
!>--- set all to active and restore saved weights        
        self%calcs(i)%weight = self%weightbackup(i)
        self%calcs(i)%active = .true.
        self%eweight(i) = self%weightbackup(i)
     enddo
     deallocate(self%weightbackup)
  end subroutine calc_set_active_restore
 

!=========================================================================================!
end module calc_type
