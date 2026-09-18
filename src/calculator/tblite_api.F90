!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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
! module tblite_api
! An interface of CREST to tblite
!====================================================!

module tblite_api
  use crest_parameters
  use strucrd
#ifdef WITH_TBLITE
  use mctc_env,only:error_type
  use mctc_io,only:structure_type,new
  use tblite_context_type,only:tblite_ctx => context_type
  use tblite_wavefunction_type,only:wavefunction_type,new_wavefunction
  use tblite_wavefunction,only:sad_guess,eeq_guess,shell_partition
  use tblite_xtb,xtb_calculator => xtb_calculator
  use tblite_xtb_calculator,only:new_xtb_calculator
#ifdef WITH_GXTB
  use tblite_xtb,only:new_gxtb_calculator
#endif
  use tblite_param,only:param_record
  use tblite_results,only:tblite_resultstype => results_type
  use tblite_wavefunction_mulliken,only:get_molecular_dipole_moment
  use tblite_ceh_singlepoint,only:ceh_singlepoint
  use tblite_ceh_ceh,only:new_ceh_calculator
  use tblite_spin,only:spin_polarization,new_spin_polarization
  use tblite_container,only:container_type
#endif
  use wiberg_mayer
  implicit none
  private

#ifndef WITH_TBLITE
  !> these are placeholders if no tblite is used!
  type :: wavefunction_type
    integer :: id = 0
  end type wavefunction_type
  type :: xtb_calculator
    integer :: id = 0
  end type xtb_calculator
  type :: tblite_ctx
    integer :: unit = stdout
    integer :: verbosity = 0
  end type tblite_ctx
  type :: tblite_resultstype
    integer :: id = 0
  end type tblite_resultstype
  type :: tblite_solvation_type
    integer :: id = 0
  end type tblite_solvation_type
#endif

!>--- tblite calculator bundle
  type :: tblite_data
    integer  :: lvl = 0
    real(wp) :: accuracy = 1.0_wp
    character(len=:),allocatable :: paramfile
    type(wavefunction_type)              :: wfn
    type(wavefunction_type),allocatable  :: wfn_aux
    type(xtb_calculator)        :: calc
    type(tblite_ctx)            :: ctx
    type(tblite_resultstype)    :: res
    logical :: spin_polarized = .false.
  end type tblite_data
  public :: tblite_data

  !> Type enumerator
  type :: enum_tblite_method
    integer :: unknown = 0
    integer :: gfn1 = 1
    integer :: gfn2 = 2
    integer :: ipea1 = 3
    !> the guesses can be used for charges, but NOT for e+grd!
    integer :: eeq = 4
    integer :: ceh = 5
    integer :: param = 6
    integer :: gxtb = 7
  end type enum_tblite_method
  type(enum_tblite_method),parameter,public :: xtblvl = enum_tblite_method()

#ifdef WITH_GXTB
  logical,parameter,public :: have_gxtb = .true.
#else
  logical,parameter,public :: have_gxtb = .false.
#endif

  !> Conversion factor from Kelvin to Hartree
  real(wp),parameter :: ktoau = 3.166808578545117e-06_wp

  public :: wavefunction_type,xtb_calculator
  public :: tblite_ctx,tblite_resultstype
  public :: tblite_setup,tblite_singlepoint,tblite_addsettings
  public :: tblite_getwbos
  public :: tblite_add_solv
  public :: tblite_add_efield
  public :: tblite_getcharges
  public :: tblite_getdipole
  public :: tblite_quick_ceh_q

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine tblite_setup(mol,chrg,uhf,lvl,etemp,tblite,ceh_guess)
!*****************************************************************
!* subroutine tblite_setup initializes the tblite object which is
!* passed between the CREST calculators and this module
!*****************************************************************
#ifdef WITH_TBLITE
#ifdef WITH_GXTB
    use multicharge,only:get_charges
#endif
#endif
    implicit none
    type(coord),intent(in)  :: mol
    integer,intent(in)      :: chrg
    integer,intent(in)      :: uhf
    type(tblite_data),intent(inout) :: tblite
    integer,intent(in)      :: lvl
    real(wp),intent(in)     :: etemp
    logical,intent(in),optional :: ceh_guess
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(error_type),allocatable :: error
    type(param_record) :: param

    real(wp) :: etemp_au,energy
    real(wp),allocatable :: grad(:,:)
    logical :: pr
    integer :: io,nspin

    pr = (tblite%ctx%verbosity > 0)

!>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

!>--- select parametrization and set up calculator
    tblite%lvl = lvl
    select case (tblite%lvl)
    case (xtblvl%gfn1)
      if (pr) call tblite%ctx%message("tblite> Setting up GFN1-xTB calculation")
      call new_gfn1_calculator(tblite%calc,mctcmol,error)
    case (xtblvl%gfn2)
      if (pr) call tblite%ctx%message("tblite> Setting up GFN2-xTB calculation")
      call new_gfn2_calculator(tblite%calc,mctcmol,error)
    case (xtblvl%ipea1)
      if (pr) call tblite%ctx%message("tblite> Setting up IPEA1-xTB calculation")
      call new_ipea1_calculator(tblite%calc,mctcmol,error)
    case (xtblvl%ceh)
      if (pr) call tblite%ctx%message("tblite> Setting up CEH calculation")
      call new_ceh_calculator(tblite%calc,mctcmol,error)
    case (xtblvl%eeq)
      if (pr) call tblite%ctx%message("tblite> Setting up D4 EEQ charges calculation")
      call new_ceh_calculator(tblite%calc,mctcmol,error) !> doesn't matter but needs initialization
    case (xtblvl%param)
      if (pr) call tblite%ctx%message("tblite> Setting up xtb calculator from parameter file")
      if (allocated(tblite%paramfile)) then
        call tblite_read_param_record(tblite%paramfile,param,io)
        call new_xtb_calculator(tblite%calc,mctcmol,param,error)
        if (allocated(error)) then
          write (stdout,*) 'Could not read tblite parameter file '//tblite%paramfile
          error stop
        end if
      else
        if (pr) call tblite%ctx%message("tblite> parameter file does not exist, defaulting to GFN2-xTB")
        call new_gfn2_calculator(tblite%calc,mctcmol,error)
      end if
#ifdef WITH_GXTB
    case (xtblvl%gxtb)
      if (pr) call tblite%ctx%message("tblite> Setting up g-xTB calculation")
      call new_gxtb_calculator(tblite%calc,mctcmol,error)
#else
    case (xtblvl%gxtb)
      write (stdout,'(a)') 'Error: g-xTB via tblite not available (compiled without WITH_GXTB).'
      write (stdout,'(a)') 'This code path should not be reached — use the xtb binary fallback.'
      error stop
#endif
    case default
      call tblite%ctx%message("Error: Unknown method in tblite!")
      error stop
    end select
    if (pr) call tblite%ctx%message('')

!>-- setup wavefunction object
    etemp_au = etemp*ktoau
    nspin = merge(2,1,tblite%spin_polarized)
    call new_wavefunction(tblite%wfn,mol%nat,tblite%calc%bas%nsh,  &
    &              tblite%calc%bas%nao,nspin,etemp_au)
#ifdef WITH_GXTB
    if (tblite%lvl == xtblvl%gxtb) then
      call sad_guess(mctcmol,tblite%calc,tblite%wfn)
    end if
#endif
    if (ceh_guess) then
      call tblite_internal_ceh_guess(mctcmol,tblite)
    end if

!>--- spin-polarization setup (spGFN2-xTB etc.)
    if (tblite%spin_polarized) then
      block
        class(container_type),allocatable :: cont
        type(spin_polarization),allocatable :: spin
        real(wp),allocatable :: wll(:,:,:)
        allocate (spin)
        call get_spin_constants(wll,mctcmol,tblite%calc%bas)
        call new_spin_polarization(spin,mctcmol,wll,tblite%calc%bas%nsh_id)
        call move_alloc(spin,cont)
        call tblite%calc%push_back(cont)
      end block
    end if

!>--- for methods with an auxiliary charge model (e.g., gxTB), pre-allocate wfn_aux.
!>--- Charges are updated at each singlepoint call (geometry-dependent).
#ifdef WITH_GXTB
    if (allocated(tblite%calc%charge_model)) then
      if (allocated(tblite%wfn_aux)) deallocate (tblite%wfn_aux)
      allocate (tblite%wfn_aux)
      call new_wavefunction(tblite%wfn_aux,mctcmol%nat,tblite%calc%bas%nsh,0,1,0.0_wp,.true.)
    end if
#endif

#else /* WITH_TBLITE */
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_setup

!========================================================================================!

  subroutine tblite_add_solv(mol,chrg,uhf,tblite,smodel,solvent)
!***********************************************************************
!* This subroutine sets up the implicit solvation container for tblite
!***********************************************************************
#ifdef WITH_TBLITE
    use tblite_container,only:container_type
    use tblite_solvation,only:new_solvation,tblite_solvation_type => solvation_type, &
    &                         solvent_data,get_solvent_data,solvation_input,  &
    &                         ddx_input,ddx_solvation_model,alpb_input,alpb_solvation, &
    &                         cds_input,new_solvation_cds,shift_input,new_solvation_shift
    use tblite_features,only:get_tblite_feature
#endif
    implicit none
    type(coord),intent(in)  :: mol
    integer,intent(in)      :: chrg
    integer,intent(in)      :: uhf
    type(tblite_data),intent(inout) :: tblite
    character(len=:),allocatable,intent(in) :: smodel
    character(len=:),allocatable,intent(in) :: solvent
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(error_type),allocatable :: error

    class(container_type),allocatable :: cont
    class(tblite_solvation_type),allocatable :: solv
    type(solvation_input),allocatable :: solv_inp
    type(solvent_data) :: solv_data
    type(alpb_input)  :: alpb_tmp
    type(cds_input)   :: cds_tmp
    type(shift_input) :: shift_tmp
    character(len=:),allocatable :: str,solvdum,method
    logical :: pr

    if (.not.allocated(smodel).or..not.allocated(solvent)) then
      return
    end if
    pr = (tblite%ctx%verbosity > 0)

!>--- some tblite calculators have nothing to do with implicit solvation
    if (tblite%lvl > 3.and.tblite%lvl .ne. xtblvl%param) then
      if (pr) call tblite%ctx%message("tblite> skipping implicit solvation setup for this potential")
      return
    end if
    select case (tblite%lvl)
    case (xtblvl%gfn1)
      method = 'gfn1'
    case (xtblvl%gfn2)
      method = 'gfn2'
    end select

!>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

    if (pr) call tblite%ctx%message("tblite> setting up tblite implicit solvation")
!>--- generat solvation parametrization
    if (solvent == 'h2o') then !> special case: tblite doesn't know 'h2o', only 'water' ...
      solvdum = 'water'
    else
      solvdum = solvent
    end if
    solv_data = get_solvent_data(solvdum)
    if (solv_data%eps <= 0.0_wp) then
      write (stdout,'(a)') 'Error: unknown solvent "'//solvdum//'" for tblite implicit solvation!'
      error stop
    end if
    allocate (solv_inp)
    select case (trim(smodel))
    case ('gbsa')
      if (pr) call tblite%ctx%message("tblite> using GBSA/"//solvdum)
      alpb_tmp%dielectric_const = solv_data%eps
      alpb_tmp%alpb = .false.
      !alpb_tmp%method=method
      alpb_tmp%solvent = solv_data%solvent
      !alpb_tmp%xtb=.true.
      allocate (solv_inp%alpb,source=alpb_tmp)
      cds_tmp%alpb = .false.
      cds_tmp%solvent = solv_data%solvent
      !cds_tmp%method=method
      allocate (solv_inp%cds,source=cds_tmp)
      shift_tmp%alpb = .false.
      shift_tmp%solvent = solv_data%solvent
      !shift_tmp%method=method
      allocate (solv_inp%shift,source=shift_tmp)
    case ('cpcm','cosmo','pcm')
      !> Continuum solvation models are provided via the ddX library since
      !> tblite 0.6.x and require tblite compiled with ddX support (ddx=true)
      if (.not.get_tblite_feature('ddx')) then
        write (stdout,'(a)') 'Error: "'//trim(smodel)//'" solvation requires '// &
        &                    'tblite compiled with ddX support!'
        error stop
      end if
      if (pr) call tblite%ctx%message("tblite> using "//trim(smodel)//"/"//solvdum)
      allocate (solv_inp%ddx)
      select case (trim(smodel))
      case ('cosmo')
        solv_inp%ddx = ddx_input(ddx_solvation_model%cosmo,solv_data%eps)
      case ('pcm')
        solv_inp%ddx = ddx_input(ddx_solvation_model%pcm,solv_data%eps)
      case default !> cpcm
        solv_inp%ddx = ddx_input(ddx_solvation_model%cpcm,solv_data%eps)
      end select
    case ('alpb')
      if (pr) call tblite%ctx%message("tblite> using ALPB/"//solvdum)
      alpb_tmp%dielectric_const = solv_data%eps
      alpb_tmp%alpb = .true.
      !alpb_tmp%method=method
      alpb_tmp%solvent = solv_data%solvent
      !alpb_tmp%xtb=.true.
      allocate (solv_inp%alpb,source=alpb_tmp)
      cds_tmp%alpb = .true.
      cds_tmp%solvent = solv_data%solvent
      !cds_tmp%method=method
      allocate (solv_inp%cds,source=cds_tmp)
      shift_tmp%alpb = .true.
      shift_tmp%solvent = solv_data%solvent
      !shift_tmp%method=method
      allocate (solv_inp%shift,source=shift_tmp)
    case default
      write (stdout,'(a)') 'Error: unknown tblite implicit solvation model "'//trim(smodel)//'"!'
      error stop
    end select

!>--- add electrostatic (Born part) to calculator
    call new_solvation(solv,mctcmol,solv_inp,error,method)
    if (allocated(error)) then
      write (stdout,'(a)') 'Error: failed to set up tblite implicit solvation: '//error%message
      error stop
    end if
    call move_alloc(solv,cont)
    call tblite%calc%push_back(cont)
!>--- add hbond and dispersion part to calculator
    if (allocated(solv_inp%cds)) then
      block
        class(tblite_solvation_type),allocatable :: cds
        call new_solvation_cds(cds,mctcmol,solv_inp,error,method)
        if (allocated(error)) return
        call move_alloc(cds,cont)
        call tblite%calc%push_back(cont)
      end block
    end if
!>--- add gsolv shift to calculator
    if (allocated(solv_inp%shift)) then
      block
        class(tblite_solvation_type),allocatable :: shift
        call new_solvation_shift(shift,solv_inp,error,method)
        if (allocated(error)) return
        call move_alloc(shift,cont)
        call tblite%calc%push_back(cont)
      end block
    end if

    deallocate (solv_inp)

#else /* WITH_TBLITE */
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_add_solv

!========================================================================================!

  subroutine tblite_singlepoint(mol,chrg,uhf,tblite,energy,gradient,iostatus)
!**************************************************
!* The actual calculator call.
!* The tblite object must be set up at this point
!**************************************************
#ifdef WITH_TBLITE
#ifdef WITH_GXTB
    use multicharge,only:get_charges
#endif
#endif
    implicit none
    type(coord),intent(in)   :: mol
    integer,intent(in)       :: chrg
    integer,intent(in)       :: uhf
    type(tblite_data),intent(inout) :: tblite
    real(wp),intent(out)     :: energy
    real(wp),intent(out)     :: gradient(3,mol%nat)
    integer,intent(out)      :: iostatus
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    type(error_type),allocatable :: error
    real(wp) :: sigma(3,3)
    logical :: pr
    integer :: verbosity

    iostatus = 0
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    pr = (tblite%ctx%verbosity > 0)
    if (tblite%ctx%verbosity > 1) then
      verbosity = tblite%ctx%verbosity
    else
      verbosity = 0
    end if

!>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)

!>--- update geometry-dependent EEQ-BC charges in wfn_aux (allocated once in tblite_setup)
#ifdef WITH_GXTB
    if (allocated(tblite%wfn_aux)) then
      call get_charges(tblite%calc%charge_model,mctcmol,error,tblite%wfn_aux%qat(:,1), &
      &                dqdr=tblite%wfn_aux%dqatdr(:,:,:,1),dqdL=tblite%wfn_aux%dqatdL(:,:,:,1))
      if (allocated(error)) then
        if (pr) call tblite%ctx%message("tblite> auxiliary charge model failed: "//error%message)
        iostatus = 1
        return
      end if
    end if
#endif

!>--- call the singlepoint routine
    select case (tblite%lvl)
    case default
#ifdef WITH_GXTB
      if (allocated(tblite%wfn_aux)) then
        call xtb_singlepoint(tblite%ctx,mctcmol,tblite%calc,tblite%wfn,tblite%accuracy, &
        &                    energy,gradient,sigma,verbosity,results=tblite%res,wfn_aux=tblite%wfn_aux)
      else
        call xtb_singlepoint(tblite%ctx,mctcmol,tblite%calc,tblite%wfn,tblite%accuracy, &
        &                    energy,gradient,sigma,verbosity,results=tblite%res)
      end if
#else
      call xtb_singlepoint(tblite%ctx,mctcmol,tblite%calc,tblite%wfn,tblite%accuracy, &
      &                    energy,gradient,sigma,verbosity,results=tblite%res)
#endif
    case (xtblvl%ceh)
      call ceh_singlepoint(tblite%ctx,tblite%calc,mctcmol,tblite%wfn, &
      &              tblite%accuracy,verbosity)
    case (xtblvl%eeq)
      call eeq_guess(mctcmol,tblite%calc,tblite%wfn,error)
    end select

    if (tblite%ctx%failed()) then
      !> Tear down the error stack to send the actual error messages back
      if (pr) call tblite%ctx%message("tblite> Singlepoint calculation failed")
      iostatus = 1
    end if

#else /* WITH_TBLITE */
    iostatus = 0
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_singlepoint

!========================================================================================!

#ifdef WITH_TBLITE
  subroutine tblite_mol2mol(mol,chrg,uhf,mctcmol)
!*************************************************************************
!* tblite uses its own molecule type thats different from our coord type
!* This routine does the minimal conversion
!*************************************************************************
    implicit none
    !> input & output
    type(coord) :: mol
    integer,intent(in) :: chrg
    integer,intent(in) :: uhf
    type(structure_type),intent(out) :: mctcmol
    !> locals
    real(wp) :: fchrg

    fchrg = real(chrg,wp)

    !>--- make an mctcmol object from mol
    if (.not.allocated(mol%lat)) then
      call new(mctcmol,mol%at,mol%xyz,charge=fchrg,uhf=uhf)
    else
      call new(mctcmol,mol%at,mol%xyz,charge=fchrg,uhf=uhf,lattice=mol%lat)
    end if

  end subroutine tblite_mol2mol
#endif

!========================================================================================!

  subroutine tblite_addsettings(tblite,maxscc,rdwbo,saveint,accuracy)
!**********************************************************
!* tblite_addsettings is used to add other settings from
!* CRESTs calculation object to the xtb_calculator
!**********************************************************
    implicit none
    type(tblite_data),intent(inout) :: tblite
    integer,intent(in) :: maxscc
    logical,intent(in) :: rdwbo
    logical,intent(in) :: saveint
    real(wp),intent(in) :: accuracy
#ifdef WITH_TBLITE
    tblite%calc%max_iter = maxscc
    tblite%calc%save_integrals = (rdwbo.or.saveint)
    tblite%accuracy = accuracy
#endif
  end subroutine tblite_addsettings

  subroutine tblite_add_efield(tblite,efield)
!**********************************************************
!* tblite_add_efield
!* if efield is allocated, add it to the tblite calculator
!**********************************************************
#ifdef WITH_TBLITE
    use tblite_container,only:container_type
    use tblite_external_field,only:electric_field
#endif
    implicit none
    type(tblite_data),intent(inout) :: tblite
    real(wp),intent(in),allocatable :: efield(:)
    class(container_type),allocatable :: cont
    logical :: pr
    character(len=90) :: str
#ifdef WITH_TBLITE
    pr = (tblite%ctx%verbosity > 0)
    if (allocated(efield)) then
      if (pr) then
        write (str,'(a,3(es10.3),a)') "tblite> Calculation includes the following electric field:"
        call tblite%ctx%message(trim(str))
        write (str,'(8x, a,3(es15.5,1x),a)') "[",efield,"] V/Å"
        call tblite%ctx%message(trim(str))
        call tblite%ctx%message('')
      end if
      cont = electric_field(efield*vatoau)
      call tblite%calc%push_back(cont)
    end if
#endif
  end subroutine tblite_add_efield

!========================================================================================!

  subroutine tblite_getwbos(tblite,nat,wbo)
!***************************
!* obtain wbos from tblite
!***************************
    implicit none
    type(tblite_data),intent(inout) :: tblite
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat)
    real(wp),allocatable :: S(:,:)
    integer :: nao,i
    real(wp),allocatable :: focca(:),foccb(:)
    real(wp),allocatable :: Pa(:,:),Pb(:,:)
    wbo = 0.0_wp
#ifdef WITH_TBLITE
    select case (tblite%lvl)
    case default

      nao = tblite%calc%bas%nao
      allocate (Pa(nao,nao),Pb(nao,nao))
      call split_foccab(nao,tblite%wfn%focc,tblite%wfn%nel(1),tblite%wfn%nel(2), &
      & focca,foccb)
      call density_matrix(nao,focca,tblite%wfn%coeff(:,:,1),Pa)
      call density_matrix(nao,foccb,tblite%wfn%coeff(:,:,1),Pb)
      call get_wbo(nat,nao,Pa,Pb,tblite%res%overlap,tblite%calc%bas%ao2at,wbo)

    case (xtblvl%ceh)
      !> no external access to the overlap in CEH, hence use the Wiberg BO with S=I
      nao = tblite%calc%bas%nao
      allocate (S(nao,nao),source=0.0_wp)
      do i = 1,nao
        S(i,i) = 1.0_wp
      end do
      call get_wbo_rhf(nat,tblite%calc%bas%nao,tblite%wfn%density, &
      &                S,tblite%calc%bas%ao2at,wbo)
      wbo = wbo*2.0_wp !> somehow this is much better

    case (xtblvl%eeq)
      wbo = 0.0_wp
    end select
#endif
  end subroutine tblite_getwbos

!========================================================================================!

  subroutine tblite_getcharges(mol,tblite,q)
!**************************************
!* obtain molecular dipole from tblite
!**************************************
    implicit none
    type(coord) :: mol
    type(tblite_data),intent(inout) :: tblite
    real(wp),intent(out) :: q(mol%nat)
#ifdef WITH_TBLITE
    q = 0.0_wp
    q(:) = tblite%wfn%qat(:,1)
#else
    q = 0.0_wp
#endif
  end subroutine tblite_getcharges

!========================================================================================!

  subroutine tblite_getdipole(mol,chrg,uhf,tblite,dipole)
!**************************************
!* obtain molecular dipole from tblite
!**************************************
    implicit none
    type(coord) :: mol
    integer :: chrg,uhf
    type(tblite_data),intent(inout) :: tblite
    real(wp),intent(out) :: dipole(3)
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    dipole = 0.0_wp
    call tblite_mol2mol(mol,chrg,uhf,mctcmol)
    !> dipole moment is obtained from molecular charges and atomic dipole moments
    call get_molecular_dipole_moment(mctcmol,tblite%wfn%qat(:,1), &
    &    tblite%wfn%dpat(:,:,1),dipole)
#else
    dipole = 0.0_wp
#endif
  end subroutine tblite_getdipole

!========================================================================================!

#ifdef WITH_TBLITE
  subroutine tblite_read_param_record(paramfile,param,io)
    use tomlf
    implicit none
    character(len=*),intent(in) :: paramfile
    type(param_record),intent(out) :: param
    integer,intent(out) :: io
    type(error_type),allocatable :: error
    type(toml_table),allocatable :: table
    type(toml_error),allocatable :: terror
    type(toml_context) :: tcontext
    logical,parameter :: color = .true.

    io = 1

    call toml_load(table,paramfile,error=terror,context=tcontext, &
    & config=toml_parser_config(color=color))
    if (allocated(terror)) then
      io = 1
      return
    end if

    call param%load_from_toml(table,error)

    if (allocated(error)) then
      io = 1
    else
      io = 0
    end if
    if (allocated(table)) deallocate (table)

  end subroutine tblite_read_param_record
#endif

!========================================================================================!

#ifdef WITH_TBLITE
  subroutine tblite_internal_ceh_guess(mctcmol,tblite)
    !*********************************************************
    !* Init the tblite calculator with a set of CEH charges
    !*********************************************************
    implicit none
    type(tblite_data),intent(inout) :: tblite
    type(structure_type),intent(in) :: mctcmol
    !> LOCAL
    type(wavefunction_type) :: wfn_ceh
    type(xtb_calculator)    :: calc_ceh
    type(error_type),allocatable :: error
    integer :: verbosity
    logical :: pr
    real(wp),parameter :: etemp_guess_au = 4000.0_wp*ktoau

    !> if we only do a eeq or ceh calc, we don't need this, so return
    select case (tblite%lvl)
    case default
      continue
    case (xtblvl%ceh,xtblvl%eeq)
      return
    end select

    pr = (tblite%ctx%verbosity > 0)
    if (tblite%ctx%verbosity > 1) then
      verbosity = tblite%ctx%verbosity
    else
      verbosity = 0
    end if

    !> ceh guess calculator and wavefunction
    call new_ceh_calculator(calc_ceh,mctcmol,error)
    if (allocated(error)) return
    call new_wavefunction(wfn_ceh,mctcmol%nat,calc_ceh%bas%nsh, &
    &                     calc_ceh%bas%nao,1,etemp_guess_au)

    !> TODO ceh guess efield

    call ceh_singlepoint(tblite%ctx,calc_ceh,mctcmol,wfn_ceh, &
    &              tblite%accuracy,verbosity)

    if (tblite%ctx%failed()) then
      if (pr) then
        call tblite%ctx%get_error(error)
        call tblite%ctx%message("CEH singlepoint calculation failed")
        call tblite%ctx%message("-> "//error%message)
      end if
      return
    end if

    !> pass on to actual calculator
    tblite%wfn%qat(:,1) = wfn_ceh%qat(:,1)
    call shell_partition(mctcmol,tblite%calc,tblite%wfn)

  end subroutine tblite_internal_ceh_guess
#endif

!========================================================================================!

  subroutine tblite_quick_ceh_q(mol,q,chrg,uhf,pr,prch)
    !*********************************************************
    !* Calculate CEH charges
    !*********************************************************
    implicit none
    type(coord),intent(in) :: mol
    integer,intent(in) :: chrg
    real(wp),intent(out),allocatable :: q(:)
    integer,intent(in),optional :: uhf
    logical,intent(in),optional :: pr
    integer,intent(in),optional :: prch
#ifdef WITH_TBLITE
    type(structure_type) :: mctcmol
    !> LOCAL
    type(wavefunction_type) :: wfn_ceh
    type(xtb_calculator)    :: calc_ceh
    type(tblite_ctx)        :: ctx
    type(error_type),allocatable :: error
#endif
    integer :: verbosity,uhf_loc
    logical :: pr_loc
    real(wp),parameter :: etemp_guess_au = 4000.0_wp*ktoau
    real(wp),parameter :: accuracy = 1.0_wp

    pr_loc = .false.
    if (present(pr)) pr_loc = pr
    verbosity = 0
    if (pr_loc) verbosity = 2

    allocate (q(mol%nat),source=0.0_wp)

#ifdef WITH_TBLITE
    uhf_loc = 0
    if (present(uhf)) uhf_loc = uhf
    if (present(prch)) ctx%unit = prch

    !>--- make an mctcmol object from mol
    call tblite_mol2mol(mol,chrg,uhf_loc,mctcmol)

    !> ceh guess calculator and wavefunction
    call new_ceh_calculator(calc_ceh,mctcmol,error)
    if (allocated(error)) return
    call new_wavefunction(wfn_ceh,mctcmol%nat,calc_ceh%bas%nsh, &
    &                     calc_ceh%bas%nao,1,etemp_guess_au)

    !> TODO ceh guess efield

    call ceh_singlepoint(ctx,calc_ceh,mctcmol,wfn_ceh, &
    &              accuracy,verbosity)

    if (ctx%failed()) then
      if (pr_loc) then
        call ctx%get_error(error)
        call ctx%message("CEH singlepoint calculation failed")
        call ctx%message("-> "//error%message)
      end if
      return
    end if

    !> pass on the charges
    q(:) = wfn_ceh%qat(:,1)
#else /* WITH_TBLITE */
    write (stdout,*) 'Error: Compiled without tblite support!'
    write (stdout,*) 'Use -DWITH_TBLITE=true in the setup to enable this function'
    error stop
#endif
  end subroutine tblite_quick_ceh_q

! ══════════════════════════════════════════════════════════════════════════════
#ifdef WITH_TBLITE
  subroutine get_spin_constants(wll,mol,bas)
    use tblite_basis_type,only:basis_type
    use tblite_data_spin,only:get_spin_constant
    real(wp),allocatable,intent(out) :: wll(:,:,:)
    type(structure_type),intent(in) :: mol
    type(basis_type),intent(in) :: bas

    integer :: izp,ish,jsh,il,jl

    allocate (wll(bas%nsh,bas%nsh,mol%nid),source=0.0_wp)

    do izp = 1,mol%nid
      do ish = 1,bas%nsh_id(izp)
        il = bas%cgto(ish,izp)%ang
        do jsh = 1,bas%nsh_id(izp)
          jl = bas%cgto(jsh,izp)%ang
          wll(jsh,ish,izp) = get_spin_constant(jl,il,mol%num(izp))
        end do
      end do
    end do
  end subroutine get_spin_constants
#endif

!========================================================================================!
!========================================================================================!
end module tblite_api
