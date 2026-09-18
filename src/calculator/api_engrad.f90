!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2025 Philipp Pracht
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

!> module api_engrad
!> a collection of engrad calls for different APIs
!> this builds the communication between CRESTs
!> "calculation_settings" and the respective API setups

module api_engrad

  use iso_fortran_env,only:wp => real64,stdout => output_unit
  use strucrd
  use calc_type
  use iomod,only:makedir,directory_exist,remove,dump_array_to_tmp
  use omp_lib
  !> API modules
  use api_helpers
  use tblite_api
  use gfn0_api
  use gfnff_api
  use crest_solvation,only:solvation_setup,solvation_core
  use crest_electrostatic,only:electrostatic_core
  use libpvol_api
  use lj
  use approxg_module
  use penalty_module
  use mlip_sc
  implicit none
!>--- private module variables and parameters
  private

  public :: tblite_engrad
  public :: gfn0_engrad,gfn0occ_engrad
  public :: gfnff_engrad
  public :: libpvol_engrad
  public :: lj_engrad   !> RE-EXPORT
  public :: modelhessian_engrad
  public :: rmsd_engrad
  public :: mlip_engrad
  public :: preinit_mlip_parallel
  public :: solvation_engrad
  public :: electrostatic_engrad

!=========================================================================================!
!=========================================================================================!
contains    !> MODULE PROCEDURES START HERE
!=========================================================================================!
!=========================================================================================!

  subroutine tblite_engrad(mol,calc,energy,grad,iostatus)
!******************************************************
!* Interface singlepoint call between CREST and tblite
!******************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    logical :: loadnew,pr

    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.

!>--- setup system call information
    !$omp critical
    call tblite_init(calc,loadnew)
!>--- tblite printout handling
    call api_handle_output(calc,'tblite.out',mol,pr)
    if (pr.or.calc%prstdout) then
      !> tblite uses its context (ctx) type, rather than calc%prch
      calc%tblite%ctx%unit = calc%prch
      calc%tblite%ctx%verbosity = 1
      if (calc%prstdout) then
        !> special case, fwd to stdout (be carefule with this!)
        calc%tblite%ctx%unit = stdout
        calc%tblite%ctx%verbosity = 2
      end if
    else
      calc%tblite%ctx%verbosity = 0
    end if

!>-- populate parameters and wavefunction
    if (loadnew) then
      call tblite_setup(mol,calc%chrg,calc%uhf,calc%tblitelvl,calc%etemp,calc%tblite,calc%ceh_guess)

      call tblite_addsettings(calc%tblite,calc%maxscc,calc%rdwbo,calc%saveint,calc%accuracy)

      call tblite_add_efield(calc%tblite,calc%efield)

      call tblite_add_solv(mol,calc%chrg,calc%uhf,calc%tblite, &
      &    calc%solvmodel,calc%solvent)
    end if
    !$omp end critical

!>--- do the engrad call
    call initsignal()
    call tblite_singlepoint(mol,calc%chrg,calc%uhf,calc%tblite, &
    &                       energy,grad,iostatus)
    if (iostatus /= 0) return
    if (.not.calc%prstdout) &
    & call api_print_e_grd(pr,calc%prch,mol,energy,grad)

!>--- postprocessing, getting other data
    !$omp critical
    call tblite_properties(calc,mol,iostatus)
    !$omp end critical

    return
  end subroutine tblite_engrad

!========================================================================================!

  subroutine solvation_engrad(mol,calc,energy,grad,iostatus)
!*********************************************************************
!* Interface singlepoint call for the composite solvation calculator
!*********************************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    iostatus = 0
    if (.not.allocated(calc%solv)) then
      iostatus = 1
      return
    end if

!>--- lazily resolve dielectric constant + GFN2/ALPB parameters
    call solvation_setup(mol,calc%solv,iostatus)
    if (iostatus /= 0) return

!>--- stitched energy + gradient
    call solvation_core(mol,calc%chrg,calc%solv,energy,grad,iostatus)
    if (.not.calc%prstdout) call api_print_e_grd(.false.,calc%prch,mol,energy,grad)

    return
  end subroutine solvation_engrad

!========================================================================================!

  subroutine electrostatic_engrad(mol,calc,energy,grad,iostatus)
!*********************************************************************
!* Interface singlepoint call for the charge-equilibration calculator
!*********************************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus
    character(len=:),allocatable :: model

    iostatus = 0
    if (.not.allocated(calc%eeq)) then
      iostatus = 1
      return
    end if
    model = 'eeqbc'
    if (allocated(calc%eeq%charge_model)) model = calc%eeq%charge_model

!>--- energy, gradient and atomic charges (kept on the settings)
    if (.not.allocated(calc%qat)) allocate (calc%qat(mol%nat))
    call electrostatic_core(mol,calc%chrg,model,energy,grad,calc%qat,iostatus)
    if (iostatus /= 0) return
    if (.not.calc%prstdout) call api_print_e_grd(.false.,calc%prch,mol,energy,grad)

    return
  end subroutine electrostatic_engrad

!========================================================================================!

  subroutine gfn0_engrad(mol,calc,g0calc,energy,grad,iostatus)
!************************************************
!* Interface singlepoint call between CREST and
!* the GFN0 engrad standard implementation
!************************************************
    implicit none
    !> INPUT
    type(coord) :: mol
    type(calculation_settings) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    !> OUTPUT
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus
    !> LOCAL
    type(gfn0_results) :: res
    logical :: loadnew
    logical :: pr

    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information
    !$omp critical
    call gfn0_init(calc,g0calc,loadnew)
!>--- printout handling
    call api_handle_output(calc,'gfn0.out',mol,pr)
!>-- populate parameters and wavefunction
    if (loadnew) then
      call gfn0_setup(mol,calc%chrg,calc%uhf,g0calc)
      call gfn0_init2(mol,calc,g0calc)
    end if
    call gfn0_init3(mol,calc,g0calc)
    !$omp end critical

!>--- do the engrad call
    call initsignal()
    call gfn0_sp(mol,calc%chrg,calc%uhf,g0calc,energy,grad,iostatus,res)
    if (iostatus /= 0) return
    if (pr) then
      call gfn0_print(calc%prch,g0calc,res)
      if (.not.calc%prstdout) &
      & call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data
    !$omp critical
    call gfn0_properties(calc,calc%g0calc,mol,iostatus)
    !$omp end critical

    return
  end subroutine gfn0_engrad

!========================================================================================!

  subroutine gfn0occ_engrad(mol,calc,g0calc,energy,grad,iostatus)
!************************************************
!* Interface singlepoint call between CREST and
!* the GFN0 multi-occupation implementation
!************************************************
    implicit none
    !> INPUT
    type(coord) :: mol
    type(calculation_settings) :: calc
    type(gfn0_data),intent(inout),allocatable  :: g0calc
    !> OUTPUT
    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(:,:)
    integer,intent(out) :: iostatus
    !> LOCAL
    type(gfn0_results) :: res
    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information
    !$omp critical
    call gfn0occ_init(calc,g0calc,loadnew)
!>--- printout handling
    call api_handle_output(calc,'gfn0.out',mol,pr)
!>--- populate parameters and wavefunction
    if (loadnew) then
      call gfn0_setup(mol,calc%chrg,calc%uhf,g0calc)
      call gfn0occ_init2(mol,calc,g0calc)
    end if
    call gfn0occ_init3(mol,calc,g0calc)
    !$omp end critical

!>--- do the engrad call
    call initsignal()
    call gfn0_sp_occ(mol,calc%chrg,calc%uhf,calc%occ,g0calc, &
    &    energy,grad,iostatus,res)
    if (iostatus /= 0) return
    if (pr) then
      call gfn0_print(calc%prch,g0calc,res)
      if (.not.calc%prstdout) &
      & call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data
    !$omp critical
    call gfn0_properties(calc,g0calc,mol,iostatus)
    !$omp end critical

    return
  end subroutine gfn0occ_engrad

!========================================================================================!

  subroutine gfnff_engrad(mol,calc,energy,grad,iostatus)
!******************************************************************
!* Interface singlepoint call between CREST and GFN-FF force field
!******************************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    real(wp),allocatable :: q(:)
    iostatus = 0
    pr = .false.
!>--- setup calculation data
    !$omp critical
    call gfnff_init(calc,loadnew)
!>--- printout handling
    call api_handle_output(calc,'gfnff.out',mol,pr)

!>--- populate parameters and neighbourlists
    if (loadnew) then
      if (calc%ceh_guess) then
        if (pr) then
          write (calc%prch,'(/,a)') 'Initializing atomic reference charges from CEH model'
        end if
        !> Compute CEH atomic charges and pass them to GFN-FF (in-memory) as the
        !> reference topology charges (qa) used in the force-field construction.
        call tblite_quick_ceh_q(mol,q,calc%chrg,pr=pr,prch=calc%prch)
        call gfnff_set_refcharges(calc%ff_dat,q)
        deallocate (q)
      end if

      !> hand any user-defined fragmentation to GFN-FF before topology setup
      if (allocated(calc%gff_fragments)) then
        call gfnff_set_fragments(calc%ff_dat,mol%nat,calc%gff_fragments)
      end if

      call gfnff_api_setup(mol,calc%chrg,calc%ff_dat,iostatus,pr,calc%prch)
    end if
    !$omp end critical
    if (iostatus /= 0) return

!>--- do the engrad call
    call initsignal()
    call gfnff_sp(mol,calc%ff_dat,energy,grad,iostatus)
    if (iostatus /= 0) return

!>--- printout
    if (pr) then
      call gfnff_printout(calc%prch,calc%ff_dat)
      if (.not.calc%prstdout) &
      & call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data
    !$omp critical
    call gfnff_properties(calc,mol,iostatus)
    !$omp end critical

    return
  end subroutine gfnff_engrad

!========================================================================================!

  subroutine libpvol_engrad(mol,calc,energy,grad,iostatus)
!***************************************************************
!* Interface singlepoint call between CREST and XHC force field
!***************************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information
    !$omp critical
    call libpvol_initcheck(calc,loadnew)
!>--- printout handling
    call api_handle_output(calc,'libpvol.out',mol,pr)
!>--- populate parameters
    if (loadnew) then
      !> call libpvol with verbosity turned off
      call libpvol_setup(mol,calc%libpvol,calc%extpressure,calc%pvmodel, &
      &        calc%ngrid,calc%proberad,calc%vdwset,calc%pvradscal,pr,calc%prch,iostatus)
    end if
    !$omp end critical
    if (iostatus /= 0) return

!>--- do the engrad call
    call initsignal()
    call libpvol_sp(mol,calc%libpvol,energy,grad,iostatus)
    if (iostatus /= 0) return

!>--- printout
    if (pr) then
      !> the libpvol_sp call includes the printout within libpvol-lib
      if (.not.calc%prstdout) &
      & call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data

    return
  end subroutine libpvol_engrad

!========================================================================================!

  subroutine modelhessian_engrad(mol,calc,energy,grad,iostatus)
!***************************************************************
!* Interface singlepoint call between CREST and XHC force field
!***************************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings) :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io,n3
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information
    !$omp critical
!>--- printout handling
    call api_handle_output(calc,'modh.out',mol,pr)
!>--- populate parameters
    n3 = mol%nat*3
    if (calc%ag%dim .ne. mol%nat) then
      calc%ag%pr = calc%prstdout.and..not.calc%numgrad

      if (allocated(calc%ag%hess)) deallocate (calc%ag%hess)
      allocate (calc%ag%hess(n3,n3),source=0.0_wp)

      if (allocated(calc%ag%h)) deallocate (calc%ag%h)
      allocate (calc%ag%h(n3*(n3+1)/2),source=0.0_wp)

      if (allocated(calc%ag%freq)) deallocate (calc%ag%freq)
      allocate (calc%ag%freq(n3),source=0.0_wp)

      if (allocated(calc%ag%xyz)) deallocate (calc%ag%xyz)
      allocate (calc%ag%xyz(3,mol%nat),source=0.0_wp)

      calc%ag%dim = mol%nat
    else
      calc%ag%hess(:,:) = 0.0_wp
      calc%ag%h(:) = 0.0_wp
    end if
    !$omp end critical
    if (iostatus /= 0) return

!>--- do the engrad call
    call initsignal()
    call modh_engrad(mol,calc%ag,energy,grad,iostatus)
    if (iostatus /= 0) return

!>--- printout
    if (pr) then
      !> the libpvol_sp call includes the printout within libpvol-lib
      if (.not.calc%prstdout) &
      & call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data

    return
  end subroutine modelhessian_engrad

!========================================================================================!

  subroutine rmsd_engrad(mol,calc,energy,grad,iostatus)
!**************************************************************************
!* Interface singlepoint to add RMSD penalty function (as in metadynamics)
!**************************************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings),target :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io,nall
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- setup system call information

    if (.not.associated(calc%penalty%biaslist)) then
      if (allocated(calc%penalty%biasfile)) then
        !$omp critical
        call rdensemble(calc%penalty%biasfile,nall,calc%penalty%biastmp)
        calc%penalty%biaslist => calc%penalty%biastmp
        !$omp end critical
      else
        return
      end if
    end if
    !$omp critical
!>--- printout handling
    call api_handle_output(calc,'rmsd_penalty.out',mol,pr)
!>--- populate parameters
    if (.not.allocated(calc%penalty%gradtmp)) then
      allocate (calc%penalty%gradtmp(3,mol%nat),source=0.0_wp)
      call calc%penalty%ccache%allocate(mol%nat)
    else
      calc%penalty%gradtmp(:,:) = 0.0_wp
    end if
    !$omp end critical
    if (iostatus /= 0) return

!>--- do the engrad call
    call initsignal()
    call rmsd_penalty_engrad(mol,calc%penalty,energy,grad,iostatus)
    if (iostatus /= 0) return

!>--- printout
    if (pr) then
      !> the libpvol_sp call includes the printout within libpvol-lib
      if (.not.calc%prstdout) &
      & call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data

    return
  end subroutine rmsd_engrad

!========================================================================================!

  subroutine preinit_mlip_parallel(calculations,T)
!***********************************************************************
!* Serially start one fmlip-relay server instance per OMP thread.
!* Must be called before the OMP parallel region to avoid fork() inside
!* a live thread team (triggers OMP Warning #191).
!* Input:  calculations - per-thread calcdata array (size >= T)
!*         T            - number of OMP threads / instances to start
!***********************************************************************
    implicit none
    type(calcdata),intent(inout) :: calculations(:)
    integer,intent(in)           :: T
    integer :: i,j
    do i = 1,T
      do j = 1,calculations(i)%ncalculations
        if (calculations(i)%calcs(j)%id == jobtype%mlip) then
          calculations(i)%calcs(j)%MPAR%max_threads = calculations(i)%calcs(j)%threads
          call fmlip_relay_init(calculations(i)%calcs(j)%MPAR,i)
        end if
      end do
    end do
  end subroutine preinit_mlip_parallel

!========================================================================================!

  subroutine mlip_engrad(mol,calc,energy,grad,iostatus)
!**************************************************************************
!* MLIP singlepoint through persistent python socket
!**************************************************************************
    implicit none
    type(coord) :: mol
    type(calculation_settings),target :: calc

    real(wp),intent(inout) :: energy
    real(wp),intent(inout) :: grad(3,mol%nat)
    integer,intent(out) :: iostatus

    logical :: loadnew,pr
    integer :: i,j,k,l,ich,och,io,iid
    logical :: ex
    iostatus = 0
    pr = .false.
!>--- each OpenMP thread owns one server instance; init is a no-op if running
    iid = OMP_GET_THREAD_NUM()+1
    !$omp critical
    !> forward the per-level core reservation as the server's thread cap so
    !> every instance (also ones started serially) is spawned with the same size
    calc%MPAR%max_threads = calc%threads
    call fmlip_relay_init(calc%MPAR,iid)
!>--- printout handling
    call api_handle_output(calc,'mlip.out',mol,pr)
    !$omp end critical
    if (iostatus /= 0) return

!>--- do the engrad call
    call initsignal()
    !> fmlip-relay expects the spin multiplicity (2S+1), not uhf = Nα-Nβ
    call calc%sync_multiplicity()
    call mlip_engrad_core(mol,calc%MPAR,energy,grad,iostatus, &
      &                   charge=calc%chrg,spin=calc%multiplicity,iid=iid)
    if (iostatus /= 0) return

!>--- printout
    if (pr) then
      !> the libpvol_sp call includes the printout within libpvol-lib
      if (.not.calc%prstdout) &
      & call api_print_e_grd(pr,calc%prch,mol,energy,grad)
    end if

!>--- postprocessing, getting other data

    return
  end subroutine mlip_engrad

!========================================================================================!
!########################################################################################!
!========================================================================================!
end module api_engrad
