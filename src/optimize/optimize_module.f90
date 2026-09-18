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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

!> This module wrapps the different optimization algorithms,
!> i.e., this is what can be called for geometry opt.

module optimize_module
  use iso_fortran_env,only:wp => real64
  use crest_parameters
  use crest_calculator
  use strucrd
  use ancopt_module
  use gradientdescent_module
  use rfo_module
  use lbfgs_module
  use optimize_utils
  use thermochem_module
  use hessian_reconstruct
  use newton_raphson_module
  use hr_utils
  implicit none
  private

  public :: optimize_geometry
  public :: print_opt_data

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine optimize_geometry(mol,molnew,calc,etot,grd,pr,wr,iostatus,logfile)
    !**********************************************************************
    !* Driver that dispatches to the selected geometry optimizer engine.
    !*
    !* logfile - optional name for the step-by-step trajectory logfile;
    !*           defaults to 'crestopt.log.xyz' when absent. Only written
    !*           when the engine's wr flag is set.
    !**********************************************************************
    implicit none
    !> Input
    type(coord)    :: mol
    type(calcdata) :: calc
    logical,intent(in)        :: pr
    logical,intent(in)        :: wr
    character(len=*),intent(in),optional :: logfile
    !> Output
    type(coord)   :: molnew
    integer,intent(out)       :: iostatus
    real(wp),intent(inout)    :: etot
    real(wp),intent(inout)    :: grd(3,mol%nat)
    real(wp),allocatable :: H_init(:,:),freq(:)
    integer :: nat3,io,idx,nrt,engine
    real(wp),allocatable :: hess(:),g_hess(:), g_hess_full(:,:), int_temps(:)
    logical :: pr2
    character(len=:),allocatable :: logfile_l


    !write(stdout,*) "RUNNING AN OPT"

    !> resolve the logfile name (default if not provided)
    if (present(logfile)) then
      logfile_l = logfile
    else
      logfile_l = 'crestopt.log.xyz'
    end if

    iostatus = -1
    !> do NOT overwrite original geometry
    !$omp critical
    molnew%at = mol%at
    molnew%xyz = mol%xyz
    molnew%nat = mol%nat
    molnew%wrextxyz = .true.
    if (allocated(mol%lat)) molnew%lat = mol%lat
    !$omp end critical
    nat3 = 3*mol%nat

    !> Check for optimization-individual calculation setup
    if (calc%optnewinit) then
      !$omp critical
      call calc%dealloc_params()
      !$omp end critical
    end if

    !> Check if Hessian Reconstruct is called and initialize the type
    if (calc%do_HR .or. calc%deform_opt_hess) then
      allocate (calc%chess)
      allocate (H_init(nat3,nat3))
      call calc%chess%alloc(mol%nat,calc%hu_steps,calc%chess_id_guess,calc%initialize_hr_type, calc%hr_hu_type)
    end if

    !> initial singlepoint
    if (calc%do_HR .or. calc%deform_opt_hess) calc%chess%track_step = .false. !this is not tracked to avoid duplicate
    call engrad(molnew,calc,etot,grd,iostatus)
    if (calc%do_HR .or. calc%deform_opt_hess) calc%chess%track_step = .true.
    !> optimization
    engine = calc%opt_engine
    !> periodic systems: the model-Hessian/ANC engines project out global
    !> rotation (invalid against a fixed lattice), so fall back to L-BFGS
    if (allocated(molnew%lat).and.engine /= 1) then
      if (pr) write (stdout,'(a)') '> periodic system detected: using L-BFGS optimizer'
      engine = 1
    end if
    select case (engine)
    case (0)
      call ancopt(molnew,calc,etot,grd,pr,wr,iostatus,logfile_l)
    case (1)
      !> l-bfgs goes here
      !write(stdout,'(a)') 'L-BFGS currently not implemented'
      !stop
      call lbfgs_optimize(molnew,calc,etot,grd,pr,wr,iostatus,logfile_l)
    case (2)
      !> rfo goes here
      call rfopt(molnew,calc,etot,grd,pr,wr,iostatus,logfile_l)
    case (3)
      !> newton-raphson step goes here, this is a newton step with updated hessians, i.e. quasi Newton
      call newton_raphson(molnew,calc,etot,grd,pr,wr,iostatus,logfile_l)
    case (-1)
      call gradientdescent(molnew,calc,etot,grd,pr,wr,iostatus,logfile_l)
    case default
      write (stdout,'(a)') 'Unknown optimization engine!'
      stop
    end select
    if (allocated(mol%lat).and..not.allocated(molnew%lat)) molnew%lat = mol%lat
    molnew%energy = etot

    if (calc%do_HR  .and. iostatus .eq. 0) then !> Hessian reconstruction and post-processing happen here, only do it if geometry relaxation successful
      if (calc%full_HR) then

        write (stdout,*)
        write (stdout,*) "THERMO FROM BFGS" !> This is here for full hessian reconstruct
        write (stdout,*)

        call calc_thermo_from_hess(molnew,calc%chess%H,pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,etot)

      else
        
        idx = minloc(calc%chess%order,1)
        if (minval(calc%chess%order) .eq. 0) idx = 1
        
        call initialize_hessian(calc,calc%chess%initialize_type,calc%chess%coords(idx,:,:),molnew%nat,molnew%at,calc%chess%hess(:),calc%chess%hguess,pr)
        call dhtosq(nat3,H_init,calc%chess%hess) 
        write(stdout,*)                                                                              
        write(stdout,*)"THERMO FROM INITIALIZED HESSIAN:"
        write(stdout,*) 
        call calc_thermo_from_hess(molnew,H_init,pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,etot)

        call calc%chess%construct_hessian()

        write (stdout,*)
        write (stdout,*) "THERMO FROM RECONSTRUCTED HESSIAN:" 
        write (stdout,*)

        call calc_thermo_from_hess(molnew,calc%chess%H(:,:),pr, &
        & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
        & calc%ht,calc%gt,calc%stot,etot)
      end if

      call calc%chess%dealloc()
      deallocate (calc%chess)
    end if

    !write(stdout,*) calc%g_sampling
    if (calc%g_sampling) then 
      pr2 = .false.
      !write(stdout,*) "Running gs"
      !write(stdout,*) "Energy pre correction", etot
      !$omp critical
      allocate(g_hess(nat3*(nat3+1)/2),g_hess_full(nat3,nat3))
      !$omp end critical
      call initialize_hessian(calc,calc%gs_hess_type,molnew%xyz,molnew%nat,molnew%at,g_hess,calc%chess%hguess,pr2)
      !write(*,*) "Hess Initialized"
      call dhtosq(nat3,g_hess_full,g_hess)
      !write(*,*) "Hess diagonalized"
      call calc_thermo_from_hess(molnew,g_hess_full,pr2, &
      & calc%nt,calc%temperatures,calc%ithr,calc%fscal,calc%sthr,calc%et, &
      & calc%ht,calc%gt,calc%stot,etot)
      !write(*,*) "Thermo calculated"

      !$omp critical
      allocate (int_temps(calc%nt))
      !$omp end critical

      int_temps = abs(calc%temperatures-298.15_wp)
      nrt = minloc(int_temps(:),1)
      etot = etot+calc%gt(nrt)
      !write(stdout,*) "Energy post correction", etot
    endif

    return
  end subroutine optimize_geometry

!========================================================================================!

  subroutine print_opt_data(calc,ich,natoms,tag)
    implicit none
    type(calcdata) :: calc
    integer,intent(in) :: ich
    integer,intent(in),optional :: natoms
    character(len=*),intent(in),optional :: tag
    integer :: tight,nat
    real(wp) :: ethr,gthr
    character(len=:),allocatable :: ttag
    if (present(tag)) then
      ttag = tag
    else
      ttag = ' '
    end if
    if (present(natoms)) then
      nat = natoms
    else
      nat = 0
    end if

    write (ich,'(a,a)',advance='no') ttag,'Optimization engine: '
    select case (calc%opt_engine)
    case (0)
      write (ich,'(a)') 'ANCOPT'
    case (1)
      write (ich,'(a)') 'L-BFGS'
    case (2)
      write (ich,'(a)') 'RFO'
    case (-1)
      write (ich,'(a)') 'Gradient Descent'
    case default
      write (ich,'(a)') 'Unknown'
    end select
    if (calc%opt_engine >= 0) then
      write (ich,'(a,a)',advance='no') ttag,'Hessian update type: '
      select case (calc%iupdat)
      case (0)
        write (ich,'(a)') 'BFGS'
      case (1)
        write (ich,'(a)') 'Powell'
      case (2)
        write (ich,'(a)') 'SR1'
      case (3)
        write (ich,'(a)') 'Bofill'
      case (4)
        write (ich,'(a)') 'Farkas-Schlegel'
      end select
    end if

    tight = calc%optlev
    call get_optthr(nat,tight,calc,ethr,gthr)
    write (ich,'(a,a,e10.3,a,e10.3,a)') ttag,'E/G convergence criteria: ',&
    & ethr,' Eh,',gthr,' Eh/a0'

    write (ich,'(a,a,i0)') ttag,'maximum optimization steps: ',calc%maxcycle
  end subroutine print_opt_data

!========================================================================================!
!========================================================================================!
end module optimize_module
