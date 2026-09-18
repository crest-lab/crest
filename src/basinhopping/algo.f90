!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2024 Philipp Pracht, David Wales
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

module bh_algo_interface
  implicit none
  interface
    subroutine single_basinhopping_core(env,mol,calc,structuredump)
      use crest_data
      use crest_calculator
      use strucrd
      implicit none
      type(systemdata),intent(inout) :: env
      type(coord),intent(inout) :: mol
      type(coord),allocatable,intent(inout) :: structuredump(:)
      type(calcdata),intent(inout) :: calc
    end subroutine single_basinhopping_core
    subroutine parallel_basinhopping_core(env,mol,calc,structuredump)
      use crest_data
      use crest_calculator
      use strucrd
      implicit none
      !> INPUT/OUTPUT
      type(systemdata),intent(inout) :: env
      type(coord),intent(inout) :: mol
      type(coord),allocatable,intent(inout) :: structuredump(:)
      type(calcdata),intent(in) :: calc
    end subroutine parallel_basinhopping_core
  end interface
end module bh_algo_interface

!================================================================================!
!================================================================================!
!================================================================================!

subroutine crest_basinhopping(env,tim)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use cregen_interface,only:unionizeEnsembles,cregen_irmsd_sort
  use optimize_module
  use bh_module
  use bh_algo_interface
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn,mciter
  logical :: pr,wr
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  integer :: nall
  type(coord),allocatable :: structuredump(:)
  integer,allocatable :: groups(:)
  logical :: parallel
  character(len=80) :: atmp
  character(len=*),parameter :: trjf = 'crest_quenched.xyz'
!========================================================================================!
  write (stdout,*)
  write (stdout,*) " ____            _       _   _                   _             "
  write (stdout,*) "| __ )  __ _ ___(_)_ __ | | | | ___  _ __  _ __ (_)_ __   __ _ "
  write (stdout,*) "|  _ \ / _` / __| | '_ \| |_| |/ _ \| '_ \| '_ \| | '_ \ / _` |"
  write (stdout,*) "| |_) | (_| \__ \ | | | |  _  | (_) | |_) | |_) | | | | | (_| |"
  write (stdout,*) "|____/ \__,_|___/_|_| |_|_| |_|\___/| .__/| .__/|_|_| |_|\__, |"
  write (stdout,*) "                                    |_|   |_|            |___/ "
  write (stdout,*) ""
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()

  calc = env%calc
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!
!>--- print calculation info
  call calc%info(stdout)
  write (stdout,'(a)') '> Geometry optimization settings:'
  call print_opt_data(calc,stdout,natoms=mol%nat,tag=' : ')
  write (stdout,*)

!>--- singlepoint of input structure
  allocate (grad(3,mol%nat),source=0.0_wp)
  call engrad(mol,calc,energy,grad,io)
  mol%energy = energy  !> we need this to start the Markov-chain

!==========================================================================================!
  parallel = .false.
  if (allocated(env%bh_ref)) then
    parallel = env%bh_ref%parallel
  end if

!=========================================================================================!
  call tim%start(14,'Basin-Hopping (BH)')

  if (parallel) then
    call parallel_basinhopping_core(env,mol,calc,structuredump)
  else
    call single_basinhopping_core(env,mol,calc,structuredump)
  end if
!>--- dump saved minima
  nall = size(structuredump,1)
  open (newunit=ich,file=trjf)
  call wrensemble(ich,nall,structuredump)
  close (ich)

  if (io == 0) then
    write (stdout,*)
    write (stdout,*) 'BH run completed successfully'
    write (stdout,*) 'Successfull quenches written to ',trjf
  else
    write (stdout,*) 'WARNING: BH run terminated ABNORMALLY'
    env%iostatus_meta = status_failed
  end if
  call tim%stop(14)

  write (stdout,*)
  call smallhead('Final Ensemble Sorting (iRMSD)')
  allocate (groups(nall),source=0)
  env%confgo = .true.
  call cregen_irmsd_sort(env,nall,structuredump,groups,allcanon=.false.,printlvl=2)

  if (allocated(groups)) deallocate (groups)
  if (allocated(structuredump)) deallocate (structuredump)
  return
end subroutine crest_basinhopping

subroutine single_basinhopping_core(env,mol,calc,structuredump)
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use cregen_interface,only:unionizeEnsembles
  use optimize_module
  use molbuilder_classify
  use bh_module
  implicit none
  type(systemdata),intent(inout) :: env
  type(coord),intent(inout) :: mol
  type(coord),allocatable,intent(inout) :: structuredump(:)
  type(calcdata),intent(inout) :: calc
!========================================================================================!
  integer :: i,j,k,l,io,ich,T,Tn,mciter
  logical :: pr,wr
  type(bh_class) :: bh
  integer :: nall
!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
!========================================================================================!
!>--- actual basin hopping
  if (allocated(env%bh_ref)) then
    bh = env%bh_ref
    call bh%init()
  else
    call bh%init(300.0_wp,200,20)
    bh%stepsize(1) = 1.0_wp
  end if
  select case (bh%steptype)
  case (1,2) !> internals, dihedral only
    write (stdout,'(a)') '> Setting up internal coordinates for input molecule:'
    call setup_classify(mol,bh%molc)
    call functional_group_classify(bh%molc)
    call bh%molc%get_zmat(.true.)
    call bh%molc%print_zmat(stdout)
    write (stdout,*)
    call bh%molc%check_dihedrals()
  end select
  bh%id = 0
  if (allocated(env%refine_queue)) then
    bh%refine_queue = env%refine_queue
  end if

  nall = 0
  do mciter = 1,bh%maxiter
    if (bh%maxiter > 1) call printiter3('Basin-Hopping Epoch',mciter)
    call bh%newiter()
    call mc(calc,mol,bh,io,verbosity=2)

    if (io .eq. 0) then
      write (stdout,'(a)') 'New structures will be appended to memory ...'
      call unionizeEnsembles(nall,structuredump,bh%saved,bh%structures, &
      &                      ethr=bh%ethr,rthr=bh%rthr)
    else
      write (stdout,'(a)') 'Skipping run with failed initial quench ...'
    end if
    write (stdout,'(a,i0,a)') 'Currently ',nall,' structures saved!'
  end do
  return
end subroutine single_basinhopping_core

subroutine parallel_basinhopping_core(env,mol,calc,structuredump)
!**************************************************************************
!* subroutine parallel_basinhopping_core
!* Perform multiple independent BH runs from a single given starting point
!* Ensembles are unified at the end and returned via structuedump
!**************************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use cregen_interface,only:unionizeEnsembles
  use optimize_module
  use bh_module
  use iomod,only:is_terminal
  use molbuilder_classify
  implicit none
  !> INPUT/OUTPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(inout) :: mol
  type(coord),allocatable,intent(inout) :: structuredump(:)
  type(calcdata),intent(in) :: calc
!========================================================================================!
  !> LOCAL
  integer :: i,j,k,l,io,ich,T,Tn,mciter
  logical :: pr,wr
  type(calcdata),allocatable :: calcp(:)
  type(bh_class),allocatable :: bhp(:)
  type(coord),allocatable    :: mols(:)
  real(wp) :: energy
  integer :: nall,verbose,iostatus
  character(len=128) :: tag
  type(mollist),allocatable :: dumplist(:)

  call new_ompautoset(env,'auto',0,T,Tn)
  !======================================================================================!
  !> THIS IS THE PARALLEL IMPORTANT BIT
  !======================================================================================!
!>--- allocate temporary spaces for parallel usage
  allocate (mols(T),source=mol)
  allocate (bhp(T))
  allocate (calcp(T),source=calc)
  allocate (dumplist(T))
  if (allocated(env%bh_ref)) then
    do K = 1,T
      bhp(K) = env%bh_ref
      call bhp(K)%init()
    end do
  else
    do K = 1,T
      call bhp(K)%init(300.0_wp,200,20)
      bhp(K)%stepsize(1) = 1.0_wp
    end do
  end if
  do K = 1,T
    call calcp(K)%copy(calc)
    bhp(K)%id = K-1
    !$omp critical
    select case (bhp(K)%steptype)
    case (1,2) !> internals, dihedral only
      if (K == 1) write (stdout,'(a)') '> Setting up internal coordinates for input molecule:'
      call setup_classify(mol,bhp(K)%molc)
      call functional_group_classify(bhp(K)%molc)
      call bhp(K)%molc%get_zmat(.true.)
      if (K == 1) call bhp(K)%molc%print_zmat(stdout)
      if (K == 1) write (stdout,*)
      call bhp(K)%molc%check_dihedrals()
    end select
    if (allocated(env%refine_queue)) then
      bhp(K)%refine_queue = env%refine_queue
    end if
    !$omp end critical
  end do

!>--- pre-start server-based calculators before forking OMP threads
  call preinit_mlip_parallel(calcp,T)

  write (stdout,'(a)') '> Starting parallel Basin-Hopping execution'
  write (stdout,*)

  do mciter = 1,bhp(1)%maxiter
    !$omp critical
    write (tag,'(a,i0,a)') 'Basin-Hopping Epoch'
    if (bhp(1)%maxiter > 1) call printiter3(trim(tag),mciter)
    !$omp end critical
    !$omp parallel do default(shared) private(K, mciter, iostatus) schedule(dynamic)
    do K = 1,T
      call bhp(K)%newiter()
      call mc(calcp(K),mols(K),bhp(K),iostatus,verbosity=1)
      if (iostatus .eq. 0) then
        !$omp critical
        write (stdout,'(a)') 'New structures will be appended to memory ...'
        !$omp end critical
        call unionizeEnsembles(dumplist(K)%nall,dumplist(K)%structure, &
        &                      bhp(K)%saved,bhp(K)%structures, &
        &                      ethr=bhp(K)%ethr,rthr=bhp(K)%rthr)
      else
        !$omp critical
        write (stdout,'(a)') 'Skipping run with failed initial quench ...'
        !$omp end critical
      end if
      !$omp critical
      write (stdout,'(a,i0,a,i0,a)') 'Currently ',dumplist(K)%nall, &
      &      ' structures saved (BH[',bhp(K)%id,'])!'
      !$omp end critical
    end do
    !$omp end parallel do

    !> Do things here (?)
  end do

  write (stdout,*)
  write (stdout,'(a)') 'Parallel BH runs done!'
  write (stdout,'(a)') 'Collecting structures in one ensemble ...'
  nall = 0
  do K = 1,T
    call unionizeEnsembles(nall,structuredump, &
    &                      dumplist(K)%nall,dumplist(K)%structure, &
    &                      ethr=bhp(K)%ethr,rthr=bhp(K)%rthr)
  end do
  write (stdout,'(a,i0,a)') 'Total of ',nall,' structures remain.'
!=======================================================================================!
!> PARALLEL BIT END
!=======================================================================================!
  return
end subroutine parallel_basinhopping_core
