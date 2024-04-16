module test_optimization
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use testmol
  use optimize_module
  implicit none
  private

  public :: collect_optimization

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
  real(wp),parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for using geometry optimization routines in CREST
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_optimization(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
#ifdef WITH_GFNFF
    new_unittest("Compiled gfnff subproject     ",test_compiled_gfnff), &
    new_unittest("optimization (ANCOPT)         ",test_ancopt), &
    new_unittest("optimization (ANCOPT,sspevx)  ",test_ancoptsmall), &
    new_unittest("optimization (grad. descent)  ",test_gradientdescent), &
    new_unittest("optimization (RFO)            ",test_rfo) &
#else
    new_unittest("Compiled gfnff subproject",test_compiled_gfnff,should_fail=.true.) &
#endif
    ]
!&>
  end subroutine collect_optimization

  subroutine test_compiled_gfnff(error)
    type(error_type),allocatable,intent(out) :: error
#ifndef WITH_GFNFF
    write(*,'("       ...")') 'gfnff not compiled, expecting fail.'
    allocate (error)
#endif
  end subroutine test_compiled_gfnff

!========================================================================================!

  subroutine test_ancopt(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol,molnew
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
    logical :: wr,pr
!&<
    real(wp),parameter :: e_ref = -4.677663337455959_wp
!&>

    !> setup
    call sett%create('gfnff')
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    wr = .false.
    pr = .false.
    call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-6_wp)
    if (allocated(error)) return

    deallocate (grad)
  end subroutine test_ancopt

!========================================================================================!

  subroutine test_ancoptsmall(error)
!*****************************************************
!* Test ANCOPT with a small molecule to trigger exact
!* Hessian eigenvalue calculation via LAPACK's sspevx
!*****************************************************
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol,molnew
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
    logical :: wr,pr
!&<
    real(wp),parameter :: e_ref = -0.630873310757319_wp
!&>

    !> setup
    call sett%create('gfnff')
    call calc%add(sett)
    call get_testmol('methane',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    wr = .false.
    pr = .false.
    call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-6_wp)
    if (allocated(error)) return

    deallocate (grad)
  end subroutine test_ancoptsmall


!========================================================================================!
  subroutine test_gradientdescent(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol,molnew
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
    logical :: wr,pr
!&<
    real(wp),parameter :: e_ref = -4.677587929227879_wp
!&>

    !> setup
    call sett%create('gfnff')
    call calc%add(sett)
    calc%opt_engine = -1
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    wr = .false.
    pr = .false.
    call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-6_wp)
    if (allocated(error)) return

    deallocate (grad)
  end subroutine test_gradientdescent

!========================================================================================!
  subroutine test_rfo(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol,molnew
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
    logical :: wr,pr
!&<
    real(wp),parameter :: e_ref = -4.677662006957390_wp
!&>

    !> setup
    call sett%create('gfnff')
    call calc%add(sett)
    calc%opt_engine = 2
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    wr = .false.
    pr = .false.
    call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-6_wp)
    if (allocated(error)) return

    deallocate (grad)
  end subroutine test_rfo

!========================================================================================!
!========================================================================================!
end module test_optimization
