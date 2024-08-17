module test_helpers
  use testdrive,only:error_type,check,test_failed
  use crest_parameters
  use crest_testmol
  use crest_calculator
  use strucrd
  implicit none
  private

  public :: test_e,test_g

!========================================================================================!
!========================================================================================!
contains  !> Helper routines for testing
!========================================================================================!
!========================================================================================!

  subroutine test_e(mol,calc,error,ref)
!***********************************************
!* calculate energy and cross-check with
!* a given reference
!***********************************************
    !> Error handling
    type(error_type),allocatable,intent(out) :: error
    !> calculation data
    type(calcdata),intent(inout) :: calc
    !> molecule
    type(coord),intent(inout) :: mol
    !> reference energy
    real(wp),intent(in) :: ref
    !> LOCAL
    real(wp),parameter :: thr = sqrt(epsilon(1.0_wp))
    real(wp),allocatable :: gradient(:,:)
    integer :: io
    real(wp) :: energy

    energy=0.0_wp

    !> calculation
    allocate (gradient(3,mol%nat),source=0.0_wp)
    call engrad(mol,calc,energy,gradient,io)

    if (abs(energy-ref) > thr) then
      !call test_failed(error,"energy does not match reference")
      call calc%info(stdout)
      call check(error,energy,ref,thr=thr)
    end if

  end subroutine test_e

  subroutine test_g(mol,calc,error)
!***********************************************
!* calculate numerical gradient and cross-check
!* with the one returned by calculator
!***********************************************
    !> Error handling
    type(error_type),allocatable,intent(out) :: error
    !> calculation data
    type(calcdata),intent(inout) :: calc
    !> molecule
    type(coord),intent(inout) :: mol
    !> LOCAL
    real(wp),parameter :: thr = 1.0E-6_wp
    real(wp),allocatable :: gradient(:,:),numg(:,:)
    integer :: io
    real(wp) :: energy

    !> calculation
    allocate (gradient(3,mol%nat),source=0.0_wp)
    call engrad(mol,calc,energy,gradient,io)

    call numgrad(mol,calc,numgradient=numg)

    if (any(abs(gradient-numg) > thr)) then
      call test_failed(error,"Gradient does not match")
      call calc%info(stdout)
      print '(3es20.13)',gradient
      print '(a)',"---"
      print '(3es20.13)',numg
      print '(a)',"---"
      print '(3es20.13)',gradient-numg
    end if
  end subroutine test_g

!========================================================================================!
!========================================================================================!
end module test_helpers
