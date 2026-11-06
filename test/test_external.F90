module test_external
!****************************************************************************
!* Unit tests AND showcase for the externally supplied potential interface.
!*
!* These tests demonstrate how a host program that links CREST as a library
!* can inject its OWN energy+gradient routine into the calculator via the
!* jobtype%external mechanism, without CREST knowing the implementation at
!* compile time. Both registration styles and the optional opaque user
!* context (userdata) are exercised, and the host potential is driven
!* through both a single engrad call and the geometry optimizer.
!*
!* The host potential used here is a simple analytic harmonic well, so the
!* energies and gradients can be checked against a closed-form reference.
!****************************************************************************
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use crest_testmol
  use optimize_module
  implicit none
  private

  public :: collect_external

  real(wp),parameter :: thr  = 1.0e-10_wp
  real(wp),parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))

  !> force constant for the stateless origin-anchored showcase potential
  real(wp),parameter :: k_origin = 0.25_wp

  !**********************************************************************
  !* Host-side context object. This is the kind of data a host program
  !* would smuggle through CREST untouched: parameters, a reference
  !* geometry, plus a call counter to prove the callback was invoked.
  !**********************************************************************
  type :: harmonic_ctx
    real(wp) :: k = 1.0_wp           !> force constant
    real(wp),allocatable :: ref(:,:) !> reference geometry (Bohr)
    integer  :: ncalls = 0           !> number of times the callback ran
  end type harmonic_ctx

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for the external (host-supplied) potential interface
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_external(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("external SP (set_external+userdata)",test_external_sp), &
    new_unittest("external SP (manual, no userdata)  ",test_external_nodata), &
    new_unittest("external potential optimization    ",test_external_opt) &
    ]
!&>
  end subroutine collect_external

!========================================================================================!

  subroutine harmonic_engrad(nat,at,xyz,chrg,uhf,energy,gradient,iostatus,userdata)
!****************************************************************
!* A host-supplied potential matching engrad_interface.
!* Implements a harmonic well around a reference geometry that
!* is carried in via the opaque userdata context:
!*    E      = 1/2 k sum_i (x_i - ref_i)^2
!*    dE/dx  = k (x_i - ref_i)
!****************************************************************
    implicit none
    integer,intent(in)              :: nat
    integer,intent(in)              :: at(nat)
    real(wp),intent(in)             :: xyz(3,nat)
    integer,intent(in)              :: chrg
    integer,intent(in)              :: uhf
    real(wp),intent(out)            :: energy
    real(wp),intent(out)            :: gradient(3,nat)
    integer,intent(out)             :: iostatus
    class(*),intent(inout),optional :: userdata
    integer  :: i,j
    real(wp) :: d

    iostatus = 0
    energy = 0.0_wp
    gradient(:,:) = 0.0_wp

    !> this potential REQUIRES its context; signal failure otherwise
    if (.not.present(userdata)) then
      iostatus = 1
      return
    end if

    select type (u => userdata)
    type is (harmonic_ctx)
      if (.not.allocated(u%ref).or.size(u%ref,2) /= nat) then
        iostatus = 2
        return
      end if
      do i = 1,nat
        do j = 1,3
          d = xyz(j,i)-u%ref(j,i)
          energy = energy+0.5_wp*u%k*d*d
          gradient(j,i) = u%k*d
        end do
      end do
      !> prove to the test that the host context is live and mutable
      u%ncalls = u%ncalls+1
    class default
      iostatus = 3
    end select
  end subroutine harmonic_engrad

!========================================================================================!

  subroutine origin_engrad(nat,at,xyz,chrg,uhf,energy,gradient,iostatus,userdata)
!****************************************************************
!* A second host-supplied potential, this one STATELESS: an
!* isotropic harmonic well anchored at the origin. It ignores
!* the (absent) userdata, showcasing the no-context path.
!****************************************************************
    implicit none
    integer,intent(in)              :: nat
    integer,intent(in)              :: at(nat)
    real(wp),intent(in)             :: xyz(3,nat)
    integer,intent(in)              :: chrg
    integer,intent(in)              :: uhf
    real(wp),intent(out)            :: energy
    real(wp),intent(out)            :: gradient(3,nat)
    integer,intent(out)             :: iostatus
    class(*),intent(inout),optional :: userdata

    iostatus = 0
    energy = 0.5_wp*k_origin*sum(xyz**2)
    gradient(:,:) = k_origin*xyz(:,:)
  end subroutine origin_engrad

!========================================================================================!

  subroutine test_external_sp(error)
!****************************************************************
!* Register the harmonic host potential via %set_external,
!* passing an opaque context, then validate one engrad call
!* against the analytic reference and confirm the context was
!* used (ncalls incremented).
!****************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    type(harmonic_ctx),target :: ctx
    class(*),pointer :: ctxp
    real(wp) :: energy,e_ref
    real(wp),allocatable :: grad(:,:),g_ref(:,:),dvec(:,:)
    integer :: io,i,j

    !> reference structure and a deterministic displacement off it
    call get_testmol('methane',mol)
    allocate (grad(3,mol%nat),g_ref(3,mol%nat),dvec(3,mol%nat))
    do i = 1,mol%nat
      do j = 1,3
        dvec(j,i) = 0.01_wp*real((i-1)*3+j,wp)
      end do
    end do

    !> set up the host context: ref = undisplaced geometry
    ctx%k = 0.7_wp
    ctx%ref = mol%xyz
    ctx%ncalls = 0
    ctxp => ctx

    !> displace the molecule by dvec -> expected E and grad are analytic
    mol%xyz = mol%xyz+dvec
    e_ref = 0.5_wp*ctx%k*sum(dvec**2)
    g_ref = ctx%k*dvec

    !> register the external potential and run a single point
    call sett%set_external(harmonic_engrad,ctxp)
    call calc%add(sett)

    call engrad(mol,calc,energy,grad,io)
    call check(error,io,0)
    if (allocated(error)) return

    !> the callback must actually have been called (context is live)
    call check(error,ctx%ncalls,1)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=thr)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"External gradient does not match analytic reference")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
    end if

    deallocate (grad,g_ref,dvec)
  end subroutine test_external_sp

!========================================================================================!

  subroutine test_external_nodata(error)
!****************************************************************
!* Register a STATELESS host potential the "manual" way (set the
!* id and the procedure pointer directly, no userdata), and
!* validate against the analytic origin-anchored reference. This
!* also exercises the dispatch branch that omits userdata.
!****************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy,e_ref
    real(wp),allocatable :: grad(:,:),g_ref(:,:)
    integer :: io

    call get_testmol('methane',mol)
    allocate (grad(3,mol%nat),g_ref(3,mol%nat))

    !> analytic reference for the origin-anchored harmonic well
    e_ref = 0.5_wp*k_origin*sum(mol%xyz**2)
    g_ref = k_origin*mol%xyz

    !> manual registration: no %set_external, no userdata
    sett%id = jobtype%external
    sett%ext_engrad => origin_engrad
    call calc%add(sett)

    call engrad(mol,calc,energy,grad,io)
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=thr)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Stateless external gradient does not match reference")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
    end if

    deallocate (grad,g_ref)
  end subroutine test_external_nodata

!========================================================================================!

  subroutine test_external_opt(error)
!****************************************************************
!* Showcase: drive the CREST geometry optimizer with a host
!* potential. Starting from a displaced geometry, the harmonic
!* well must be minimized back to its reference (E -> 0), proving
!* the optimizer transparently consumes the external potential
!* (and that the calcdata copy carries the callback through).
!****************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol,molnew
    type(harmonic_ctx),target :: ctx
    class(*),pointer :: ctxp
    real(wp) :: energy,e_start
    real(wp),allocatable :: grad(:,:),ref(:,:)
    integer :: io,i,j
    logical :: wr,pr

    call get_testmol('methane',mol)
    allocate (grad(3,mol%nat),ref(3,mol%nat))

    !> reference = undisplaced geometry; context anchors the well there
    ref = mol%xyz
    ctx%k = 0.5_wp
    ctx%ref = ref
    ctx%ncalls = 0
    ctxp => ctx

    !> displace the start geometry (small internal distortion)
    do i = 1,mol%nat
      do j = 1,3
        mol%xyz(j,i) = mol%xyz(j,i)+0.02_wp*real(mod((i-1)*3+j,5)-2,wp)
      end do
    end do
    e_start = 0.5_wp*ctx%k*sum((mol%xyz-ref)**2)

    !> register external potential; use the gradient-descent engine so the
    !> exact analytic minimum (E=0 at ref) is reachable
    call sett%set_external(harmonic_engrad,ctxp)
    call calc%add(sett)
    calc%opt_engine = -1

    wr = .false.
    pr = .false.
    call optimize_geometry(mol,molnew,calc,energy,grad,pr,wr,io)
    call check(error,io,0)
    if (allocated(error)) return

    !> the optimizer must have repeatedly called the host potential
    if (ctx%ncalls < 2) then
      call test_failed(error,"Optimizer did not drive the external potential")
      return
    end if

    !> energy must have decreased and reached the analytic minimum (~0)
    if (energy >= e_start) then
      call test_failed(error,"External potential energy did not decrease")
      return
    end if
    call check(error,energy,0.0_wp,thr=1.0e-6_wp)
    if (allocated(error)) return

    !> optimized geometry must match the reference well
    if (any(abs(molnew%xyz-ref) > 1.0e-3_wp)) then
      call test_failed(error,"Optimized geometry did not reach the reference")
    end if

    deallocate (grad,ref)
  end subroutine test_external_opt

!========================================================================================!
!========================================================================================!
end module test_external
