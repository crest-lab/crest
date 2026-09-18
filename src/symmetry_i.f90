!> symmetry_i.f90
!> Brute force symmetry analyzer - Fortran module
!>
!> Fortran conversion of the original C code, 2026 Philipp Pracht
!> 
!> Original C code: (C) 1996, 2003 S. Patchkovskii
!>
!>
!> This program is free software; you can redistribute it and/or modify
!> it under the terms of the GNU General Public License as published by
!> the Free Software Foundation; either version 2 of the License, or
!> (at your option) any later version.

module symmetry_i
  use iso_fortran_env,only:wp => real64
  implicit none
  private

  ! Public interface
  public :: schoenflies
  public :: getsym
  !public :: symmetry_element,atom_t,symmetry_state_t

  !> Mathematical constants
  real(wp),parameter :: PI = 3.14159265358979323846d0
  integer,parameter :: DIMENSION = 3
  integer,parameter :: MAXPARAM = 7

  !> Atom type
  type :: atom_t
    integer :: atom_type
    real(wp) :: x(DIMENSION)
  end type atom_t

  !> Symmetry element type
  type :: symmetry_element
    integer :: transform_type  ! 1=mirror, 2=invert, 3=rotate, 4=rotate_reflect
    integer,allocatable :: transform(:)
    integer :: order
    integer :: nparam
    real(wp) :: maxdev
    real(wp) :: distance
    real(wp) :: normal(DIMENSION)
    real(wp) :: direction(DIMENSION)
  end type symmetry_element

  !> Point group type
  type :: point_group
    character(len=8) :: group_name
    character(len=64) :: symmetry_code
  end type point_group

  !> Number of point groups in the lookup table
  integer,parameter :: PointGroupsCount = 60

  !> Compile-time point group lookup table
  type(point_group),parameter :: PointGroups(PointGroupsCount) = [                &
  &  point_group("C1",""),                                                        &
  &  point_group("Cs","(sigma) "),                                                &
  &  point_group("Ci","(i) "),                                                    &
  &  point_group("C2","(C2) "),                                                   &
  &  point_group("C3","(C3) "),                                                   &
  &  point_group("C4","(C4) (C2) "),                                              &
  &  point_group("C5","(C5) "),                                                   &
  &  point_group("C6","(C6) (C3) (C2) "),                                         &
  &  point_group("C7","(C7) "),                                                   &
  &  point_group("C8","(C8) (C4) (C2) "),                                         &
  &  point_group("D2","3*(C2) "),                                                 &
  &  point_group("D3","(C3) 3*(C2) "),                                            &
  &  point_group("D4","(C4) 5*(C2) "),                                            &
  &  point_group("D5","(C5) 5*(C2) "),                                            &
  &  point_group("D6","(C6) (C3) 7*(C2) "),                                       &
  &  point_group("D7","(C7) 7*(C2) "),                                            &
  &  point_group("D8","(C8) (C4) 9*(C2) "),                                       &
  &  point_group("C2v","(C2) 2*(sigma) "),                                        &
  &  point_group("C3v","(C3) 3*(sigma) "),                                        &
  &  point_group("C4v","(C4) (C2) 4*(sigma) "),                                   &
  &  point_group("C5v","(C5) 5*(sigma) "),                                        &
  &  point_group("C6v","(C6) (C3) (C2) 6*(sigma) "),                              &
  &  point_group("C7v","(C7) 7*(sigma) "),                                        &
  &  point_group("C8v","(C8) (C4) (C2) 8*(sigma) "),                              &
  &  point_group("C2h","(i) (C2) (sigma) "),                                      &
  &  point_group("C3h","(C3) (S3) (sigma) "),                                     &
  &  point_group("C4h","(i) (C4) (C2) (S4) (sigma) "),                            &
  &  point_group("C5h","(C5) (S5) (sigma) "),                                     &
  &  point_group("C6h","(i) (C6) (C3) (C2) (S6) (S3) (sigma) "),                  &
  &  point_group("C7h","(C7) (S7) (sigma) "),                                     &
  &  point_group("C8h","(i) (C8) (C4) (C2) (S8) (S4) (sigma) "),                  &
  &  point_group("D2h","(i) 3*(C2) 3*(sigma) "),                                  &
  &  point_group("D3h","(C3) 3*(C2) (S3) 4*(sigma) "),                            &
  &  point_group("D4h","(i) (C4) 5*(C2) (S4) 5*(sigma) "),                        &
  &  point_group("D5h","(C5) 5*(C2) (S5) 6*(sigma) "),                            &
  &  point_group("D6h","(i) (C6) (C3) 7*(C2) (S6) (S3) 7*(sigma) "),              &
  &  point_group("D7h","(C7) 7*(C2) (S7) 8*(sigma) "),                            &
  &  point_group("D8h","(i) (C8) (C4) 9*(C2) (S8) (S4) 9*(sigma) "),              &
  &  point_group("D2d","3*(C2) (S4) 2*(sigma) "),                                 &
  &  point_group("D3d","(i) (C3) 3*(C2) (S6) 3*(sigma) "),                        &
  &  point_group("D4d","(C4) 5*(C2) (S8) 4*(sigma) "),                            &
  &  point_group("D5d","(i) (C5) 5*(C2) (S10) 5*(sigma) "),                       &
  &  point_group("D6d","(C6) (C3) 7*(C2) (S12) (S4) 6*(sigma) "),                 &
  &  point_group("D7d","(i) (C7) 7*(C2) (S14) 7*(sigma) "),                       &
  &  point_group("D8d","(C8) (C4) 9*(C2) (S16) 8*(sigma) "),                      &
  &  point_group("S4","(C2) (S4) "),                                              &
  &  point_group("S6","(i) (C3) (S6) "),                                          &
  &  point_group("S8","(C4) (C2) (S8) "),                                         &
  &  point_group("T","4*(C3) 3*(C2) "),                                           &
  &  point_group("Th","(i) 4*(C3) 3*(C2) 4*(S6) 3*(sigma) "),                     &
  &  point_group("Td","4*(C3) 3*(C2) 3*(S4) 6*(sigma) "),                         &
  &  point_group("O","3*(C4) 4*(C3) 9*(C2) "),                                    &
  &  point_group("Oh","(i) 3*(C4) 4*(C3) 9*(C2) 4*(S6) 3*(S4) 9*(sigma) "),       &
  &  point_group("Cinfv","(Cinf) (sigma) "),                                      &
  &  point_group("Dinfh","(i) (Cinf) (C2) 2*(sigma) "),                           &
  &  point_group("I","6*(C5) 10*(C3) 15*(C2) "),                                  &
  &  point_group("Ih","(i) 6*(C5) 10*(C3) 15*(C2) 6*(S10) 10*(S6) 15*(sigma) "),  &
  &  point_group("Kh","(i) (Cinf) (sigma) "),                                     &
  &  point_group("",""),                                                          &
  &  point_group("","")]

  !> All symmetry-analysis state collected in one derived type
  type,public :: symmetry_state_t
    ! Tolerance / control
    real(wp) :: ToleranceSame = 1.0d-3
    real(wp) :: TolerancePrimary = 5.0d-2
    real(wp) :: ToleranceFinal = 1.0d-4
    real(wp) :: MaxOptStep = 5.0d-1
    real(wp) :: MinOptStep = 1.0d-7
    real(wp) :: GradientStep = 1.0d-7
    real(wp) :: OptChangeThreshold = 1.0d-10
    integer  :: verbose = 0
    integer  :: MaxAxisOrder = 20
    integer  :: MaxOptCycles = 200
    integer  :: OptChangeHits = 5
    ! Geometry / working data
    real(wp)              :: CenterOfSomething(3) = 0.0_wp
    real(wp),allocatable  :: DistanceFromCenter(:)
    integer               :: AtomsCount = 0
    type(atom_t),allocatable :: Atoms(:)
    ! Symmetry elements
    integer                            :: PlanesCount = 0
    type(symmetry_element),allocatable :: Planes(:)
    type(symmetry_element)             :: MolecularPlane
    logical                            :: MolecularPlaneExists = .false.
    integer                            :: InversionCentersCount = 0
    type(symmetry_element),allocatable :: InversionCenters(:)
    integer                            :: NormalAxesCount = 0
    type(symmetry_element),allocatable :: NormalAxes(:)
    integer                            :: ImproperAxesCount = 0
    type(symmetry_element),allocatable :: ImproperAxes(:)
    integer,allocatable                :: NormalAxesCounts(:)
    integer,allocatable                :: ImproperAxesCounts(:)
    integer                            :: BadOptimization = 0
    character(len=256)                 :: SymmetryCode = ""
    character(len=8)                   :: MaxRotAxis = ""
    ! Statistics
    integer(8) :: StatTotal = 0
    integer(8) :: StatEarly = 0
    integer(8) :: StatPairs = 0
    integer(8) :: StatDups = 0
    integer(8) :: StatOrder = 0
    integer(8) :: StatOpt = 0
    integer(8) :: StatAccept = 0
  end type symmetry_state_t

! ══════════════════════════════════════════════════════════════════════════════
contains    !> MODULE PROCEDURES START HERE
! ══════════════════════════════════════════════════════════════════════════════

  !> Initialise (reset) a symmetry_state_t to its defaults
  subroutine init_symmetry_state(state)
    type(symmetry_state_t),intent(out) :: state
    ! intent(out) resets all scalar fields to their type-definition defaults
    ! and deallocates every allocatable component (transitively).
    ! Explicitly destroy the non-allocatable MolecularPlane's inner allocatable.
    call destroy_symmetry_element(state%MolecularPlane)
  end subroutine init_symmetry_state

  !> Square function
  pure real(wp) function pow2(x)
    real(wp),intent(in) :: x
    pow2 = x*x
  end function pow2

  !> Allocate a symmetry element
  subroutine alloc_symmetry_element(state,elem)
    type(symmetry_state_t),intent(in) :: state
    type(symmetry_element),intent(out) :: elem
    integer :: i

    allocate (elem%transform(state%AtomsCount))
    do i = 1,state%AtomsCount
      elem%transform(i) = state%AtomsCount+1  ! Impossible value
    end do
    elem%order = 0
    elem%nparam = 0
    elem%maxdev = 0.0d0
    elem%distance = 0.0d0
    elem%normal = 0.0d0
    elem%direction = 0.0d0
    elem%transform_type = 0
  end subroutine alloc_symmetry_element

  !> Deallocate a symmetry element
  subroutine destroy_symmetry_element(elem)
    type(symmetry_element),intent(inout) :: elem
    if (allocated(elem%transform)) deallocate (elem%transform)
  end subroutine destroy_symmetry_element

  !> Mirror an atom through a plane
  subroutine mirror_atom(plane,from_atom,to_atom)
    type(symmetry_element),intent(in) :: plane
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom
    integer :: i
    real(wp) :: r

    r = plane%distance
    do i = 1,DIMENSION
      r = r-from_atom%x(i)*plane%normal(i)
    end do

    to_atom%atom_type = from_atom%atom_type
    do i = 1,DIMENSION
      to_atom%x(i) = from_atom%x(i)+2.0d0*r*plane%normal(i)
    end do
  end subroutine mirror_atom

  !> Invert an atom through a center
  subroutine invert_atom(center,from_atom,to_atom)
    type(symmetry_element),intent(in) :: center
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom
    integer :: i

    to_atom%atom_type = from_atom%atom_type
    do i = 1,DIMENSION
      to_atom%x(i) = 2.0d0*center%distance*center%normal(i)-from_atom%x(i)
    end do
  end subroutine invert_atom

  !> Rotate an atom around an axis
  subroutine rotate_atom(axis,from_atom,to_atom)
    type(symmetry_element),intent(in) :: axis
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom
    real(wp) :: x(3),y(3),a(3),b(3),c(3)
    real(wp) :: angle,a_sin,a_cos,dot_val
    integer :: i

    if (axis%order /= 0) then
      angle = 2.0d0*PI/dble(axis%order)
    else
      angle = 1.0d0
    end if
    a_sin = sin(angle)
    a_cos = cos(angle)

    do i = 1,3
      x(i) = from_atom%x(i)-axis%distance*axis%normal(i)
    end do

    dot_val = 0.0d0
    do i = 1,3
      dot_val = dot_val+x(i)*axis%direction(i)
    end do

    do i = 1,3
      a(i) = axis%direction(i)*dot_val
    end do

    do i = 1,3
      b(i) = x(i)-a(i)
    end do

    c(1) = b(2)*axis%direction(3)-b(3)*axis%direction(2)
    c(2) = b(3)*axis%direction(1)-b(1)*axis%direction(3)
    c(3) = b(1)*axis%direction(2)-b(2)*axis%direction(1)

    do i = 1,3
      y(i) = a(i)+b(i)*a_cos+c(i)*a_sin
    end do

    do i = 1,3
      to_atom%x(i) = y(i)+axis%distance*axis%normal(i)
    end do
    to_atom%atom_type = from_atom%atom_type
  end subroutine rotate_atom

  !> Rotate and reflect an atom (improper rotation)
  subroutine rotate_reflect_atom(axis,from_atom,to_atom)
    type(symmetry_element),intent(in) :: axis
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom
    real(wp) :: x(3),y(3),a(3),b(3),c(3)
    real(wp) :: angle,a_sin,a_cos,dot_val
    integer :: i

    angle = 2.0d0*PI/dble(axis%order)
    a_sin = sin(angle)
    a_cos = cos(angle)

    do i = 1,3
      x(i) = from_atom%x(i)-axis%distance*axis%normal(i)
    end do

    dot_val = 0.0d0
    do i = 1,3
      dot_val = dot_val+x(i)*axis%direction(i)
    end do

    do i = 1,3
      a(i) = axis%direction(i)*dot_val
    end do

    do i = 1,3
      b(i) = x(i)-a(i)
    end do

    c(1) = b(2)*axis%direction(3)-b(3)*axis%direction(2)
    c(2) = b(3)*axis%direction(1)-b(1)*axis%direction(3)
    c(3) = b(1)*axis%direction(2)-b(2)*axis%direction(1)

    do i = 1,3
      y(i) = -a(i)+b(i)*a_cos+c(i)*a_sin
    end do

    do i = 1,3
      to_atom%x(i) = y(i)+axis%distance*axis%normal(i)
    end do
    to_atom%atom_type = from_atom%atom_type
  end subroutine rotate_reflect_atom

  !> Transform atom based on element type
  subroutine transform_atom(elem,from_atom,to_atom)
    type(symmetry_element),intent(in) :: elem
    type(atom_t),intent(in) :: from_atom
    type(atom_t),intent(out) :: to_atom

    select case (elem%transform_type)
    case (1)
      call mirror_atom(elem,from_atom,to_atom)
    case (2)
      call invert_atom(elem,from_atom,to_atom)
    case (3)
      call rotate_atom(elem,from_atom,to_atom)
    case (4)
      call rotate_reflect_atom(elem,from_atom,to_atom)
    case default
      to_atom = from_atom
    end select
  end subroutine transform_atom

  !> Establish pairs of atoms related by symmetry
  function establish_pairs(state,elem) result(status)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(inout) :: elem
    integer :: status
    integer :: i,j,k,best_j
    logical,allocatable :: atom_used(:)
    real(wp) :: distance,best_distance
    type(atom_t) :: symmetric

    status = 0
    allocate (atom_used(state%AtomsCount))
    atom_used = .false.

    do i = 1,state%AtomsCount
      if (elem%transform(i) > state%AtomsCount) then
        call transform_atom(elem,state%Atoms(i),symmetric)
        best_j = i
        best_distance = 2.0d0*state%TolerancePrimary

        do j = 1,state%AtomsCount
          if (state%Atoms(j)%atom_type /= symmetric%atom_type.or.atom_used(j)) cycle

          distance = 0.0d0
          do k = 1,DIMENSION
            distance = distance+pow2(symmetric%x(k)-state%Atoms(j)%x(k))
          end do
          distance = sqrt(distance)

          if (distance < best_distance) then
            best_j = j
            best_distance = distance
          end if
        end do

        if (best_distance > state%TolerancePrimary) then
          deallocate (atom_used)
          status = -1
          return
        end if

        elem%transform(i) = best_j
        atom_used(best_j) = .true.
      end if
    end do

    deallocate (atom_used)
  end function establish_pairs

  !> Check if transformation order is correct
  function check_transform_order(state,elem) result(status)
    type(symmetry_state_t),intent(in) :: state
    type(symmetry_element),intent(in) :: elem
    integer :: status
    integer :: i,j,k

    status = 0

    do i = 1,state%AtomsCount
      if (elem%transform(i) == i) cycle

      if (elem%transform_type == 4) then  ! rotate_reflect
        j = elem%transform(i)
        if (elem%transform(j) == i) cycle
      end if

      k = elem%transform(i)
      do j = elem%order-1,1,-1
        if (k == i) then
          status = -1
          return
        end if
        k = elem%transform(k)
      end do

      if (k /= i.and.elem%transform_type == 4) then
        do j = elem%order,1,-1
          if (k == i) then
            status = -1
            return
          end if
          k = elem%transform(k)
        end do
      end if

      if (k /= i) then
        status = -1
        return
      end if
    end do
  end function check_transform_order

  !> Check if two transforms are the same
  function same_transform(state,a,b) result(is_same)
    type(symmetry_state_t),intent(in) :: state
    type(symmetry_element),intent(in) :: a,b
    logical :: is_same
    integer :: i,j,code

    is_same = .false.

    if (a%order /= b%order.or.a%nparam /= b%nparam.or. &
        a%transform_type /= b%transform_type) return

    code = 1
    do i = 1,state%AtomsCount
      if (a%transform(i) /= b%transform(i)) then
        code = 0
        exit
      end if
    end do

    if (code == 0.and.a%order > 2) then
      do i = 1,state%AtomsCount
        j = a%transform(i)
        if (b%transform(j) /= i) return
      end do
      is_same = .true.
      return
    end if

    is_same = (code == 1)
  end function same_transform

  !> Check transform quality
  function check_transform_quality(state,elem) result(status)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(inout) :: elem
    integer :: status
    integer :: i,j,k
    type(atom_t) :: symmetric
    real(wp) :: r,max_r

    status = 0
    max_r = 0.0d0

    do i = 1,state%AtomsCount
      j = elem%transform(i)
      call transform_atom(elem,state%Atoms(i),symmetric)

      r = 0.0d0
      do k = 1,DIMENSION
        r = r+pow2(symmetric%x(k)-state%Atoms(j)%x(k))
      end do
      r = sqrt(r)

      if (r > state%ToleranceFinal) then
        status = -1
        return
      end if
      if (r > max_r) max_r = r
    end do

    elem%maxdev = max_r
  end function check_transform_quality

  !> Evaluate optimization target function
  function eval_optimization_target_function(state,elem,finish) result(target)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(inout) :: elem
    logical,intent(out),optional :: finish
    real(wp) :: target
    integer :: i,j,k
    type(atom_t) :: symmetric
    real(wp) :: r,maxr

    ! Normalize normal vector
    if (elem%nparam >= 4) then
      r = 0.0d0
      do k = 1,DIMENSION
        r = r+elem%normal(k)*elem%normal(k)
      end do
      r = sqrt(r)
      if (r < state%ToleranceSame) then
        write (*,*) "Normal collapsed!"
        stop
      end if
      elem%normal = elem%normal/r
      if (elem%distance < 0.0d0) then
        elem%distance = -elem%distance
        elem%normal = -elem%normal
      end if
    end if

    ! Normalize direction vector
    if (elem%nparam >= 7) then
      r = 0.0d0
      do k = 1,DIMENSION
        r = r+elem%direction(k)*elem%direction(k)
      end do
      r = sqrt(r)
      if (r < state%ToleranceSame) then
        write (*,*) "Direction collapsed!"
        stop
      end if
      elem%direction = elem%direction/r
    end if

    target = 0.0d0
    maxr = 0.0d0

    do i = 1,state%AtomsCount
      call transform_atom(elem,state%Atoms(i),symmetric)
      j = elem%transform(i)

      r = 0.0d0
      do k = 1,DIMENSION
        r = r+pow2(state%Atoms(j)%x(k)-symmetric%x(k))
      end do
      if (r > maxr) maxr = r
      target = target+r
    end do

    if (present(finish)) then
      finish = (sqrt(maxr) < state%ToleranceFinal)
    end if
  end function eval_optimization_target_function

  !> Get parameters from element
  subroutine get_params(elem,values)
    type(symmetry_element),intent(in) :: elem
    real(wp),intent(out) :: values(MAXPARAM)

    values(1) = elem%distance
    values(2:4) = elem%normal(1:3)
    if (elem%nparam >= 7) then
      values(5:7) = elem%direction(1:3)
    end if
  end subroutine get_params

  !> Set parameters to element
  subroutine set_params(elem,values)
    type(symmetry_element),intent(inout) :: elem
    real(wp),intent(in) :: values(MAXPARAM)

    elem%distance = values(1)
    elem%normal(1:3) = values(2:4)
    if (elem%nparam >= 7) then
      elem%direction(1:3) = values(5:7)
    end if
  end subroutine set_params

  !> Optimize transformation parameters
  subroutine optimize_transformation_params(state,elem)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(inout) :: elem
    real(wp) :: values(MAXPARAM),grad(MAXPARAM),force(MAXPARAM),step(MAXPARAM)
    real(wp) :: f,fold,fnew,fnew2,fdn,fup,snorm
    real(wp) :: a,b,x
    integer :: vars,cycle,i,hits
    logical :: finish

    values = 0.0_wp
    grad = 0.0_wp
    force = 0.0_wp
    step = 0.0_wp

    vars = elem%nparam
    if (vars > MAXPARAM) then
      write (*,*) "Catastrophe in optimize_transformation_params!"
      stop
    end if

    f = 0.0d0
    cycle = 0
    hits = 0

    do
      fold = f
      f = eval_optimization_target_function(state,elem,finish)

      if (finish) exit

      if (cycle > 0) then
        if (abs(f-fold) > state%OptChangeThreshold) then
          hits = 0
        else
          hits = hits+1
        end if
        if (hits >= state%OptChangeHits) exit
      end if

      call get_params(elem,values)

      ! Calculate gradient and force constants
      do i = 1,vars
        values(i) = values(i)-state%GradientStep
        call set_params(elem,values)
        fdn = eval_optimization_target_function(state,elem)

        values(i) = values(i)+2.0d0*state%GradientStep
        call set_params(elem,values)
        fup = eval_optimization_target_function(state,elem)

        values(i) = values(i)-state%GradientStep
        grad(i) = (fup-fdn)/(2.0d0*state%GradientStep)
        force(i) = (fup+fdn-2.0d0*f)/(state%GradientStep*state%GradientStep)
      end do

      ! Quasi-Newton step
      snorm = 0.0d0
      do i = 1,vars
        if (force(i) < 0.0d0) force(i) = -force(i)
        if (force(i) < 1.0d-3) force(i) = 1.0d-3
        if (force(i) > 1.0d3) force(i) = 1.0d3
        step(i) = -grad(i)/force(i)
        snorm = snorm+step(i)*step(i)
      end do
      snorm = sqrt(snorm)

      if (snorm > state%MaxOptStep) then
        step = step*state%MaxOptStep/snorm
        snorm = state%MaxOptStep
      end if

      do while (snorm > state%MinOptStep)
        values = values+step
        call set_params(elem,values)
        fnew = eval_optimization_target_function(state,elem)

        if (fnew < f) exit

        values = values-step
        step = step/2.0d0
        call set_params(elem,values)
        snorm = snorm/2.0d0
      end do

      ! Quadratic interpolation
      if (snorm > state%MinOptStep.and.snorm < state%MaxOptStep/2.0d0) then
        values = values+step
        call set_params(elem,values)
        fnew2 = eval_optimization_target_function(state,elem)
        values = values-2.0d0*step

        a = (4.0d0*f-fnew2-3.0d0*fnew)/2.0d0
        b = (f+fnew2-2.0d0*fnew)/2.0d0

        if (b > 0.0d0) then
          x = -a/(2.0d0*b)
          if (x > 0.2d0.and.x < 1.8d0) then
            values = values+x*step
          else
            b = 0.0d0
          end if
        end if

        if (b <= 0.0d0) then
          if (fnew2 < fnew) then
            values = values+2.0d0*step
          else
            values = values+step
          end if
        end if
        call set_params(elem,values)
      end if

      cycle = cycle+1
      if (snorm <= state%MinOptStep.or.cycle >= state%MaxOptCycles) exit
    end do

    f = eval_optimization_target_function(state,elem)
    if (cycle >= state%MaxOptCycles) state%BadOptimization = 1
  end subroutine optimize_transformation_params

  !> Refine symmetry element
  function refine_symmetry_element(state,elem,build_table) result(status)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(inout) :: elem
    logical,intent(in) :: build_table
    integer :: status
    integer :: i

    status = 0

    if (build_table) then
      if (establish_pairs(state,elem) < 0) then
        state%StatPairs = state%StatPairs+1
        status = -1
        return
      end if
    end if

    ! Check for duplicates
    do i = 1,state%PlanesCount
      if (same_transform(state,state%Planes(i),elem)) then
        state%StatDups = state%StatDups+1
        status = -1
        return
      end if
    end do

    do i = 1,state%InversionCentersCount
      if (same_transform(state,state%InversionCenters(i),elem)) then
        state%StatDups = state%StatDups+1
        status = -1
        return
      end if
    end do

    do i = 1,state%NormalAxesCount
      if (same_transform(state,state%NormalAxes(i),elem)) then
        state%StatDups = state%StatDups+1
        status = -1
        return
      end if
    end do

    do i = 1,state%ImproperAxesCount
      if (same_transform(state,state%ImproperAxes(i),elem)) then
        state%StatDups = state%StatDups+1
        status = -1
        return
      end if
    end do

    if (check_transform_order(state,elem) < 0) then
      state%StatOrder = state%StatOrder+1
      status = -1
      return
    end if

    call optimize_transformation_params(state,elem)

    if (check_transform_quality(state,elem) < 0) then
      state%StatOpt = state%StatOpt+1
      status = -1
      return
    end if

    state%StatAccept = state%StatAccept+1
  end function refine_symmetry_element

  !> Initialize mirror plane
  subroutine init_mirror_plane(state,i,j,plane,success)
    type(symmetry_state_t),intent(inout) :: state
    integer,intent(in) :: i,j
    type(symmetry_element),intent(out) :: plane
    logical,intent(out) :: success
    real(wp) :: dx(DIMENSION),midpoint(DIMENSION),rab,r
    integer :: k

    success = .false.
    state%StatTotal = state%StatTotal+1

    call alloc_symmetry_element(state,plane)
    plane%transform_type = 1  ! mirror
    plane%order = 2
    plane%nparam = 4

    rab = 0.0d0
    do k = 1,DIMENSION
      dx(k) = state%Atoms(i)%x(k)-state%Atoms(j)%x(k)
      midpoint(k) = (state%Atoms(i)%x(k)+state%Atoms(j)%x(k))/2.0d0
      rab = rab+dx(k)*dx(k)
    end do
    rab = sqrt(rab)

    if (rab < state%ToleranceSame) then
      call destroy_symmetry_element(plane)
      return
    end if

    r = 0.0d0
    do k = 1,DIMENSION
      plane%normal(k) = dx(k)/rab
      r = r+midpoint(k)*plane%normal(k)
    end do

    if (r < 0.0d0) then
      r = -r
      plane%normal = -plane%normal
    end if
    plane%distance = r

    if (refine_symmetry_element(state,plane,.true.) < 0) then
      call destroy_symmetry_element(plane)
      return
    end if

    success = .true.
  end subroutine init_mirror_plane

  !> Initialize ultimate (whole-molecule) plane
  subroutine init_ultimate_plane(state,plane,success)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(out) :: plane
    logical,intent(out) :: success
    real(wp) :: d0(DIMENSION),d1(DIMENSION),d2(DIMENSION),p(DIMENSION)
    real(wp) :: r,s0,s1,s2
    real(wp),pointer :: d(:)
    integer :: i,j,k

    success = .false.
    state%StatTotal = state%StatTotal+1

    call alloc_symmetry_element(state,plane)
    plane%transform_type = 1
    plane%order = 1
    plane%nparam = 4

    d0 = 0.0d0; d1 = 0.0d0; d2 = 0.0d0
    d0(1) = 1.0d0; d1(2) = 1.0d0; d2(3) = 1.0d0

    do i = 2,state%AtomsCount
      do j = 1,i-1
        r = 0.0d0
        do k = 1,DIMENSION
          p(k) = state%Atoms(i)%x(k)-state%Atoms(j)%x(k)
          r = r+p(k)*p(k)
        end do
        r = sqrt(r)

        s0 = 0.0d0; s1 = 0.0d0; s2 = 0.0d0
        do k = 1,DIMENSION
          p(k) = p(k)/r
          s0 = s0+p(k)*d0(k)
          s1 = s1+p(k)*d1(k)
          s2 = s2+p(k)*d2(k)
        end do

        do k = 1,DIMENSION
          d0(k) = d0(k)-s0*p(k)
          d1(k) = d1(k)-s1*p(k)
          d2(k) = d2(k)-s2*p(k)
        end do
      end do
    end do

    s0 = sum(d0)
    s1 = sum(d1)
    s2 = sum(d2)

    if (s0 >= s1.and.s0 >= s2) then
      plane%normal = d0
    else if (s1 >= s0.and.s1 >= s2) then
      plane%normal = d1
    else
      plane%normal = d2
    end if

    r = sqrt(sum(plane%normal**2))
    if (r > 0.0d0) then
      plane%normal = plane%normal/r
    else
      plane%normal = [1.0d0,0.0d0,0.0d0]
    end if

    r = dot_product(state%CenterOfSomething,plane%normal)
    plane%distance = r

    do k = 1,state%AtomsCount
      plane%transform(k) = k
    end do

    if (refine_symmetry_element(state,plane,.false.) < 0) then
      call destroy_symmetry_element(plane)
      return
    end if

    success = .true.
  end subroutine init_ultimate_plane

  !> Initialize inversion center
  subroutine init_inversion_center(state,center,success)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(out) :: center
    logical,intent(out) :: success
    real(wp) :: r
    integer :: k

    success = .false.
    state%StatTotal = state%StatTotal+1

    call alloc_symmetry_element(state,center)
    center%transform_type = 2  ! invert
    center%order = 2
    center%nparam = 4

    r = sqrt(sum(state%CenterOfSomething**2))

    if (r > 0.0d0) then
      center%normal = state%CenterOfSomething/r
    else
      center%normal = [1.0d0,0.0d0,0.0d0]
    end if
    center%distance = r

    if (refine_symmetry_element(state,center,.true.) < 0) then
      call destroy_symmetry_element(center)
      return
    end if

    success = .true.
  end subroutine init_inversion_center

  !> Initialize ultimate (infinity) axis
  subroutine init_ultimate_axis(state,axis,success)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: dir(DIMENSION),rel(DIMENSION),s
    integer :: i,k

    success = .false.
    state%StatTotal = state%StatTotal+1

    call alloc_symmetry_element(state,axis)
    axis%transform_type = 3  ! rotate
    axis%order = 0
    axis%nparam = 7

    dir = 0.0d0
    do i = 1,state%AtomsCount
      s = 0.0d0
      do k = 1,DIMENSION
        rel(k) = state%Atoms(i)%x(k)-state%CenterOfSomething(k)
        s = s+rel(k)*dir(k)
      end do
      if (s >= 0.0d0) then
        dir = dir+rel
      else
        dir = dir-rel
      end if
    end do

    s = sqrt(sum(dir**2))
    if (s > 0.0d0) then
      axis%direction = dir/s
    else
      axis%direction = [1.0d0,0.0d0,0.0d0]
    end if

    s = sqrt(sum(state%CenterOfSomething**2))
    if (s > 0.0d0) then
      axis%normal = state%CenterOfSomething/s
    else
      axis%normal = [1.0d0,0.0d0,0.0d0]
    end if
    axis%distance = s

    do k = 1,state%AtomsCount
      axis%transform(k) = k
    end do

    if (refine_symmetry_element(state,axis,.false.) < 0) then
      call destroy_symmetry_element(axis)
      return
    end if

    success = .true.
  end subroutine init_ultimate_axis

  !> Initialize axis parameters from three points
  subroutine init_axis_parameters(state,a,b,c,axis,success)
    type(symmetry_state_t),intent(inout) :: state
    real(wp),intent(in) :: a(3),b(3),c(3)
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: ra,rb,rc,rab,rbc,rac,r,angle
    integer :: i,order,sign_val

    success = .false.

    ra = sqrt(sum(a**2))
    rb = sqrt(sum(b**2))
    rc = sqrt(sum(c**2))

    if (abs(ra-rb) > state%TolerancePrimary.or. &
        abs(ra-rc) > state%TolerancePrimary.or. &
        abs(rb-rc) > state%TolerancePrimary) then
      state%StatEarly = state%StatEarly+1
      return
    end if

    rab = sqrt(sum((a-b)**2))
    rac = sqrt(sum((a-c)**2))
    rbc = sqrt(sum((c-b)**2))

    if (abs(rab-rbc) > state%TolerancePrimary) then
      state%StatEarly = state%StatEarly+1
      return
    end if

    if (rab <= state%ToleranceSame.or.rbc <= state%ToleranceSame.or. &
        rac <= state%ToleranceSame) then
      state%StatEarly = state%StatEarly+1
      return
    end if

    rab = (rab+rbc)/2.0d0
    angle = PI-2.0d0*asin(rac/(2.0d0*rab))

    if (abs(angle) <= PI/(state%MaxAxisOrder+1)) then
      state%StatEarly = state%StatEarly+1
      return
    end if

    order = nint((2.0d0*PI)/angle)
    if (order <= 2.or.order > state%MaxAxisOrder) then
      state%StatEarly = state%StatEarly+1
      return
    end if

    call alloc_symmetry_element(state,axis)
    axis%order = order
    axis%nparam = 7

    r = sqrt(sum(state%CenterOfSomething**2))
    if (r > 0.0d0) then
      axis%normal = state%CenterOfSomething/r
    else
      axis%normal = [1.0d0,0.0d0,0.0d0]
    end if
    axis%distance = r

    ! Cross product for direction
    axis%direction(1) = (b(2)-a(2))*(c(3)-b(3))-(b(3)-a(3))*(c(2)-b(2))
    axis%direction(2) = (b(3)-a(3))*(c(1)-b(1))-(b(1)-a(1))*(c(3)-b(3))
    axis%direction(3) = (b(1)-a(1))*(c(2)-b(2))-(b(2)-a(2))*(c(1)-b(1))

    ! Select direction so first non-zero component is positive
    sign_val = 0
    if (axis%direction(1) < 0.0d0) then
      sign_val = 1
    else if (axis%direction(1) == 0.0d0) then
      if (axis%direction(2) < 0.0d0) then
        sign_val = 1
      else if (axis%direction(2) == 0.0d0) then
        if (axis%direction(3) < 0.0d0) sign_val = 1
      end if
    end if

    if (sign_val == 1) axis%direction = -axis%direction

    r = sqrt(sum(axis%direction**2))
    axis%direction = axis%direction/r

    success = .true.
  end subroutine init_axis_parameters

  !> Initialize C2 axis
  subroutine init_c2_axis(state,i,j,support,axis,success)
    type(symmetry_state_t),intent(inout) :: state
    integer,intent(in) :: i,j
    real(wp),intent(in) :: support(DIMENSION)
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: ris,rjs,r,center(DIMENSION)
    integer :: k

    success = .false.
    state%StatTotal = state%StatTotal+1

    ! Quick sanity check
    ris = 0.0d0
    rjs = 0.0d0
    do k = 1,DIMENSION
      ris = ris+pow2(state%Atoms(i)%x(k)-support(k))
      rjs = rjs+pow2(state%Atoms(j)%x(k)-support(k))
    end do
    ris = sqrt(ris)
    rjs = sqrt(rjs)

    if (abs(ris-rjs) > state%TolerancePrimary) then
      state%StatEarly = state%StatEarly+1
      return
    end if

    call alloc_symmetry_element(state,axis)
    axis%transform_type = 3  ! rotate
    axis%order = 2
    axis%nparam = 7

    r = sqrt(sum(state%CenterOfSomething**2))
    if (r > 0.0d0) then
      axis%normal = state%CenterOfSomething/r
    else
      axis%normal = [1.0d0,0.0d0,0.0d0]
    end if
    axis%distance = r

    r = 0.0d0
    do k = 1,DIMENSION
      center(k) = (state%Atoms(i)%x(k)+state%Atoms(j)%x(k))/2.0d0-support(k)
      r = r+center(k)*center(k)
    end do
    r = sqrt(r)

    if (r <= state%TolerancePrimary) then
      ! C2 is underdefined
      if (state%MolecularPlaneExists) then
        axis%direction = state%MolecularPlane%normal
      else
        do k = 1,DIMENSION
          center(k) = state%Atoms(i)%x(k)-state%Atoms(j)%x(k)
        end do
        if (abs(center(3))+abs(center(2)) > state%ToleranceSame) then
          axis%direction = [0.0d0,center(3),-center(2)]
        else
          axis%direction = [-center(3),0.0d0,center(1)]
        end if
        r = sqrt(sum(axis%direction**2))
        axis%direction = axis%direction/r
      end if
    else
      axis%direction = center/r
    end if

    if (refine_symmetry_element(state,axis,.true.) < 0) then
      call destroy_symmetry_element(axis)
      return
    end if

    success = .true.
  end subroutine init_c2_axis

  !> Initialize higher-order axis
  subroutine init_higher_axis(state,ia,ib,ic,axis,success)
    type(symmetry_state_t),intent(inout) :: state
    integer,intent(in) :: ia,ib,ic
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: a(DIMENSION),b(DIMENSION),c(DIMENSION)
    integer :: i

    success = .false.
    state%StatTotal = state%StatTotal+1

    do i = 1,DIMENSION
      a(i) = state%Atoms(ia)%x(i)-state%CenterOfSomething(i)
      b(i) = state%Atoms(ib)%x(i)-state%CenterOfSomething(i)
      c(i) = state%Atoms(ic)%x(i)-state%CenterOfSomething(i)
    end do

    call init_axis_parameters(state,a,b,c,axis,success)
    if (.not.success) return

    axis%transform_type = 3  ! rotate

    if (refine_symmetry_element(state,axis,.true.) < 0) then
      call destroy_symmetry_element(axis)
      success = .false.
      return
    end if

    success = .true.
  end subroutine init_higher_axis

  !> Initialize improper axis
  subroutine init_improper_axis(state,ia,ib,ic,axis,success)
    type(symmetry_state_t),intent(inout) :: state
    integer,intent(in) :: ia,ib,ic
    type(symmetry_element),intent(out) :: axis
    logical,intent(out) :: success
    real(wp) :: a(DIMENSION),b(DIMENSION),c(DIMENSION)
    real(wp) :: centerpoint(DIMENSION),r
    integer :: i

    success = .false.
    state%StatTotal = state%StatTotal+1

    do i = 1,DIMENSION
      a(i) = state%Atoms(ia)%x(i)-state%CenterOfSomething(i)
      b(i) = state%Atoms(ib)%x(i)-state%CenterOfSomething(i)
      c(i) = state%Atoms(ic)%x(i)-state%CenterOfSomething(i)
    end do

    r = 0.0d0
    do i = 1,DIMENSION
      centerpoint(i) = a(i)+c(i)+2.0d0*b(i)
      r = r+centerpoint(i)*centerpoint(i)
    end do
    r = sqrt(r)

    if (r <= state%ToleranceSame) then
      state%StatEarly = state%StatEarly+1
      return
    end if

    centerpoint = centerpoint/r
    r = dot_product(centerpoint,b)
    b = 2.0d0*r*centerpoint-b

    call init_axis_parameters(state,a,b,c,axis,success)
    if (.not.success) return

    axis%transform_type = 4  ! rotate_reflect

    if (refine_symmetry_element(state,axis,.true.) < 0) then
      call destroy_symmetry_element(axis)
      success = .false.
      return
    end if

    success = .true.
  end subroutine init_improper_axis

  !> Find center of something (centroid)
  subroutine find_center_of_something(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i,j
    real(wp) :: coord_sum(DIMENSION),r

    coord_sum = 0.0d0
    do i = 1,state%AtomsCount
      coord_sum = coord_sum+state%Atoms(i)%x
    end do
    state%CenterOfSomething = coord_sum/dble(state%AtomsCount)

    if (allocated(state%DistanceFromCenter)) deallocate (state%DistanceFromCenter)
    allocate (state%DistanceFromCenter(state%AtomsCount))

    do i = 1,state%AtomsCount
      r = 0.0d0
      do j = 1,DIMENSION
        r = r+pow2(state%Atoms(i)%x(j)-state%CenterOfSomething(j))
      end do
      state%DistanceFromCenter(i) = r
    end do
  end subroutine find_center_of_something

  !> Add plane to planes array
  subroutine add_plane(state,plane)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(in) :: plane
    type(symmetry_element),allocatable :: temp(:)

    state%PlanesCount = state%PlanesCount+1
    if (allocated(state%Planes)) then
      allocate (temp(state%PlanesCount))
      temp(1:state%PlanesCount-1) = state%Planes
      temp(state%PlanesCount) = plane
      call move_alloc(temp,state%Planes)
    else
      allocate (state%Planes(1))
      state%Planes(1) = plane
    end if
  end subroutine add_plane

  !> Add normal axis to array
  subroutine add_normal_axis(state,axis)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(in) :: axis
    type(symmetry_element),allocatable :: temp(:)

    state%NormalAxesCount = state%NormalAxesCount+1
    if (allocated(state%NormalAxes)) then
      allocate (temp(state%NormalAxesCount))
      temp(1:state%NormalAxesCount-1) = state%NormalAxes
      temp(state%NormalAxesCount) = axis
      call move_alloc(temp,state%NormalAxes)
    else
      allocate (state%NormalAxes(1))
      state%NormalAxes(1) = axis
    end if
  end subroutine add_normal_axis

  !> Add improper axis to array
  subroutine add_improper_axis(state,axis)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element),intent(in) :: axis
    type(symmetry_element),allocatable :: temp(:)

    state%ImproperAxesCount = state%ImproperAxesCount+1
    if (allocated(state%ImproperAxes)) then
      allocate (temp(state%ImproperAxesCount))
      temp(1:state%ImproperAxesCount-1) = state%ImproperAxes
      temp(state%ImproperAxesCount) = axis
      call move_alloc(temp,state%ImproperAxes)
    else
      allocate (state%ImproperAxes(1))
      state%ImproperAxes(1) = axis
    end if
  end subroutine add_improper_axis

  !> Find planes of symmetry
  subroutine find_planes(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i,j
    type(symmetry_element) :: plane
    logical :: success

    call init_ultimate_plane(state,plane,success)
    if (success) then
      state%MolecularPlane = plane
      state%MolecularPlaneExists = .true.
      call add_plane(state,plane)
    end if

    do i = 2,state%AtomsCount
      do j = 1,i-1
        if (state%Atoms(i)%atom_type /= state%Atoms(j)%atom_type) cycle

        call init_mirror_plane(state,i,j,plane,success)
        if (success) call add_plane(state,plane)
      end do
    end do
  end subroutine find_planes

  !> Find inversion centers
  subroutine find_inversion_centers(state)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element) :: center
    logical :: success

    call init_inversion_center(state,center,success)
    if (success) then
      state%InversionCentersCount = 1
      allocate (state%InversionCenters(1))
      state%InversionCenters(1) = center
    end if
  end subroutine find_inversion_centers

  !> Find infinity axis
  subroutine find_infinity_axis(state)
    type(symmetry_state_t),intent(inout) :: state
    type(symmetry_element) :: axis
    logical :: success

    call init_ultimate_axis(state,axis,success)
    if (success) call add_normal_axis(state,axis)
  end subroutine find_infinity_axis

  !> Find C2 axes
  subroutine find_c2_axes(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i,j,k,l,m
    real(wp) :: center(DIMENSION),r
    real(wp),allocatable :: distances(:)
    type(symmetry_element) :: axis
    logical :: success

    allocate (distances(state%AtomsCount))

    do i = 2,state%AtomsCount
      do j = 1,i-1
        if (state%Atoms(i)%atom_type /= state%Atoms(j)%atom_type) cycle
        if (abs(state%DistanceFromCenter(i)-state%DistanceFromCenter(j)) > &
            state%TolerancePrimary) cycle

        ! Try using CenterOfSomething
        r = 0.0d0
        do k = 1,DIMENSION
          center(k) = (state%Atoms(i)%x(k)+state%Atoms(j)%x(k))/2.0d0
          r = r+pow2(center(k)-state%CenterOfSomething(k))
        end do
        r = sqrt(r)

        if (r > 5.0d0*state%TolerancePrimary) then
          call init_c2_axis(state,i,j,state%CenterOfSomething,axis,success)
          if (success) call add_normal_axis(state,axis)
          cycle
        end if

        ! Try through atoms
        do k = 1,state%AtomsCount
          call init_c2_axis(state,i,j,state%Atoms(k)%x,axis,success)
          if (success) call add_normal_axis(state,axis)
        end do

        ! Calculate distances for prescreening
        do k = 1,state%AtomsCount
          r = 0.0d0
          do l = 1,DIMENSION
            r = r+pow2(state%Atoms(k)%x(l)-center(l))
          end do
          distances(k) = sqrt(r)
        end do

        ! Try through midpoints of atom pairs
        do k = 1,state%AtomsCount
          do l = 1,state%AtomsCount
            if (state%Atoms(k)%atom_type /= state%Atoms(l)%atom_type) cycle
            if (abs(state%DistanceFromCenter(k)-state%DistanceFromCenter(l)) > &
                state%TolerancePrimary.or. &
                abs(distances(k)-distances(l)) > state%TolerancePrimary) cycle

            do m = 1,DIMENSION
              center(m) = (state%Atoms(k)%x(m)+state%Atoms(l)%x(m))/2.0d0
            end do

            call init_c2_axis(state,i,j,center,axis,success)
            if (success) call add_normal_axis(state,axis)
          end do
        end do
      end do
    end do

    deallocate (distances)
  end subroutine find_c2_axes

  !> Find higher-order axes
  subroutine find_higher_axes(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i,j,k
    type(symmetry_element) :: axis
    logical :: success

    do i = 1,state%AtomsCount
      do j = i+1,state%AtomsCount
        if (state%Atoms(i)%atom_type /= state%Atoms(j)%atom_type) cycle
        if (abs(state%DistanceFromCenter(i)-state%DistanceFromCenter(j)) > &
            state%TolerancePrimary) cycle

        do k = 1,state%AtomsCount
          if (state%Atoms(i)%atom_type /= state%Atoms(k)%atom_type) cycle
          if (abs(state%DistanceFromCenter(i)-state%DistanceFromCenter(k)) > &
              state%TolerancePrimary.or. &
              abs(state%DistanceFromCenter(j)-state%DistanceFromCenter(k)) > &
              state%TolerancePrimary) cycle

          call init_higher_axis(state,i,j,k,axis,success)
          if (success) call add_normal_axis(state,axis)
        end do
      end do
    end do
  end subroutine find_higher_axes

  !> Find improper axes
  subroutine find_improper_axes(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i,j,k
    type(symmetry_element) :: axis
    logical :: success

    do i = 1,state%AtomsCount
      do j = i+1,state%AtomsCount
        do k = 1,state%AtomsCount
          call init_improper_axis(state,i,j,k,axis,success)
          if (success) call add_improper_axis(state,axis)
        end do
      end do
    end do
  end subroutine find_improper_axes

  !> Find all symmetry elements
  subroutine find_symmetry_elements(state)
    type(symmetry_state_t),intent(inout) :: state
    call find_center_of_something(state)
    call find_inversion_centers(state)
    call find_planes(state)
    call find_infinity_axis(state)
    call find_c2_axes(state)
    call find_higher_axes(state)
    call find_improper_axes(state)
  end subroutine find_symmetry_elements

  !> Compare axes for sorting
  function compare_axes(a,b) result(cmp)
    type(symmetry_element),intent(in) :: a,b
    integer :: cmp
    integer :: order_a,order_b

    order_a = a%order
    order_b = b%order
    if (order_a == 0) order_a = 10000
    if (order_b == 0) order_b = 10000

    cmp = order_b-order_a
    if (cmp /= 0) return

    if (a%maxdev > b%maxdev) then
      cmp = -1
    else if (a%maxdev < b%maxdev) then
      cmp = 1
    else
      cmp = 0
    end if
  end function compare_axes

  !> Sort symmetry elements (simple bubble sort)
  subroutine sort_symmetry_elements(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i,j
    type(symmetry_element) :: temp

    ! Sort planes
    do i = 1,state%PlanesCount-1
      do j = i+1,state%PlanesCount
        if (compare_axes(state%Planes(i),state%Planes(j)) < 0) then
          temp = state%Planes(i)
          state%Planes(i) = state%Planes(j)
          state%Planes(j) = temp
        end if
      end do
    end do

    ! Sort normal axes
    do i = 1,state%NormalAxesCount-1
      do j = i+1,state%NormalAxesCount
        if (compare_axes(state%NormalAxes(i),state%NormalAxes(j)) < 0) then
          temp = state%NormalAxes(i)
          state%NormalAxes(i) = state%NormalAxes(j)
          state%NormalAxes(j) = temp
        end if
      end do
    end do

    ! Sort improper axes
    do i = 1,state%ImproperAxesCount-1
      do j = i+1,state%ImproperAxesCount
        if (compare_axes(state%ImproperAxes(i),state%ImproperAxes(j)) < 0) then
          temp = state%ImproperAxes(i)
          state%ImproperAxes(i) = state%ImproperAxes(j)
          state%ImproperAxes(j) = temp
        end if
      end do
    end do
  end subroutine sort_symmetry_elements

  !> Summarize symmetry elements
  subroutine summarize_symmetry_elements(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i

    if (allocated(state%NormalAxesCounts)) deallocate (state%NormalAxesCounts)
    if (allocated(state%ImproperAxesCounts)) deallocate (state%ImproperAxesCounts)

    allocate (state%NormalAxesCounts(0:state%MaxAxisOrder))
    allocate (state%ImproperAxesCounts(0:state%MaxAxisOrder))

    state%NormalAxesCounts = 0
    state%ImproperAxesCounts = 0

    do i = 1,state%NormalAxesCount
      state%NormalAxesCounts(state%NormalAxes(i)%order) = &
        state%NormalAxesCounts(state%NormalAxes(i)%order)+1
    end do

    do i = 1,state%ImproperAxesCount
      state%ImproperAxesCounts(state%ImproperAxes(i)%order) = &
        state%ImproperAxesCounts(state%ImproperAxes(i)%order)+1
    end do
  end subroutine summarize_symmetry_elements

  !> Report symmetry elements brief
  subroutine report_symmetry_elements_brief(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i,n,tlen
    character(len=32) :: buf

    state%SymmetryCode = ""
    n = 0

    if (state%PlanesCount+state%NormalAxesCount+state%ImproperAxesCount+ &
        state%InversionCentersCount > 0) then
      if (state%InversionCentersCount > 0) then
        state%SymmetryCode(n+1:n+4) = "(i) "
        n = n+4
      end if

      if (state%NormalAxesCounts(0) == 1) then
        state%SymmetryCode(n+1:n+7) = "(Cinf) "
        n = n+7
      else if (state%NormalAxesCounts(0) > 1) then
        write (buf,'(I0,A)') state%NormalAxesCounts(0),"*(Cinf) "
        tlen = len_trim(buf)+1
        state%SymmetryCode(n+1:n+tlen) = buf(1:tlen)
        n = n+tlen
      end if

      do i = state%MaxAxisOrder,2,-1
        if (state%NormalAxesCounts(i) == 1) then
          write (buf,'(A,I0,A)') "(C",i,") "
          tlen = len_trim(buf)+1
          state%SymmetryCode(n+1:n+tlen) = buf(1:tlen)
          n = n+tlen
        else if (state%NormalAxesCounts(i) > 1) then
          write (buf,'(I0,A,I0,A)') state%NormalAxesCounts(i),"*(C",i,") "
          tlen = len_trim(buf)+1
          state%SymmetryCode(n+1:n+tlen) = buf(1:tlen)
          n = n+tlen
        end if
      end do

      do i = state%MaxAxisOrder,2,-1
        if (state%ImproperAxesCounts(i) == 1) then
          write (buf,'(A,I0,A)') "(S",i,") "
          tlen = len_trim(buf)+1
          state%SymmetryCode(n+1:n+tlen) = buf(1:tlen)
          n = n+tlen
        else if (state%ImproperAxesCounts(i) > 1) then
          write (buf,'(I0,A,I0,A)') state%ImproperAxesCounts(i),"*(S",i,") "
          tlen = len_trim(buf)+1
          state%SymmetryCode(n+1:n+tlen) = buf(1:tlen)
          n = n+tlen
        end if
      end do

      if (state%PlanesCount == 1) then
        state%SymmetryCode(n+1:n+8) = "(sigma) "
        n = n+8
      else if (state%PlanesCount > 1) then
        write (buf,'(I0,A)') state%PlanesCount,"*(sigma) "
        tlen = len_trim(buf)+1
        state%SymmetryCode(n+1:n+tlen) = buf(1:tlen)
        n = n+tlen
      end if
    end if
  end subroutine report_symmetry_elements_brief

  !> Report highest rotation axis only
  subroutine report_symmetry_elements_brief_conly(state)
    type(symmetry_state_t),intent(inout) :: state
    integer :: i
    character(len=8) :: buf

    state%MaxRotAxis = ""

    if (state%PlanesCount+state%NormalAxesCount+state%ImproperAxesCount+ &
        state%InversionCentersCount > 0) then
      do i = state%MaxAxisOrder,2,-1
        if (state%NormalAxesCounts(i) >= 1) then
          write (buf,'(A,I0)') "C",i
          state%MaxRotAxis = trim(buf)
          return
        end if
      end do
    end if
  end subroutine report_symmetry_elements_brief_conly

  !> Identify point group
  function identify_point_group(state) result(last_matching)
    type(symmetry_state_t),intent(in) :: state
    integer :: last_matching
    integer :: i,matching_count

    last_matching = -1
    matching_count = 0

    do i = 1,PointGroupsCount
      if (len_trim(PointGroups(i)%group_name) == 0) cycle
      if (trim(state%SymmetryCode) == trim(PointGroups(i)%symmetry_code)) then
        last_matching = i
        matching_count = matching_count+1
      end if
    end do

    if (matching_count == 0) then
      last_matching = -1
    else if (matching_count > 1) then
      last_matching = -1
    end if
  end function identify_point_group

  !> Main entry point: determine Schoenflies symbol
  subroutine schoenflies(natoms,attype,coord,symbol,paramar)
    integer,intent(in)  :: natoms
    integer,intent(in)  :: attype(natoms)
    real(wp),intent(in) :: coord(3,natoms)
    character(len=*),intent(out) :: symbol
    real(wp),intent(in),optional :: paramar(11)
    type(symmetry_state_t) :: state
    integer :: last_pg,i

    call init_symmetry_state(state)

    ! Set parameters if provided
    if (present(paramar)) then
      state%verbose = nint(paramar(1))
      state%MaxAxisOrder = nint(paramar(2))
      state%MaxOptCycles = nint(paramar(3))
      state%ToleranceSame = paramar(4)
      state%TolerancePrimary = paramar(5)
      state%ToleranceFinal = paramar(6)
      state%MaxOptStep = paramar(7)
      state%MinOptStep = paramar(8)
      state%GradientStep = paramar(9)
      state%OptChangeThreshold = paramar(10)
      state%OptChangeHits = nint(paramar(11))
    end if

    ! Set up atoms
    state%AtomsCount = natoms
    allocate (state%Atoms(state%AtomsCount))

    do i = 1,state%AtomsCount
      state%Atoms(i)%atom_type = attype(i)
      state%Atoms(i)%x(1) = coord(1,i)
      state%Atoms(i)%x(2) = coord(2,i)
      state%Atoms(i)%x(3) = coord(3,i)
    end do

    ! Find and analyze symmetry
    call find_symmetry_elements(state)
    call sort_symmetry_elements(state)
    call summarize_symmetry_elements(state)
    call report_symmetry_elements_brief(state)

    last_pg = identify_point_group(state)

    if (last_pg >= 1) then
      symbol = trim(PointGroups(last_pg)%group_name)
    else
      call report_symmetry_elements_brief_conly(state)
      if (len_trim(state%MaxRotAxis) == 0) then
        symbol = "C1"
      else
        symbol = trim(state%MaxRotAxis)
      end if
    end if
  end subroutine schoenflies

! ══════════════════════════════════════════════════════════════════════════════

  subroutine getsym(pr,iunit,n,iat,xyz,sfsym,symthr,maxatdesy)
    use iso_fortran_env,only:wp => real64
    implicit none
    logical,intent(in)           :: pr
    integer,intent(in)           :: iunit
    integer,intent(in)           :: n
    integer,intent(in)           :: iat(n)
    real(wp),intent(in)          :: xyz(3,n)
    character(len=*),intent(out) :: sfsym
    real(wp),intent(in),optional :: symthr    ! default 0.1
    integer,intent(in),optional  :: maxatdesy ! default 200
    real(wp) :: thr
    integer  :: maxat
    character(len=8) :: atmp
    real(wp) :: paramar(11)

    thr = 0.1_wp; if (present(symthr)) thr = symthr
    maxat = 200; if (present(maxatdesy)) maxat = maxatdesy

    if (n > maxat) then
      if (pr) write (iunit,*) 'symmetry recognition skipped because # atoms >',maxat
      sfsym = 'none'
      return
    end if

    if (pr) write (iunit,'(a)')
    paramar(1) = -1       ! verbose
    paramar(2) = 10       ! MaxAxisOrder
    paramar(3) = 100      ! MaxOptCycles
    paramar(4) = 0.001d0  ! ToleranceSame
    paramar(5) = 0.5d0    ! TolerancePrimary
    paramar(6) = thr      ! ToleranceFinal
    paramar(7) = 0.5d0    ! MaxOptStep
    paramar(8) = 1.0d-7   ! MinOptStep
    paramar(9) = 1.0d-7   ! GradientStep
    paramar(10) = 1.0d-8   ! OptChangeThreshold
    paramar(11) = 5        ! OptChangeHits

    atmp = '        '
    call schoenflies(n,iat,xyz,atmp,paramar)

    sfsym(1:3) = atmp(1:3)
    if (sfsym(1:1) == 'D') sfsym(1:1) = 'd'
    if (sfsym(1:1) == 'C') sfsym(1:1) = 'c'
    if (sfsym(1:1) == 'T') sfsym(1:1) = 't'
    if (sfsym(1:1) == 'O') sfsym(1:1) = 'o'
    if (sfsym(1:1) == 'I') sfsym(1:1) = 'i'
    if (sfsym(1:1) == 'S') sfsym(1:1) = 's'
    ! Linear molecules: fix to correct 3-char codes
    if (sfsym(1:3) == 'dih') sfsym(1:3) = 'din'
    if (sfsym(1:3) == 'civ') sfsym(1:3) = 'cin'
    if (sfsym(3:3) > 'v'.or.sfsym(3:3) < 'a') sfsym(3:3) = ' '

    if (pr) write (iunit,'(a3,'' symmetry found (for desy threshold: '',e9.2,'')'')') &
      sfsym,thr
  end subroutine getsym

! ══════════════════════════════════════════════════════════════════════════════
end module symmetry_i
