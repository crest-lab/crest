!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2025 Philipp Pracht
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

!> This module implemnts a simple L-BFGS with different coordinate choices

module lbfgs_module
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use crest_calculator
  use strucrd
  use optimize_type  !> This module provides the 'optimizer' type.
  use optimize_utils
  use coordinate_transform_module
  implicit none
  private

  public :: lbfgs_optimize

!========================================================================================!
!========================================================================================!
contains !>  MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  function lbfgs_direction(nvar,g,k,OPT,gamm) result(d)
    !*******************************************************************************
    !* Two-loop recursion routine to compute the search direction.
    !*
    !* This function uses the stored correction pairs (S and Y) and
    !* corresponding scaling factors (rho) to approximate the product
    !* H*g, where H is the inverse Hessian approximation.
    !*
    !* The algorithm proceeds in two loops:
    !* 1. The first (backward) loop computes the coefficients "alpha" and
    !*    subtracts corrections from the gradient.
    !* 2. The second (forward) loop applies the corrections in the reverse order.
    !* Finally, the result is negated to obtain the descent direction.
    !* Note, this routine is NOT called for the very first iteration (k == 0)
    !*
    !* @param nvar   Dimension of the variable space.
    !* @param k      Number of stored correction pairs (k ≤ m).
    !* @param g      Current gradient vector.
    !* @param OPT    optimizer type that stores the following variables
    !*     => S      Matrix containing the s-vectors (x_{k+1} - x_k),
    !*               of size(nvar, m).
    !*     => Y      Matrix containing the y-vectors (g_{k+1} - g_k),
    !*               of size(nvar, m).
    !*     => rho    Array of stored values 1/(y^T*s) for each correction.
    !*     => alpha  coefficients (get computed in this function)
    !*     => q      temporary workspace
    !* @param gamm  Scaling factor for the initial Hessian approximation.
    !*
    !* @return d     Computed search direction (negative approximate inverse
    !*               Hessian times g).
    !********************************************************************************
    !> INPUT
    integer,intent(in) :: nvar,k
    type(optimizer),intent(inout) :: OPT
    real(wp),intent(in) :: g(nvar)
    real(wp),intent(in) :: gamm
    !> OUTPUT
    real(wp) :: d(nvar)
    !> LOCAL
    integer :: i

    associate (S => OPT%S,Y => OPT%Y,alpha => OPT%alpha,rho => OPT%rho,q => OPT%q)

      !write(*,*) k
      !write(*,*) S(:,k)
      !write(*,*) Y(:,k)
      !> Initialize q with the current gradient.
      q = g

      !---------------------------------------------------------
      !> First loop (backward pass): for i = k downto 1,
      !> compute the coefficient alpha(i) and update q.
      !---------------------------------------------------------
      do i = k,1,-1
        alpha(i) = rho(i)*dot_product(S(1:nvar,i),q)
        q = q-alpha(i)*Y(1:nvar,i)
      end do

      !---------------------------------------------------------
      !> Apply the initial Hessian approximation.
      !> We use a scaled identity matrix H0 = gamm * I.
      !---------------------------------------------------------
      d = gamm*q

      !---------------------------------------------------------
      !> Second loop (forward pass): for i = 1 to k,
      !> compute the correction and update d.
      !---------------------------------------------------------
      do i = 1,k
        d = d+S(1:nvar,i)*(alpha(i)-rho(i)*dot_product(Y(1:nvar,i),d))
      end do

      !> The final search direction is the negative of d.
      d = -d

    end associate
  end function lbfgs_direction

!========================================================================================!

  subroutine lbfgs_optimize(mol,calc,etot,grd,pr,wr,io,logfile)
    !**************************************************************************
    !* L-BFGS Optimization Routine
    !*
    !* Performs optimization using the Limited-memory BFGS (L-BFGS) algorithm.
    !* The routine updates the coordinate vector x to approach a local minimum of the
    !* objective function. It integrates with an optimizer type (OPT) to manage the
    !* correction pairs (s and y) and related internal data using associate constructs.
    !*
    !* The main steps include:
    !*  1. Evaluating the objective function and gradient at the current x.
    !*  2. Computing the search direction via the two-loop recursion (lbfgs_direction).
    !*  3. Updating x using a fixed step (with the option to incorporate a line search).
    !*  4. Updating the correction pairs: s = x_new - x and y = g_new - g, while managing
    !*     the history using a shifting strategy when full.
    !*
    !* @param io       Integer. Output status variable (0 indicates success).
    !* @param logfile  name of the trajectory logfile (default crestopt.log.xyz).
    !**************************************************************************
    implicit none
    !> INPUT
    type(coord),intent(inout) :: mol
    type(calcdata),intent(in) :: calc
    real(wp),intent(inout) :: etot
    real(wp),intent(inout) :: grd(3,mol%nat)
    logical,intent(in)     :: pr
    logical,intent(in)     :: wr   !> write the trajectory logfile?
    character(len=*),intent(in) :: logfile
    !> OUTPUT
    integer,intent(out)    :: io
    !> LOCAL
    type(optimizer) :: OPT
    integer :: iter,k,nvar,m
    integer :: tight,max_iter
    real(wp) :: gnorm,deltaE,energy
    real(wp) :: ethr,gthr,maxerise
    logical :: econverged,gconverged,converged,Erise
    real(wp),allocatable :: x(:),g(:),d(:),g_new(:),x_new(:),gtmp(:,:)
    real(wp) :: f,f_new,gamm,step
    integer :: ilog

    !> Prepare settings
    io = 0
    nvar = compute_nvar(mol)
    m = calc%lbfgs_histsize
    gnorm = norm2(grd)
    deltaE = huge(deltaE)
    tight = calc%optlev
    call get_optthr(mol%nat,tight,calc,ethr,gthr)
    max_iter = calc%maxcycle !> automatic setting in get_optthr or by user
    maxerise = calc%maxerise
    econverged = .false.
    gconverged = .false.
    converged = .false.
    mol%wrextxyz = calc%logextxyz

!>--- only write the trajectory log when requested (wr). The logfile name is
!>    fixed, so an unconditional open/write races across threads in the parallel
!>    ensemble optimization (where wr=.false.). Mirrors the ancopt behaviour.
    if (wr) then
      open (newunit=ilog,file=logfile)
      if (calc%logextxyz) then
        mol%gradient = grd
      end if
      call mol%appendlog(ilog,etot)
    end if

    !$omp critical
    !> Allocate the vectors for position, gradient, and search direction.
    allocate (x(nvar),g(nvar),d(nvar),g_new(nvar),x_new(nvar),gtmp(3,mol%nat))
    !> Allocate matrices to store up to m correction pairs (columns correspond to each stored pair).
    call OPT%allocatelbfgs(nvar,m)
    !$omp end critical

    associate (S => OPT%S,Y => OPT%Y,rho => OPT%rho,alpha => OPT%alpha)
      S = 0.0_wp
      Y = 0.0_wp
      rho = 0.0_wp
      k = 0  ! Initially, no correction pairs are stored.

      !> First trafo
      call transform_mol('cart2v',mol,nvar,x)
      call transform_grd('cart2v',mol,grd,nvar,g)

      iter = 0
      if (pr) then
        call print_optiter(iter)
        write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') etot
        write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') 0.0_wp
        write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")') gnorm
      end if

      LBFGS_iter: do while (.not.converged.and.iter < max_iter)
        iter = iter+1
        if (pr) call print_optiter(iter)

        if (iter == 1) then
          !> First iteration: use the steepest descent direction.
          d = -g
        else
          !---------------------------------------------------------
          !> Determine the scaling factor gamm for the initial Hessian.
          !  Here we use the most recent correction pair.
          !---------------------------------------------------------
          if (k > 0) then
            gamm = dot_product(S(1:nvar,k),Y(1:nvar,k))/ &
                   dot_product(Y(1:nvar,k),Y(1:nvar,k))
          else
            gamm = 1.0_wp
          end if

          !> Compute the search direction using the two-loop recursion.
          d = lbfgs_direction(nvar,g,k,OPT,gamm)
        end if

        !---------------------------------------------------------
        !> A fixed step size could be used here for simplicity.
        !  In a full implementation, a line search could be used.
        !  If the energy rises, we reduce the stepsize iteratively
        !---------------------------------------------------------
        step = 0.2_wp

        Erise = .true.
        do while (Erise)
          !> Update the position: x_new = x + step * d.
          x_new = x+step*d

          !====================================================================!
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
          !> Evaluate the objective and gradient at the new position.
          call transform_mol('v2cart',mol,nvar,x_new)
          grd = 0.0_wp
          call engrad(mol,calc,energy,gtmp,io)
          if (wr) call mol%appendlog(ilog,energy)
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
          !====================================================================!
          if (io /= 0) exit

          call transform_grd('cart2v',mol,gtmp,nvar,g_new)
          gnorm = norm2(gtmp)
          deltaE = energy-etot

          Erise = (deltaE > maxerise)
          if (Erise) then
            step = step*0.25_wp
            if (pr) then
              write (*,'(" * energy rise detected, decreasing stepsize")')
            end if
          end if
        end do
        econverged = abs(deltaE) .lt. ethr
        gconverged = gnorm .lt. gthr

        !---------------------------------------------------------
        !> Compute the correction pair:
        !   s = x_new - x and y = g_new - g.
        !---------------------------------------------------------
        if (k < m) then
          !> If there is still room in the history, simply add the new pair.
          k = k+1
          S(1:nvar,k) = x_new-x
          Y(1:nvar,k) = g_new-g
        else
          !> When the history is full, shift the stored corrections and
          !> insert the new pair at the end.
          S(1:nvar,1:m-1) = S(1:nvar,2:m)
          Y(1:nvar,1:m-1) = Y(1:nvar,2:m)
          S(1:nvar,m) = x_new-x
          Y(1:nvar,m) = g_new-g
        end if

        !> Update the scaling factor for the new correction pair.
        rho(k) = 1.0_wp/dot_product(Y(1:nvar,k),S(1:nvar,k))

        !> Update the current position, gradient, and function value.
        x = x_new
        g = g_new
        etot = energy

        !> Optional: print iteration information.
        if (pr) then
          write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') etot
          write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') deltaE
          write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
          call print_convd(econverged,gconverged)
        end if
        converged = econverged.and.gconverged
      end do LBFGS_iter

      !> stop associating
    end associate

    !> Final trafo
    call transform_mol('v2cart',mol,nvar,x_new)
    call transform_grd('v2cart',mol,grd,nvar,g_new)

    !> close the trajectory logfile (if it was opened)
    if (wr) close (ilog)

    !> Deallocate all temporary arrays.
    deallocate (x_new,g_new,d,g,x)
  end subroutine lbfgs_optimize

!========================================================================================!
!========================================================================================!
end module lbfgs_module

