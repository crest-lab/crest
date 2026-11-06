!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Lukas Rindt
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

!> This module implements a standard NR algorithm (in Cart. coords)

module newton_raphson_module
  use iso_fortran_env,only:wp => real64,sp => real32
  use crest_calculator
  use axis_module
  use strucrd
  use irmsd_module,only:rmsd,rmsd_core_cache

  use optimize_type
  use optimize_maths
  use modelhessian_module
  use hessupdate_module
  use optimize_utils
  use hessian_reconstruct
  use hr_utils
  implicit none
  private

  public :: newton_raphson

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine newton_raphson(mol,calc,etot,grd,pr,wr,iostatus,logfile)
!*************************************************************************
!> subroutine rfopt
!> Implementation of the standard rational function optimizer (RFO)
!>
!> Input/Output:
!>      mol  - object containing the molecule,
!>             Cartesian coordinates in Bohrs.
!>             will be overwritten on output
!>     calc  - object containing calculation settings
!>             and optimization thresholds (look for calc% )
!>     etot  - on input initial energy (do a singlepoint before ancopt)
!>             on output final energy
!>      grd  - Cartesian gradient
!>       pr  - printout bool
!>       wr  - logfile bool (name given by logfile, default crestopt.log.xyz)
!>  iostatus - return status of the routine
!>             (success=0, error<0, not converged>0)
!>   logfile - name of the trajectory logfile to write when wr is set
!!***********************************************************************
    implicit none
    !> INPUT/OUTPUT
    type(coord),intent(inout) :: mol
    type(calcdata),intent(inout) :: calc
    real(wp),intent(inout) :: etot
    real(wp),intent(inout) :: grd(3,mol%nat)
    logical,intent(in) :: pr
    logical,intent(in) :: wr
    integer,intent(out) :: iostatus
    character(len=*),intent(in) :: logfile
    !> LOCAL
    integer  :: tight
    real(wp) :: eel
    real(wp) :: et
    real(wp) :: egap
    logical :: fail
    !> Local objects
    type(coord)   :: molopt
    type(optimizer)  :: OPT
    type(mhparam) :: mhset

    real(wp) :: step,amu2au,au2cm,dumi,dumj,damp,hlow,edum,s6,thr
    real(wp) :: maxdispl,gthr,ethr,hmax,energy,rij(3),t1,t0,w1,w0
    real(wp) :: rot(3),gnorm
    integer :: n3,i,j,k,l,jjj,ic,jc,ia,ja,ii,jj,info,nat3
    integer :: nvar,iter,nread,maxcycle,maxmicro,itry,maxopt,iupdat,iii
    integer :: id,ihess,error
    integer :: ilog,imax(3)
    real(wp) :: depred,echng,alp,alpold,gnold,eold,gchng,dummy,dsnrm,maxd
    real(wp),allocatable :: h(:,:)
    real(wp),allocatable :: b(:,:)
    real(wp),allocatable :: fc(:)
    real(wp),allocatable :: eig(:)
    real(wp),allocatable :: aux(:)
    real(wp),allocatable :: hess(:)
    integer,allocatable :: iwork(:)
    integer,allocatable :: totsym(:)
    real(wp),allocatable :: pmode(:,:)
    real(wp),allocatable :: grd1(:)
    real(wp),allocatable :: gold(:)
    real(wp),allocatable :: displ(:)
    integer :: nvar1,npvar,npvar1
    real(wp), allocatable :: int_hess(:), c(:)
    integer, allocatable :: IPIV(:)
    integer :: info2,info3
    type(convergence_log),allocatable :: avconv
    real(wp) :: rmsdval
    type(rmsd_core_cache) :: rcache
    integer :: modef
    logical :: ex,converged,linear,exact
    logical :: econverged,gconverged,lowered
    real(wp) :: estart,esave
    real(wp),parameter :: r4dum = 1.e-8
    integer :: unit
    real(wp), allocatable :: dx_test(:)
    !> LAPACK & BLAS
    external :: dgemv
    external :: dppsv
    external :: dspsv
    real(wp), external :: ddot
    !real(sp),external :: sdot

    iostatus = 0
    fail = .false.
    converged = .false.
    if (mol%nat .eq. 1) return
!>  defaults
    tight = calc%optlev
    modef = 0
    call get_optthr(mol%nat,tight,calc,ethr,gthr)
    iupdat = calc%iupdat
    hlow = calc%hlow_opt !> 0.01 in ancopt, 0.002 too small
    hmax = calc%hmax_opt
    maxdispl = calc%maxdispl_opt
    gnorm = 0.0_wp
    depred = 0.0_wp
    echng = 0.0_wp
    alp = 1.0_wp
    alpold = 1.0_wp
    exact = calc%exact_rf .or. tight>0

    maxmicro = 100
    maxcycle = calc%maxcycle
    if (maxcycle .lt. maxmicro) maxmicro = maxcycle

    !> check if the molecule is linear
    call axis(mol%nat,mol%at,mol%xyz,rot,dumi)
    linear = (rot(3) .lt. 1.d-10).or.(mol%nat == 2)

    !> set degrees of freedom
    nat3 = 3*mol%nat
    nvar = nat3-6
    if (linear) then
      nvar = nat3-5
    end if
    if (calc%nfreeze .gt. 0) then ! exact fixing
      nvar = nat3-3*calc%nfreeze-3
      if (nvar .le. 0) nvar = 1
    end if

    !$omp critical
    allocate (pmode(nat3,1)) ! dummy allocated
    !$omp end critical

!>--- print a summary of settings, if desired
    if (pr) then
      call print_optsummary(calc,tight,nvar,maxcycle,maxmicro, &
      &                       ethr,gthr,linear,wr)
    end if

!>--- initialize OPT object
    !$omp critical
    allocate (h(nat3,nat3),hess(nat3*(nat3+1)/2),eig(nat3))
    call OPT%allocate2(mol%nat) !> NOTE: OPT%nvar will be nat*3 !!!
    allocate (molopt%at(mol%nat),molopt%xyz(3,mol%nat))
    nvar1 = OPT%nvar+1         
    npvar = OPT%nvar*(nvar1)/2 !> packed size of Hessian (note the abuse of nvar1!)
    allocate (gold(OPT%nvar),displ(OPT%nvar),grd1(OPT%nvar),source=0.0_wp)
    allocate(int_hess(size(OPT%hess)))
    allocate(c(nat3))
    allocate(IPIV(nat3))
    !$omp end critical

!>------------------------------------------------------------------------
!>--- put the Hessian guess into the type
!>------------------------------------------------------------------------
    !k = 0
    !do i = 1,nat3
    !  do j = 1,i
    !    k = k+1
    !    if (i /= j) then
    !      OPT%hess(k) = 0.0_wp
    !    else
    !      OPT%hess(k) = calc%hguess
    !    end if
    !  end do
    !end do

    call initialize_hessian(calc,calc%hess_init,mol%xyz,mol%nat,mol%at,OPT%hess,calc%hguess,pr)

!>--- backup coordinates, and starting energy
    molopt%nat = mol%nat
    molopt%at = mol%at
    molopt%xyz = mol%xyz
    molopt%wrextxyz = calc%logextxyz
    estart = etot

!>--- initialize .log file, if desired
    ilog = 942
    if (wr) then
      open (newunit=ilog,file=logfile)
    end if

    iter = 0

!>--- start with a printout of the preceeding single point
    if (pr) call print_optiter(iter)
    gnorm = norm2(grd)
    if (pr) then
      write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') etot
      write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') 0.0_wp
      write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
      write (*,'(2x,"predicted",e18.7)',advance='no') 0.0_wp
      write (*,'(1x,"("f7.2"%)")')-0.0_wp
    end if

!>======================================================================
    NR_iter: do while (iter < maxcycle.and..not.converged)
!>======================================================================
!>--- count the step and print out
      iter = iter+1
      if (pr) call print_optiter(iter)
      gold = reshape(grd, [nat3])
      gnold = gnorm
      eold = energy

!>--- calc predicted energy change based on E = E0 + delta * G + delta^2 * H
      if (iter > 1) then
        call prdechng(OPT%nvar,gold,displ,OPT%hess,depred)
      end if

!>------------------------------------------------------------------------
!>--- SINGLEPOINT CALCULATION
!>------------------------------------------------------------------------
      grd = 0.0_wp
      call engrad(molopt,calc,energy,grd,iostatus)
      if (iostatus .ne. 0) then
        fail = .true.
        exit NR_iter
      end if
      gnorm = norm2(grd)
      grd1 = reshape(grd, [nat3])

!>--- dump to .log file
      if (wr) then
        call molopt%appendlog(ilog,energy)
      end if

      if (gnorm .gt. 500.0_wp) then
        if (pr) write (*,*) '|grad| > 500, something is totally wrong!'
        fail = .true.
        iostatus = -1
        exit NR_iter
      end if

!>--- check for convergence
      gchng = gnorm-gnold
      echng = energy-eold
      econverged = abs(echng) .lt. ethr
      gconverged = gnorm .lt. gthr
      lowered = echng .lt. 0.0_wp

!>--- optimization step printout
      if (pr) then
        write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') energy
        write (*,'(5x,"change ΔE",e18.7,1x,"Eh")') echng
        write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
        write (*,'(2x,"predicted",e18.7)',advance='no') depred
        if (iter > 1) then
          dummy = (depred-echng)/echng*100.0_wp
          if (abs(dummy) < 1000.0_wp) then
            write (*,'(1x,"("f7.2"%)")') dummy
          else
            write (*,'(1x,"(*******%)")')
          end if
        else
          write (*,'(1x,"("f7.2"%)")')-100.0_wp
        end if
      end if

!>--- dynamic scaling in dependence of grad norm
!>--- if we are close to convergence we can take larger steps
      alpold = alp

      !alp = 1.0d-0
      !if (gnorm .lt. 0.002) then ! 0.002
      !  alp = 1.5d-0 ! 1.5
      !end if
      !if (gnorm .lt. 0.0006) then
      !  alp = 2.0d-0 ! 2
      !end if
      !if (gnorm .lt. 0.0003) then
      !  alp = 3.0d-1 ! 3
      !end if

      
      !alp = alp_generate(gnorm, calc%optlev,calc%opt_engine, calc%hess_init)
      alp = 1.0_wp
      !write(stdout,*) alp

!>------------------------------------------------------------------------
!> Update the Hessian
!>------------------------------------------------------------------------
      if (iter .gt. 1) then
!>--- Hessian update, but only after first iteration (iter > 1)
        select case (iupdat)
        case (0)
          call bfgs(OPT%nvar,gnorm,grd1,gold,displ*alpold,OPT%hess)
        case (1)
          call powell(OPT%nvar,gnorm,grd1,gold,displ*alpold,OPT%hess)
        case (2)
          call sr1(OPT%nvar,gnorm,grd1,gold,displ*alpold,OPT%hess)
        case (3)
          call bofill(OPT%nvar,gnorm,grd1,gold,displ*alpold,OPT%hess)
        case (4)
          call schlegel(OPT%nvar,gnorm,grd1,gold,displ*alpold,OPT%hess)
        case default
          write (*,*) 'invalid hessian update selection'
          stop
        end select
      end if

      !allocate(calc%chess%H(nat3,nat3))
      if (calc%full_HR) then
        call dhtosq(nat3,calc%chess%H(:,:),OPT%hess(:))
      end if
!>------------------------------------------------------------------------
!>  Newton Raphson (NR) method
!>------------------------------------------------------------------------

!> Solve linear system H*dx = -g
    int_hess = OPT%hess
    c = -grd1 !> This will become the update step after system is solved
    call dppsv('U', nat3, 1, int_hess, c, nat3, info2) 
    if (info2 /= 0) then !> dppsv assumes matrix to be positive definite, fallback dspsv
        int_hess = OPT%hess
        c = -grd1
        call dspsv('U', nat3, 1, int_hess, IPIV, c, nat3, info3)
        displ = c  
    else
        displ = c
    endif
    
!>--- rescale displacement if necessary
      maxd = alp*sqrt(ddot(OPT%nvar,displ,1,displ,1))
      if (maxd > maxdispl) then
        if (pr) write (*,'(" * rescaling step by",f14.7)') maxdispl/maxd
        displ = maxdispl*displ/maxd
      end if

!>--- now some output
      dsnrm = sqrt(ddot(OPT%nvar,displ,1,displ,1))
      if (pr) then
        !> this array is currently not used and will be overwritten in next step
        gold = abs(displ)
        imax(1) = maxloc(gold,1); gold(imax(1)) = 0.0_wp
        imax(2) = maxloc(gold,1); gold(imax(2)) = 0.0_wp
        imax(3) = maxloc(gold,1)
        write (*,'(3x,"displ. norm   :",f14.7,1x,"a.u.")',advance='no') &
          dsnrm*alp
        !write (*,'(6x,"lambda   ",e18.7)') eaug(1)
        write (*,'(3x,"maximum displ.:",f14.7,1x,"a.u.")',advance='no') &
          abs(displ(imax(1)))*alp
        write (*,'(6x,"in coords ",3("#",i0,", "),"...")') imax
      end if

!>------------------------------------------------------------------------
!>--- new coordinates
!>------------------------------------------------------------------------
      molopt%xyz = molopt%xyz+reshape(displ, [3,molopt%nat])*alp

!>--- converged ?
      econverged = abs(echng) .lt. ethr
      gconverged = gnorm .lt. gthr
      lowered = echng .lt. 0.0_wp
      converged = econverged.and.gconverged.and.lowered
      if (pr) then
        call print_convd(econverged,gconverged)
      end if
      if (converged) then
        converged = .true.
        etot = energy
        exit NR_iter
      end if

!>======================================================================
    end do NR_iter
!>======================================================================

!>--- close .log file
    if (wr) then
      close (ilog)
    end if

    if (converged) then
!>--- if the relaxation converged properly do this
      iostatus = 0
      if (pr) then
        call rcache%allocate(mol%nat)
        rmsdval = rmsd(mol,molopt,ccache=rcache)
        write (*,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
          "GEOMETRY OPTIMIZATION CONVERGED AFTER",iter,"ITERATIONS"
        write (*,'(72("-"))')
        write (*,'(1x,"total energy gain   :",F18.7,1x,"Eh",F14.4,1x,"kcal/mol")') &
          etot-estart, (etot-estart)*autokcal
        write (*,'(1x,"total RMSD          :",F18.7,1x,"a0",F14.4,1x,"Å")') &
          rmsdval,rmsdval*autoaa
        write (*,'(72("-"))')
      end if
    else if (iostatus .ne. 0) then
!>--- if iostatus =/= 0, something went wrong in the relaxation
      if (pr) then
        write (*,'(/,3x,"***",1x,a,1x,"***",/)') &
          "GEOMETRY RELAXATION FAILED"
      end if
    else
!>--- not converging in the given cycles is considered a FAILURE
      !> some iostatus>0 is selected to signal this
      iostatus = iter
      if (pr) then
        write (*,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
          "FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN",iter,"ITERATIONS"
      end if
    end if

!>--- overwrite input structure with optimized one
    mol%nat = molopt%nat
    mol%at = molopt%at
    mol%xyz = molopt%xyz

!> deallocate data
    !$omp critical
    if (allocated(gold)) deallocate (gold)
    if (allocated(displ)) deallocate (displ)
    if (allocated(grd1)) deallocate (grd1)
    if (allocated(pmode)) deallocate (pmode)
    if (allocated(h)) deallocate (h)
    if (allocated(hess)) deallocate (hess)
    if (allocated(molopt%at)) deallocate (molopt%at)
    if (allocated(molopt%xyz)) deallocate (molopt%xyz)
    call OPT%deallocate
    !$omp end critical

    return
  end subroutine newton_raphson


!========================================================================================!
!========================================================================================!
end module newton_raphson_module
