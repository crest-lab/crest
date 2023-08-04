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

module ancopt_module
  use iso_fortran_env,only:wp => real64,sp => real32
  use crest_calculator
  use axis_module
  use strucrd
  use ls_rmsd
  use testmol

  use type_anc
  use optimize_maths
  use modelhessian_module
  use hessupdate_module
  implicit none

  integer,parameter :: olev_crude = -3
  integer,parameter :: olev_sloppy = -2
  integer,parameter :: olev_loose = -1
  integer,parameter :: olev_lax = -4
  integer,parameter :: olev_normal = 0
  integer,parameter :: olev_tight = 1
  integer,parameter :: olev_vtight = 2
  integer,parameter :: olev_extreme = 3

  real(wp),parameter :: autoaa = 0.52917726_wp
  real(wp),parameter :: aatoau = 1.0_wp / autoaa
  real(wp),parameter :: autokcal = 627.509541_wp

  type :: convergence_log
    integer :: nlog
    real(wp),allocatable :: elog(:)
    real(wp),allocatable :: glog(:)
  contains
    procedure :: set_eg_log
    procedure :: get_averaged_energy
    procedure :: get_averaged_gradient
  end type convergence_log
  interface convergence_log
    module procedure new_convergence_log
  end interface convergence_log

  public :: ancopt
contains
!========================================================================================!
!> subroutine ancopt
!> Implementation of the Aproximate Normal Coordinate (ANC) optimizer
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
!>       wr  - logfile (crestopt.log) bool
!>  iostatus - return status of the routine
!>             (success=0, error<0, not converged>0)
!>-----------------------------------------------------------------------
  subroutine ancopt(mol,calc,etot,grd,pr,wr,iostatus)
    implicit none
    !> Inputs
    type(coord),intent(inout) :: mol
    type(calcdata),intent(in) :: calc
    real(wp),intent(inout) :: etot
    real(wp),intent(inout) :: grd(3,mol%nat)
    logical,intent(in) :: pr
    logical,intent(in) :: wr
    integer,intent(out) :: iostatus
    !> Local
    integer  :: tight
    real(wp) :: eel
    real(wp) :: et
    real(wp) :: egap
    logical :: fail
    !> Local objects
    type(coord)   :: molopt
    type(tb_anc)  :: anc
    type(mhparam) :: mhset

    real(wp) :: step,amu2au,au2cm,dumi,dumj,damp,hlow,edum,s6,thr
    real(wp) :: maxdispl,gthr,ethr,hmax,energy,rij(3),t1,t0,w1,w0
    real(wp) :: rot(3)
    integer :: n3,i,j,k,l,jjj,ic,jc,ia,ja,ii,jj,info,nat3
    integer :: nvar,iter,nread,maxcycle,maxmicro,itry,maxopt,iupdat,iii
    integer :: id,ihess,error
    integer :: ilog
    real(wp),allocatable :: h(:,:)
    real(wp),allocatable :: b(:,:)
    real(wp),allocatable :: fc(:)
    real(wp),allocatable :: eig(:)
    real(wp),allocatable :: aux(:)
    real(wp),allocatable :: hess(:)
    integer,allocatable :: iwork(:)
    integer,allocatable :: totsym(:)
    real(wp),allocatable :: pmode(:,:)
    real(wp),allocatable :: grmsd(:,:)
    type(convergence_log),allocatable :: avconv
    real(wp) :: U(3,3),x_center(3),y_center(3),rmsdval
    integer :: modef
    logical :: ex,converged,linear
    real(wp) :: estart,esave
    character(len=*),parameter :: scifmt = &
                                  '(10x,"│",3x,a,e22.7,1x,a,1x,"│")'
    character(len=*),parameter :: dblfmt = &
                                  '(10x,"│",3x,a,f18.7,5x,a,1x,"│")'
    character(len=*),parameter :: intfmt = &
                                  '(10x,"│",3x,a,i18,      10x,"│")'
    character(len=*),parameter :: chrfmt = &
                                  '(10x,"│",3x,a,a18,      10x,"│")'

    iostatus = 0
    fail = .false.
    converged = .false.
    if (mol%nat .eq. 1) return
!>  defaults
    tight = calc%optlev
    modef = 0
    call get_optthr(mol%nat,tight,calc,ethr,gthr)
    iupdat = calc%iupdat
    hlow = calc%hlow_opt !optset%hlow_opt !> 0.01 in ancopt, 0.002 too small
    hmax = calc%hmax_opt
    maxdispl = calc%maxdispl_opt !optset%maxdispl_opt
    s6 = mhset%s6 !> slightly better than 30 for various proteins

!> initial number of steps in relax() routine before
!> new ANC are made by model Hessian
!> increased during opt.
    maxmicro = calc%micro_opt
    maxcycle = calc%maxcycle
    if (maxcycle .lt. maxmicro) maxmicro = maxcycle

    !> check if the molecule is linear
    call axis(mol%nat,mol%at,mol%xyz,rot,dumi)
    linear = rot(3) .lt. 1.d-10

    !> set degrees of freedom
    nat3 = 3 * mol%nat
    nvar = nat3 - 6
    if (linear) then
      nvar = nat3 - 5
    end if
    !if (fixset%n .gt. 0) then ! exact fixing
    !  nvar = nat3 - 3 * fixset%n - 3
    !  if (nvar .le. 0) nvar = 1
    !end if
    
    !$omp critical
    allocate (pmode(nat3,1),grmsd(3,mol%nat)) ! dummy allocated
    !$omp end critical

!>--- print a summary of settings, if desired
    if (pr) then
      write (*,'(/,10x,"┍",49("━"),"┑")')
      write (*,'(10x,"│",11x,a,11x,"│")') "GEOMETRY OPTIMIZATION SETUP"
      write (*,'(10x,"┝",49("━"),"┥")')
      !write (*,chrfmt) "optimization level",int2optlevel(tight)
      write (*,intfmt) "optimization level",tight
      write (*,intfmt) "max. optcycles    ",maxcycle
      write (*,intfmt) "ANC micro-cycles  ",maxmicro
      write (*,intfmt) "degrees of freedom",nvar
      if (modef > 0) then
        write (*,intfmt) "# mode follow     ",modef
      end if
      write (*,'(10x,"├",49("─"),"┤")')
      if (calc%exact_rf) then
        write (*,chrfmt) "RF solver         ","spevx"
      else
        write (*,chrfmt) "RF solver         ","davidson"
      end if
      select case( iupdat )
      case( 0 )
       write (*,chrfmt) "Hessian update    ","bfgs"
      case( 1 )
       write (*,chrfmt) "Hessian update    ","powell"
      case( 2 )
       write (*,chrfmt) "Hessian update    ","sr1"
      case( 3 )
       write (*,chrfmt) "Hessian update    ","bofill"
      case( 4 )
       write (*,chrfmt) "Hessian update    ","schlegel"
      end select
      write (*,chrfmt) "write crestopt.log",bool2string(wr)
      if (linear) then
        write (*,chrfmt) "linear (good luck)",bool2string(linear)
      else
        write (*,chrfmt) "linear?           ",bool2string(linear)
      end if
      write (*,scifmt) "energy convergence",ethr,"Eh  "
      write (*,scifmt) "grad. convergence ",gthr,"Eh/α"
      write (*,dblfmt) "maximium RF displ.",maxdispl,"    "
      write (*,scifmt) "Hlow (freq-cutoff)",hlow,"    "
      write (*,dblfmt) "Hmax (freq-cutoff)",hmax,"    "
      write (*,dblfmt) "S6 in model hess. ",s6,"    "
      write (*,'(10x,"└",49("─"),"┘")')
    end if

!>--- initialize anc object
    !$omp critical
    allocate (h(nat3,nat3),hess(nat3 * (nat3 + 1) / 2),eig(nat3))
    call anc%allocate(mol%nat,nvar,hlow,hmax)
    allocate(molopt%at(mol%nat),molopt%xyz(3,mol%nat))
    !$omp end critical

!>--- backup coordinates, and starting energy
    molopt%nat = mol%nat
    molopt%at  = mol%at
    molopt%xyz = mol%xyz
    estart = etot

!>--- initialize .log file, if desired
    ilog = 942
    if (wr) then
      open (newunit=ilog,file='crestopt.log')
    end if

!>--- The ANCOPT iteration loop. "iter" is updated in relax() subroutine
    iter = 0
!>======================================================================
    ANC_microiter: do while (iter < maxcycle .and. .not. converged)
!>======================================================================
!>--- generate model Hessian
      if (pr) write (*,'(/,''generating ANC from model Hessian ...'')')
      call modhes(calc,mhset,molopt%nat,molopt%xyz,molopt%at,hess,pr)

!>--- project trans. and rot. from Hessian
      if (.not. linear) then
        call trproj(molopt%nat,nat3,molopt%xyz,hess,.false.,0,pmode,1) ! normal
      end if

!>--- ANC generation (requires blowup)
      k = 0
      do i = 1,nat3
        do j = 1,i
          k = k + 1
          h(i,j) = hess(k)
          h(j,i) = hess(k)
        end do
      end do
      call anc%new(molopt%xyz,h,pr,linear,fail)
      if (fail) then
        iostatus = -1
        exit ANC_microiter
      end if

      esave = etot !> save energy before relaxation
!>--- call the actual relaxation routine
!>    this routine will perform [maxmicro] relaxation steps
      call relax(molopt,calc,anc,iter,maxmicro,etot,grd,  &
            &      ethr,gthr,converged,                  &
            &      pr,wr,ilog,iostatus,avconv)
      if (iostatus .ne. 0) then
        if (pr) write (*,*) 'Structure relaxation failed'
        exit ANC_microiter
      end if

!>--- update max. iterations for next relax() call
      maxmicro = min(int(float(maxmicro) * 1.1_wp),2 * calc%micro_opt)

!>--- check structural change by RMSD
      call rmsd(molopt%nat,anc%xyz,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)
      if (.not. converged .and. pr) then
        write (*,'(" * RMSD in coord.:",f14.7,1x,"α")',advance='no') rmsdval
        write (*,'(6x,"energy gain",e16.7,1x,"Eh")') etot - esave
      end if
!>======================================================================
    end do ANC_microiter
!>======================================================================

!>--- close .log file
    if (wr) then
      close (ilog)
    end if

    if (converged) then
!>--- if the relaxation converged properly do this
      iostatus = 0
      if (pr) then
        call rmsd(mol%nat,mol%xyz,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)
        write (*,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
          "GEOMETRY OPTIMIZATION CONVERGED AFTER",iter,"ITERATIONS"
        write (*,'(72("-"))')
        write (*,'(1x,"total energy gain   :",F18.7,1x,"Eh",F14.4,1x,"kcal/mol")') &
          etot - estart, (etot - estart) * autokcal
        write (*,'(1x,"total RMSD          :",F18.7,1x,"a0",F14.4,1x,"Å")') &
          rmsdval,rmsdval * autoaa
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
    mol%at  = molopt%at
    mol%xyz = molopt%xyz


!> deallocate data
    !$omp critical
    if (allocated(grmsd)) deallocate (grmsd)
    if (allocated(pmode)) deallocate (pmode)
    if (allocated(h)) deallocate (h)
    if (allocated(hess)) deallocate (hess)
    if(allocated(molopt%at)) deallocate( molopt%at )
    if(allocated(molopt%xyz)) deallocate(molopt%xyz)
    call anc%deallocate
    !$omp end critical

    return
  end subroutine ancopt

!========================================================================================!
!> subroutine relax
!>
!> Implements the microiteration relaxation cycles, i.e.,
!> the update steps and diagonalizations between the
!> new ANC generation.
!>--------------------------------------------------------
  subroutine relax(mol,calc,anc,iter,maxmicro,etot,grd,  &
            &      ethr,gthr,converged,                  &
            &      pr,wr,ilog,iostatus,avconv)
    implicit none

    type(coord) :: mol
    type(calcdata) :: calc
    type(tb_anc) :: anc
    integer,intent(inout) :: iter
    integer,intent(in)    :: maxmicro
    real(wp),intent(inout) :: etot
    real(wp),intent(inout) :: grd(3,mol%nat)
    real(wp),intent(in) :: ethr
    real(wp),intent(in) :: gthr
    logical,intent(in) :: pr
    logical,intent(in) :: wr
    integer,intent(in) :: ilog
    logical,intent(out) :: converged
    integer,intent(out) :: iostatus
    type(convergence_log),intent(inout),optional :: avconv

    integer :: maxiter
    integer :: iupdat

    real(wp) :: et
    real(wp) :: egap
    real(wp) :: acc_in

    logical :: exact
    integer :: nvar1,npvar,npvar1
    logical :: fail
    logical :: econverged
    logical :: gconverged
    logical :: lowered
    integer :: i,j,ii,jj,k,lwork,info,m,idum,imax(3)
    real(wp) :: energy,dsnrm,maxdispl,t0,w0,t1,w1
    real(wp) :: lambda,gnorm,dnorm,ddot,eold,xdum,estart,acc,e_in
    real(wp) :: depred,echng,dummy,maxd,alp,gchng,gnold
    real(wp),allocatable :: gold(:)
    real(wp),allocatable :: displ(:),gint(:)
    real(sp),allocatable :: eaug(:)
    real(sp),allocatable :: Uaug(:,:)
    real(sp),allocatable :: Aaug(:)
    real(sp),parameter :: r4dum = 1.e-8
    !> LAPACK & BLAS
    external :: dgemv
    real(sp),external :: sdot

    iostatus = 0

    !$omp critical
    allocate (gold(anc%nvar),displ(anc%nvar),gint(anc%nvar),source=0.0_wp)

    gnorm = 0.0_wp
    depred = 0.0_wp
    echng = 0.0_wp
    maxdispl = calc%maxdispl_opt
    acc_in = calc%acc_opt
    energy = etot
    e_in = etot
    alp = 1.0_wp
    converged = .false.
    exact = calc%exact_rf
    iupdat = calc%iupdat

    nvar1 = anc%nvar + 1             !> dimension of RF calculation
    npvar = anc%nvar * (nvar1) / 2   !> packed size of Hessian (note the abuse of nvar1!)
    npvar1 = nvar1 * (nvar1 + 1) / 2 !> packed size of augmented Hessian
    allocate (Uaug(nvar1,1),eaug(nvar1),Aaug(npvar1))
    !$omp end critical

!! ========================================================================
    main_loop: do ii = 1,maxmicro
!! ========================================================================
      iter = iter + 1
      if (pr) &
        write (*,'(/,"┌",76("─"),"┐",/,"│",32(" ")," CYCLE",i5,1x,32(" "),"│",/,"└",76("─"),"┘")') iter

      gold = gint
      gnold = gnorm
      eold = energy
!>--- calc predicted energy change based on E = E0 + delta * G + delta^2 * H
      if (ii > 1) then
        call prdechng(anc%nvar,gold,displ,anc%hess,depred)
      end if

!>------------------------------------------------------------------------
!>--- SINGLEPOINT CALCULATION
!>------------------------------------------------------------------------
      !> get Cartestian coordinates and gradient
      call anc%get_cartesian(mol%xyz)

      grd = 0.0_wp
      call engrad(mol,calc,energy,grd,iostatus)
      if (iostatus .ne. 0) then
        fail = .true.
        exit main_loop
      end if
!>------------------------------------------------------------------------

!>--- dump to .log file
      if (wr) then
        call mol%appendlog(ilog,energy)
      end if
!>--- transform Cartesian xyz to internal gradient
      call dgemv('t',anc%n3,anc%nvar,1.0_wp,anc%B,anc%n3,grd,1,0.0_wp,gint,1)
      gnorm = norm2(gint)

      if (gnorm .gt. 500.0_wp) then
        if (pr) write (*,*) '|grad| > 500, something is totally wrong!'
        fail = .true.
        iostatus = -1
        exit main_loop
      end if

      !if (present(avconv)) then
      !  call avconv%set_eg_log(energy,gnorm)
      !  energy = avconv%get_averaged_energy()
      !  gnorm = avconv%get_averaged_gradient()
      !  if (pr) then
      !    write (*,'("av. E:",1x,f14.7,1x,"->",1x,f14.7)') &
      !      avconv%elog(avconv%nlog),energy
      !    write (*,'("av. G:",1x,f14.7,1x,"->",1x,f14.7)') &
      !      avconv%glog(avconv%nlog),gnorm
      !  end if
      !end if

!>--- check for convergence
      gchng = gnorm - gnold
      echng = energy - eold
      econverged = abs(echng) .lt. ethr
      gconverged = gnorm .lt. gthr
      lowered = echng .lt. 0.0_wp

!>--- optimization step printout
      if (pr) then
        write (*,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') energy
        write (*,'(5x,"change   ",e18.7,1x,"Eh")') echng
        write (*,'(3x,"gradient norm :",f14.7,1x,"Eh/α")',advance='no') gnorm
        write (*,'(3x,"predicted",e18.7)',advance='no') depred
        if (ii > 1) then
          dummy = (depred - echng) / echng * 100.0_wp
          if (abs(dummy) < 1000.0_wp) then
            write (*,'(1x,"("f7.2"%)")') dummy
          else
            write (*,'(1x,"(*******%)")')
          end if
        else
          write (*,'(1x,"("f7.2"%)")') - 100.0_wp
        end if
      end if

      if (gnorm .lt. 0.002) then ! 0.002
        alp = 1.5d0 ! 1.5
      elseif (gnorm .lt. 0.0006) then
        alp = 2.0d0 ! 2
      elseif (gnorm .lt. 0.0003) then
        alp = 3.0d0 ! 3
      else
        alp = 1.0d0
      end if

!>------------------------------------------------------------------------
!> Update the Hessian
!>------------------------------------------------------------------------
      if (ii .gt. 1) then
!>--- Hessian update, but only after first iteration (ii > 1)
        select case( iupdat )
        case( 0 )
           call bfgs(anc%nvar,gnorm,gint,gold,displ,anc%hess)
        case( 1 )
           call powell(anc%nvar,gnorm,gint,gold,displ,anc%hess)
        case( 2 )
           call sr1(anc%nvar,gnorm,gint,gold,displ,anc%hess)
        case( 3 )
           call bofill(anc%nvar,gnorm,gint,gold,displ,anc%hess)
        case( 4 )
           call schlegel(anc%nvar,gnorm,gint,gold,displ,anc%hess)
        case default
           write(*,*) 'invalid hessian update selection'
           stop
        end select
      end if

!>------------------------------------------------------------------------
!>  rational function (RF) method
!>------------------------------------------------------------------------
!>  To get initial guess for the displacement solve this:
!>  (Note: SINGLE PRECISION accuracy, there are some typecasts!)
!>
!>   ⎛ H  g ⎞ ⎛ dx ⎞     ⎛ dx ⎞
!>   ⎝ g  0 ⎠ ⎝  1 ⎠ = λ ⎝  1 ⎠
!>     Aaug    Uaug       Uaug

!>--- first, augment Hessian by gradient, everything packed, no blowup
      Aaug(1:npvar) = real(anc%hess(1:npvar),sp)
      Aaug(npvar + 1:npvar1 - 1) = real(gint(1:anc%nvar),sp)
      Aaug(npvar1) = 0.0_sp

!>--- choose solver
      if (exact .or. nvar1 .lt. 50) then
        call solver_sspevx(nvar1,r4dum,Aaug,Uaug,eaug,fail)
      else
        !>--- steepest decent guess for displacement
        if (ii .eq. 1) then
          Uaug(:,1) = [-real(gint(1:anc%nvar),sp),1.0_sp]
          dsnrm = sqrt(sdot(nvar1,Uaug,1,Uaug,1))
          Uaug = Uaug / real(dsnrm,sp)
        end if
        call solver_sdavidson(nvar1,r4dum,Aaug,Uaug,eaug,fail,.false.)
        !>--- if that failed, retry with better solver
        if (fail) then
          call solver_sspevx(nvar1,r4dum,Aaug,Uaug,eaug,fail)
        end if
      end if

!>--- divide by last element(=λ) to get the displacement vector dx
      if (fail .or. abs(Uaug(nvar1,1)) .lt. 1.e-10) then
        if (pr) write (*,*) "internal rational function error"
        iostatus = -1
        exit main_loop
      end if
      displ(1:anc%nvar) = Uaug(1:anc%nvar,1) / Uaug(nvar1,1)

!>--- check if step is too large, just cut off everything thats to large
      !do j = 1,anc%nvar
      !  if (abs(displ(j)) .gt. maxdispl) then
      !    if (displ(j) < 0) displ(j) = -maxdispl
      !    if (displ(j) > 0) displ(j) = maxdispl
      !  end if
      !end do
!>--- maybe more consistent version is to rescale displacement
      maxd = alp * sqrt(ddot(anc%nvar,displ,1,displ,1))
      if (maxd > maxdispl) then
        if (pr) write (*,'(" * rescaling step by",f14.7)') maxdispl / maxd
        displ = maxdispl * displ / maxd
      end if

!>--- now some output
      dsnrm = sqrt(ddot(anc%nvar,displ,1,displ,1))
      if (pr) then
        !> this array is currently not used and will be overwritten in next step
        gold = abs(displ)
        imax(1) = maxloc(gold,1); gold(imax(1)) = 0.0_wp
        imax(2) = maxloc(gold,1); gold(imax(2)) = 0.0_wp
        imax(3) = maxloc(gold,1)
        write (*,'(3x,"displ. norm   :",f14.7,1x,"α")',advance='no') &
          dsnrm * alp
        write (*,'(6x,"lambda   ",e18.7)') eaug(1)
        write (*,'(3x,"maximum displ.:",f14.7,1x,"α")',advance='no') &
          abs(displ(imax(1))) * alp
        write (*,'(6x,"in ANC''s ",3("#",i0,", "),"...")') imax
        !call prdispl(anc%nvar,displ)
      end if
!>------------------------------------------------------------------------

!>--- 2nd: exit and redo hessian (internal restart)
      if (ii .gt. 2 .and. dsnrm .gt. 2.0) then
        if (pr) write (*,*) 'exit because of too large step'
        exit main_loop
      end if

!>--- new coordinates
      anc%coord = anc%coord + displ * alp

!>--- converged ?
      econverged = abs(echng) .lt. ethr
      gconverged = gnorm .lt. gthr
      lowered = echng .lt. 0.0_wp
      converged = econverged .and. gconverged .and. lowered
      if (pr) then
        write (*,'(3x,"converged δE/grad :",1x,l," /",l)') econverged,gconverged
      end if
      if (converged) then
        !if (econverged .and. gconverged .and. lowered) then
        converged = .true.
        etot = energy
        !return
        exit main_loop
      end if
!>========================================================================
    end do main_loop
!>========================================================================
    !$omp critical
    if (allocated(Uaug)) deallocate (Uaug)
    if (allocated(eaug)) deallocate (eaug)
    if (allocated(Aaug)) deallocate (Aaug)
    !$omp end critical

    etot = energy
    call anc%get_cartesian(mol%xyz)

    return
  end subroutine relax

!========================================================================================!
!> subroutine get_optthr
!> routine to set some optimization thresholds
  subroutine get_optthr(n,olev,calc,ethr,gthr)
    implicit none
    integer,intent(in) :: n
    integer,intent(in) :: olev
    type(calcdata) :: calc
    real(wp),intent(out) :: ethr
    real(wp),intent(out) :: gthr
    integer :: maxcycle
    real(wp) :: acc
    select case (olev)
!> very approximate = crude
    case (olev_crude)
      ethr = 5.d-4
      gthr = 1.d-2
      maxcycle = n
      acc = 3.00d0
!> approximate = sloopy
    case (olev_sloppy)
      ethr = 1.d-4
      gthr = 6.d-3
      maxcycle = n
      acc = 3.00d0
!> loose
    case (olev_loose)
      ethr = 5.d-5
      gthr = 4.d-3
      maxcycle = n * 2
      acc = 2.00d0
!>  for DCOSMO-RS opts with TM i.e. between loose and normal, keyword "lax"
    case (olev_lax)
      ethr = 2.d-5
      gthr = 2.5d-3
      maxcycle = n * 2
      acc = 2.00d0
!> normal
    case default
      ethr = 5.d-6
      gthr = 1.d-3
      maxcycle = n * 3
      acc = 1.0d0
!> tight
    case (olev_tight)
      ethr = 1.d-6
      gthr = 8.d-4
      maxcycle = n * 5
      acc = 0.20d0
!> very tight
    case (olev_vtight)
      ethr = 1.d-7
      gthr = 2.d-4
      maxcycle = n * 20
      acc = 0.05d0
!> extreme
    case (olev_extreme)
      ethr = 5.d-8
      gthr = 5.d-5
      maxcycle = n * 20
      acc = 0.01d0
    end select
    maxcycle = min(maxcycle,10000)
    maxcycle = max(maxcycle,200)
    if (calc%maxcycle <= 0) then
      calc%maxcycle = maxcycle
    end if
    !calc%iupdat = 0 !0=BFGS, 1=Powell
    if (calc%tsopt) then
      calc%hlow_opt = max(calc%hlow_opt,0.250d0)
      calc%iupdat = 1
    end if
    return
  end subroutine get_optthr

!========================================================================================!
!> Purpose:
!> Calculates predicted energy change according to the second order
!> model.
!>
!> Input:
!> nat3  - 3*natoms
!> hess  - Hessian matrix stored as lower triangle
!> grad  - Gradient vector
!> displ - Displacement vector
!>
!> Output:
!> depred - Predicted energy change
!---------------------------------------------------------------------
  subroutine prdechng(nat3,grad,displ,hess,depred)
    implicit none
    !> Input:
    integer,intent(in) :: nat3
    real(wp),intent(in) :: grad(nat3)
    real(wp),intent(in) :: displ(nat3)
    real(wp),intent(in) :: hess(nat3 * (nat3 + 1) / 2)
    !> Output:
    real(wp),intent(out) :: depred
    !> Local:
    real(wp),allocatable :: hdx(:)
    real(wp) :: gtmp,htmp
    !> BLAS functions:
    real(wp),external :: ddot
    external :: dspmv
    allocate (hdx(nat3),source=0.0_wp)
    call dspmv('u',nat3,0.5d0,hess,displ,1,0.0d0,hdx,1)
    gtmp = ddot(nat3,displ,1,grad,1)
    htmp = ddot(nat3,displ,1,hdx,1)
    depred = htmp + gtmp
    return
  end subroutine prdechng

  subroutine trfp2xyz(nvar,nat3,p,xyz0,h,dspl)
    implicit none
    integer,intent(in) :: nat3
    integer,intent(in) :: nvar
    integer :: nat,icount,i,j,k
    real(wp),intent(in) :: xyz0(3,nat3 / 3)
    real(wp),intent(out) :: dspl(3,nat3 / 3)
    real(wp),intent(in) :: h(nat3,nat3)
    real(wp),intent(in) :: p(nvar)
    real(wp) :: dum

    dspl = 0.0d0
    nat = nat3 / 3

! generate cartesian displacement vector
    do i = 1,nvar
      icount = 0
      do j = 1,nat
        do k = 1,3
          icount = icount + 1
          dum = h(icount,i) * p(i)
          dspl(k,j) = dspl(k,j) + dum
        end do
      end do
    end do

    dspl = dspl + xyz0

    return
  end subroutine trfp2xyz

  subroutine prdispl(nvar,displ)
    implicit none
    integer,intent(in) :: nvar
    real(wp),intent(in) :: displ(nvar)
    real(wp),allocatable :: er(:)
    integer,allocatable :: merk(:)
    integer :: i,j,ii,k
    integer :: ihilf
    real(wp) :: pp
    allocate (er(nvar),source=0.0_wp)
    allocate (merk(nvar),source=0)

    er = abs(displ)

    do i = 1,nvar
      merk(i) = i
    end do
    do ii = 2,nvar
      i = ii - 1
      k = i
      pp = er(i)
      do j = ii,nvar
        if (er(j) .le. pp) cycle
        k = j
        pp = er(j)
      end do
      if (k .eq. i) cycle
      er(k) = er(i)
      er(i) = pp
      ihilf = merk(i)
      merk(i) = merk(k)
      merk(k) = ihilf
    end do

    write (*,'(''Largest |displ|/coords:'',5(f8.4,'' ('',i4,'')''))') &
      (er(i),merk(i),i=1,min(3,nvar))

  end subroutine prdispl

!========================================================================================!
  function bool2string(bool)
    implicit none
    character(len=:),allocatable :: bool2string
    logical :: bool
    if (bool) then
      bool2string = "True"
    else
      bool2string = "False"
    end if
  end function bool2string

!========================================================================================!
  subroutine geoconvav(nc,e,g,val,deriv)
    implicit none
    integer :: nc     !> total number of E/G points
    real(wp) :: e(*)  !> total energy in Eh
    real(wp) :: g(*)  !> norm of Cartesian gradient (in TM: |dE/dxyz|)
    real(wp) :: val   !> av. energy in Eh to be used further
    real(wp) :: deriv !> av. gradient

    integer :: low
    integer :: i,j
    integer,parameter:: nav = 5 !> average over last nav
    real(wp) :: eav,gav

    !> only apply it if sufficient number of points i.e. a "tail" can exist
    !> with the censo blockl = 8 default, this can first be effective in the second
    if (nc .lt. 3 * nav) then
      val = e(nc)
      deriv = g(nc)
      return
    end if

    low = max(1,nc - nav + 1)
    j = 0
    eav = 0
    do i = nc,low,-1
      j = j + 1
      eav = eav + e(i)
      gav = gav + g(i)
    end do
    val = eav / float(j)

    low = max(1,nc - nav + 1)
    j = 0
    gav = 0
    do i = nc,low,-1
      j = j + 1
      gav = gav + g(i)
    end do
    ! adjust the gradient norm to xtb "conventions" because e.g. a noisy
    ! DCOSMO-RS gradient for large cases can never (even on average)
    ! become lower than the "-opt normal" thresholds
    deriv = gav / float(j) / 2.d0
  end subroutine geoconvav

  pure function new_convergence_log(nmax) result(self)
    integer,intent(in) :: nmax
    type(convergence_log) :: self
    self%nlog = 0
    allocate (self%elog(nmax))
    allocate (self%glog(nmax))
  end function new_convergence_log

  pure function get_averaged_energy(self) result(val)
    class(convergence_log),intent(in) :: self
    real(wp) :: eav,val
    integer :: i,j,low
    integer,parameter :: nav = 5

    ! only apply it if sufficient number of points i.e. a "tail" can exist
    ! with the censo blockl = 8 default, this can first be effective in the second
    if (self%nlog .lt. 3 * nav) then
      val = self%elog(self%nlog)
    else
      low = max(1,self%nlog - nav + 1)
      j = 0
      eav = 0
      do i = self%nlog,low,-1
        j = j + 1
        eav = eav + self%elog(i)
      end do
      val = eav / float(j)
    end if

  end function get_averaged_energy

  pure function get_averaged_gradient(self) result(deriv)
    class(convergence_log),intent(in) :: self
    real(wp) :: gav,deriv
    integer :: i,j,low
    integer,parameter :: nav = 5

    ! only apply it if sufficient number of points i.e. a "tail" can exist
    ! with the censo blockl = 8 default, this can first be effective in the second
    if (self%nlog .lt. 3 * nav) then
      deriv = self%glog(self%nlog)
    else
      low = max(1,self%nlog - nav + 1)
      j = 0
      gav = 0
      do i = self%nlog,low,-1
        j = j + 1
        gav = gav + self%glog(i)
      end do
      ! adjust the gradient norm to xtb "conventions" because e.g. a noisy
      ! DCOSMO-RS gradient for large cases can never (even on average)
      ! become lower than the "-opt normal" thresholds
      deriv = gav / float(j) / 2.d0
    end if

  end function get_averaged_gradient

  pure subroutine set_eg_log(self,e,g)
    class(convergence_log),intent(inout) :: self
    real(wp),intent(in) :: e,g
    real(wp),allocatable :: dum(:)
    integer :: k,k2
    k = size(self%elog)
    if (self%nlog >= k) then
      k2 = k + 1
      allocate (dum(k2))
      dum(1:k) = self%elog(1:k)
      call move_alloc(dum,self%elog)
      allocate (dum(k2))
      dum(1:k) = self%glog(1:k)
      call move_alloc(dum,self%glog)
    end if
    self%nlog = self%nlog + 1
    self%elog(self%nlog) = e
    self%glog(self%nlog) = g
  end subroutine set_eg_log

!========================================================================================!
  subroutine rdhess(nat3,h,fname)
    integer,intent(in)  :: nat3
    real(wp),intent(out) :: h(nat3,nat3)
    character(len=*),intent(in) :: fname
    integer  :: iunit,i,j,mincol,maxcol
    character(len=5)  :: adum
    character(len=80) :: a80

    !     write(*,*) 'Reading Hessian <',trim(fname),'>'
    open (newunit=iunit,file=fname)
50  read (iunit,'(a)') a80
    if (index(a80,'$hessian') .ne. 0) then
      do i = 1,nat3
        maxcol = 0
200     mincol = maxcol + 1
        maxcol = min(maxcol + 5,nat3)
        read (iunit,*) (h(j,i),j=mincol,maxcol)
        if (maxcol .lt. nat3) goto 200
      end do
      close (iunit)
      goto 300
    end if
    goto 50

300 return
  end subroutine rdhess

!========================================================================================!
end module ancopt_module
