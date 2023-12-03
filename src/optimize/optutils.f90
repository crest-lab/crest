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

!> this module implements several utility routines for geometry optimization
!> including the threshold selection

module optimize_utils
  use crest_parameters
  use crest_calculator
  use axis_module
  use strucrd
  use ls_rmsd
  use testmol

  use optimize_type
  use optimize_maths
  use modelhessian_module
  use hessupdate_module
  implicit none
  public !> all public

  integer,parameter :: olev_crude = -3
  integer,parameter :: olev_sloppy = -2
  integer,parameter :: olev_loose = -1
  integer,parameter :: olev_lax = -4
  integer,parameter :: olev_normal = 0
  integer,parameter :: olev_tight = 1
  integer,parameter :: olev_vtight = 2
  integer,parameter :: olev_extreme = 3

  public :: get_optthr

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine get_optthr(n,olev,calc,ethr,gthr)
!************************************************
!* subroutine get_optthr
!* routine to set some optimization thresholds
!************************************************
    implicit none
    integer,intent(in) :: n
    integer,intent(inout) :: olev
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
      maxcycle = n*2
      acc = 2.00d0
!>  for DCOSMO-RS opts with TM i.e. between loose and normal, keyword "lax"
    case (olev_lax)
      ethr = 2.d-5
      gthr = 2.5d-3
      maxcycle = n*2
      acc = 2.00d0
!> normal
    case default
      ethr = 5.d-6
      gthr = 1.d-3
      maxcycle = n*3
      acc = 1.0d0
!> tight
    case (olev_tight)
      ethr = 1.d-6
      gthr = 8.d-4
      maxcycle = n*5
      acc = 0.20d0
!> very tight
    case (olev_vtight)
      ethr = 1.d-7
      gthr = 2.d-4
      maxcycle = n*20
      acc = 0.05d0
!> extreme
    case (olev_extreme)
      ethr = 5.d-8
      gthr = 5.d-5
      maxcycle = n*20
      acc = 0.01d0
    end select
    maxcycle = min(maxcycle,10000)
    maxcycle = max(maxcycle,200)
!> user-defined overwrites
    if (calc%maxcycle <= 0) then
      calc%maxcycle = maxcycle
    end if
    if (calc%ethr_opt > 0.0_wp) then
      ethr = calc%ethr_opt
      olev = 99
    end if
    if (calc%gthr_opt > 0.0_wp) then
      gthr = calc%gthr_opt
      olev = 99
    end if
    if (calc%tsopt) then
      calc%hlow_opt = max(calc%hlow_opt,0.250d0)
      calc%iupdat = 1
    end if

    return
  end subroutine get_optthr

!========================================================================================!

  subroutine prdechng(nat3,grad,displ,hess,depred)
!********************************************************************
!* Purpose:
!* Calculates predicted energy change according to the second order
!* model.
!*
!* Input:
!* nat3  - 3*natoms
!* hess  - Hessian matrix stored as lower triangle
!* grad  - Gradient vector
!* displ - Displacement vector
!*
!* Output:
!* depred - Predicted energy change
!******************************************************************
    implicit none
    !> INPUT
    integer,intent(in) :: nat3
    real(wp),intent(in) :: grad(nat3)
    real(wp),intent(in) :: displ(nat3)
    real(wp),intent(in) :: hess(nat3*(nat3+1)/2)
    !> OUTPUT
    real(wp),intent(out) :: depred
    !> LOCAL
    real(wp),allocatable :: hdx(:)
    real(wp) :: gtmp,htmp
    !> BLAS functions:
    real(wp),external :: ddot
    external :: dspmv
    allocate (hdx(nat3),source=0.0_wp)
    call dspmv('u',nat3,0.5d0,hess,displ,1,0.0d0,hdx,1)
    gtmp = ddot(nat3,displ,1,grad,1)
    htmp = ddot(nat3,displ,1,hdx,1)
    depred = htmp+gtmp
    return
  end subroutine prdechng

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
      i = ii-1
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
    if (nc .lt. 3*nav) then
      val = e(nc)
      deriv = g(nc)
      return
    end if

    low = max(1,nc-nav+1)
    j = 0
    eav = 0
    do i = nc,low,-1
      j = j+1
      eav = eav+e(i)
      gav = gav+g(i)
    end do
    val = eav/float(j)

    low = max(1,nc-nav+1)
    j = 0
    gav = 0
    do i = nc,low,-1
      j = j+1
      gav = gav+g(i)
    end do
    ! adjust the gradient norm to xtb "conventions" because e.g. a noisy
    ! DCOSMO-RS gradient for large cases can never (even on average)
    ! become lower than the "-opt normal" thresholds
    deriv = gav/float(j)/2.d0
  end subroutine geoconvav

!========================================================================================!
  subroutine rdhess(nat3,h,fname)
    integer,intent(in)  :: nat3
    real(wp),intent(out) :: h(nat3,nat3)
    character(len=*),intent(in) :: fname
    integer  :: iunit,i,j,mincol,maxcol
    character(len=5)  :: adum
    character(len=80) :: a80

    open (newunit=iunit,file=fname)
50  read (iunit,'(a)') a80
    if (index(a80,'$hessian') .ne. 0) then
      do i = 1,nat3
        maxcol = 0
200     mincol = maxcol+1
        maxcol = min(maxcol+5,nat3)
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

  subroutine print_optiter(iter)
    implicit none
    integer :: iter
    write (*,'(/,"┌",76("─"),"┐")')
    write (*,'(  "│",32(" ")," CYCLE",i5,1x,32(" "),"│")') iter
    write (*,'(  "└",76("─"),"┘")')
  end subroutine print_optiter

!========================================================================================!

  subroutine print_optsummary(calc,tight,nvar,maxcycle,maxmicro, &
      &                       ethr,gthr,linear,wr)
    implicit none
    !> INPUT
    type(calcdata) :: calc
    integer,intent(in) :: tight,maxcycle,maxmicro,nvar
    real(wp),intent(in) :: ethr,gthr
    logical,intent(in) :: linear,wr
    !> LOCAL
    integer :: iupdat,engine
    real(wp) :: hlow,hmax,maxdispl,s6
    type(mhparam) :: mhset
    !> FORMATS
    character(len=*),parameter :: scifmt = &
                                  '(10x,"│",3x,a,e21.7,1x,a,1x,"│")'
    character(len=*),parameter :: dblfmt = &
                                  '(10x,"│",3x,a,f18.7,5x,a,1x,"│")'
    character(len=*),parameter :: intfmt = &
                                  '(10x,"│",3x,a,i18,      10x,"│")'
    character(len=*),parameter :: chrfmt = &
                                  '(10x,"│",3x,a,a18,      10x,"│")'

!>--- set params
    engine = calc%opt_engine
    iupdat = calc%iupdat
    hlow = calc%hlow_opt !> 0.01 in ancopt, 0.002 too small
    hmax = calc%hmax_opt
    maxdispl = calc%maxdispl_opt
    s6 = mhset%s6 !> slightly better than 30 for various proteins

!>--- print a summary of settings, if desired
    write (*,'(/,10x,"┍",49("━"),"┑")')
    write (*,'(10x,"│",11x,a,11x,"│")') "GEOMETRY OPTIMIZATION SETUP"
    write (*,'(10x,"┝",49("━"),"┥")')
    select case (engine)
    case (0)
      write (*,chrfmt) "algorithm         ","          ANCOPT"
    case (1)
      write (*,chrfmt) "algorithm         ","          L-BFGS"
    case (2)
      write (*,chrfmt) "algorithm         ","Rational Function"
    case (-1)
      write (*,chrfmt) "algorithm         ","Gradient Descent"
    end select

    select case (tight)
    case (-3)
      write (*,chrfmt) "optimization level",'crude (-3)'
    case (-2)
      write (*,chrfmt) "optimization level",'vloose (-2)'
    case (-1)
      write (*,chrfmt) "optimization level",'loose (-1)'
    case (0)
      write (*,chrfmt) "optimization level",'normal (0)'
    case (1)
      write (*,chrfmt) "optimization level",'tight (1)'
    case (2)
      write (*,chrfmt) "optimization level",'vtight (2)'
    case (3)
      write (*,chrfmt) "optimization level",'extreme (3)'
    case default
      write (*,chrfmt) "optimization level",'custom'
    end select
    write (*,intfmt) "max. optcycles    ",maxcycle

    select case (engine)
    case (0)
      write (*,intfmt) "ANC micro-cycles  ",maxmicro
    end select

    write (*,intfmt) "degrees of freedom",nvar
    !if (modef > 0) then
    !  write (*,intfmt) "# mode follow     ",modef
    !end if
    write (*,'(10x,"├",49("─"),"┤")')
    select case (engine)
    case (0,2)
      if (calc%exact_rf) then
        write (*,chrfmt) "RF solver         ","spevx"
      else
        write (*,chrfmt) "RF solver         ","davidson"
      end if
    end select
    if (engine .ne. -1) then
      select case (iupdat)
      case (0)
        write (*,chrfmt) "Hessian update    ","bfgs"
      case (1)
        write (*,chrfmt) "Hessian update    ","powell"
      case (2)
        write (*,chrfmt) "Hessian update    ","sr1"
      case (3)
        write (*,chrfmt) "Hessian update    ","bofill"
      case (4)
        write (*,chrfmt) "Hessian update    ","schlegel"
      end select
    end if
    write (*,chrfmt) "write crestopt.log",bool2string(wr)
    if (linear) then
      write (*,chrfmt) "linear (good luck)",bool2string(linear)
    else
      write (*,chrfmt) "linear?           ",bool2string(linear)
    end if
    write (*,scifmt) "energy convergence",ethr,"Eh   "

    select case (engine)
    case (0)
      write (*,scifmt) "grad. convergence ",gthr,"Eh/α "
    case default
      write (*,scifmt) "grad. convergence ",gthr,"Eh/a0"
    end select

    select case (engine)
    case (0)
      write (*,dblfmt) "maximium RF displ.",maxdispl,"    "
      write (*,scifmt) "Hlow (freq-cutoff)",hlow,"     "
      write (*,dblfmt) "Hmax (freq-cutoff)",hmax,"    "
      write (*,dblfmt) "S6 in model hess. ",s6,"    "
    case (2)
      write (*,dblfmt) "maximium RF displ.",maxdispl,"    "
      write (*,'(10x,"│",3x,a,e21.7,1x,a,"│")') &
      & "Hessian guess     ",calc%hguess,"Eh/a0²"
    case (-1)
      write (*,dblfmt) "max step size     ",maxdispl,"    "
    end select
    write (*,'(10x,"└",49("─"),"┘")')

  end subroutine print_optsummary

!========================================================================================!

  subroutine print_convd(econverged,gconverged)
    logical,intent(in) :: econverged,gconverged
    write (*,'(3x,"converged δE/grad :",1x,a5," /",1x,a)') &
    & bool2string(econverged),bool2string(gconverged)
  end subroutine print_convd

!========================================================================================!
!========================================================================================!
end module optimize_utils
