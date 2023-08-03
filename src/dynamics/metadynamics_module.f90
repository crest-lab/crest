!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2023 Philipp Pracht
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

module metadynamics_module

  use crest_parameters
  use ls_rmsd
  use strucrd

  implicit none

  !======================================================================================!
  !--- private module variables and parameters
  private
  integer,parameter,public :: cv_std_mtd = 1
  integer,parameter,public :: cv_rmsd = 2
  integer,parameter :: damp_heaviside = 3
  integer,parameter :: damp_heaviside_cv = 4
  integer,parameter,public :: cv_rmsd_static = 5

  !======================================================================================!
  !data object that contains settings and trackers for a single MTD potential
  type :: mtdpot

    integer :: mtdtype = 0
    integer :: nmax = 0
    integer :: ncur = 0

    real(wp) :: kpush = 1.0_wp
    real(wp) :: alpha = 1.0_wp

    !>--- regular MTD, some CV
    real(wp),allocatable :: cv(:) !list of cv values at each timestep
    real(wp),allocatable :: cvgrd(:,:) !Cartesian gradient of CV at step t

    !>--- RMSD MTD
    integer :: cvdump = 0 !xyz dump counter
    real(wp) :: cvdump_fs = 0.0_wp !xyz dump frequency (in fs)
    integer :: cvdumpstep = 0 !xyz dump frequency (in MD steps)
    integer :: maxsave = 0
    character(len=:),allocatable :: biasfile !> specify a file from which the bias is obtained
    logical,allocatable :: atinclude(:) !specify atoms to include in RMSD potentail
    real(wp),allocatable :: cvxyz(:,:,:) !ensemble of CV structures to calculate RMSD from

    !>--- damping of the MTD potential
    integer :: damptype = 0
    real(wp) :: ramp = -1.0_wp
    real(wp) :: damp = 1.0_wp
    real(wp),allocatable :: damping(:)  !input for snapshot-specific damping

  contains
    procedure :: deallocate => mtd_deallocate
    procedure :: info => mtd_info
  end type mtdpot

  public :: mtdpot
  public :: mtd_ini
  public :: cv_dump
  public :: calc_mtd

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine mtd_ini(mol,pot,tstep,mdlength,pr)
!*************************************
!* subroutine mtd_ini
!* initialize metadynamics settings
!*************************************
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(in) :: tstep !> MD timestep in fs
    real(wp),intent(in) :: mdlength !> MD length in ps
    logical,intent(in) :: pr
    real(wp) :: dum1,dum2
    integer :: i,j
    integer :: idum1,idum2,nall,nat
    logical :: ex
    integer,allocatable :: at(:)

    if (pr) then
      write (stdout,'(">--- metadynamics parameter ---")')
    end if

    dum1 = anint((mdlength*1000.0_wp)/tstep)
    idum1 = nint(dum1)

    select case (pot%mtdtype)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    case (cv_std_mtd) !>--- "standard" MTD
      pot%nmax = idum1
      allocate (pot%cv(idum1),source=0.0_wp)
      allocate (pot%cvgrd(3,mol%nat),source=0.0_wp)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    case (cv_rmsd) !>--- RMSD MTD
      if (pot%cvdump_fs <= 0.0_wp) return !> structure dumpstep in fs must be given
      dum2 = max(1.0_wp, (pot%cvdump_fs/tstep))
      dum2 = anint(dum2)
      pot%cvdumpstep = nint(dum2) !> bias structure dump step in MD steps
      dum1 = dum1/dum2
      dum1 = floor(dum1)
      pot%nmax = nint(dum1) !> max number of bias structure dump
      if (pot%maxsave == 0) pot%maxsave = nint(dum1)
      if (allocated(pot%cvxyz)) deallocate (pot%cvxyz)
      allocate (pot%cvxyz(3,mol%nat,pot%maxsave),source=0.0_wp)
      !>--- automatic ramp parameter (acounted for both different MD time steps and CV dump steps)
      !> (should yield damp≈0.5 for cvdumpstep/2, but is at least 0.03 as in xtb)
      if (pot%ramp <= 0.0_wp) then !> only if not set by the user
        pot%ramp = log(3.0_wp)/(0.5_wp*dum2)
        pot%ramp = max(pot%ramp,0.03_wp)
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    case (cv_rmsd_static) !> static RMSD MTD / umbrella sampling with RMSD pot
      pot%cvdumpstep = huge(idum1) !> we won't modify the potential, so set the dumstep to "infinity"
      pot%cvdump_fs = huge(dum1)   !> same for _fs version
      if (allocated(pot%cvxyz)) then !> if structures were already loaded, just determine the rest
        nat = size(pot%cvxyz,2)
        nall = size(pot%cvxyz,3)
      else if (allocated(pot%biasfile)) then !> else try to read from file
        inquire (file=pot%biasfile,exist=ex)
        if (.not.ex) return !> if the file is absent, we return
        call rdensembleparam(pot%biasfile,nat,nall)
        if (nat .ne. 0.and.nall .ne. 0) then
          allocate (pot%cvxyz(3,nat,nall),source=0.0_wp)
          allocate (at(nat),source=0)
          call rdensemble(pot%biasfile,nat,nall,at,pot%cvxyz)
          deallocate (at)
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
          !>--- Important: if we read here we must convert to Bohrs
          pot%cvxyz = pot%cvxyz*aatoau
          !>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
        end if
      else  !> if both failed we must return (the potential is not set up)
        return
      end if
      if (nat .ne. mol%nat) then !> can't do that! something is wrong
        if (allocated(pot%cvxyz)) deallocate (pot%cvxyz)
        pot%mtdtype = 0
        write (*,'(1x,a)') '*WARNING* static metadynamics setup failed! Mismatch of #atoms'
        !return
        error stop
      end if
      !> for safety, perturb all bias slightly (so the potential won't explode)
      do i = 1,nall
        call rmsdcv_perturb(nat,pot%cvxyz(:,:,i))
      end do
      pot%ncur = nall    !> will not change
      pot%maxsave = nall !> won't change either
      if (pot%ramp <= 0.0_wp) then          !> only if not set by the user
        pot%ramp = (tstep/5.0_wp)*0.015_wp  !> a default derived from the entropy mode at GFN2-xTB
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    case default
      return

    end select

    !>--- printout
    if (pr) then
      call pot%info(stdout)
    end if


    return
  end subroutine mtd_ini

!========================================================================================!
  subroutine mtd_deallocate(self)
!**********************************************
!* subroutine mtd_deallocate
!* type internal procedure to deallocate data
!**********************************************
    class(mtdpot) :: self
    if (allocated(self%cvxyz)) deallocate (self%cvxyz)
    if (allocated(self%atinclude)) deallocate (self%atinclude)
    if (allocated(self%cv)) deallocate (self%cv)
    if (allocated(self%cvgrd)) deallocate (self%cvgrd)

    self%mtdtype = 0
    self%nmax = 0
    self%ncur = 0
    self%kpush = 1.0_wp
    self%alpha = 1.0_wp
    self%cvdump = 0 !xyz dump counter
    self%cvdump_fs = 0.0_wp !xyz dump frequency (in fs)
    self%cvdumpstep = 0 !xyz dump frequency (in MD steps)
    self%maxsave = 0
    self%damptype = 0
    self%ramp = -1.0_wp
    self%damp = 1.0_wp

    return
  end subroutine mtd_deallocate

!========================================================================================!
  subroutine mtd_info(self,iunit)
!**********************************************
!* subroutine mtd_info
!* print information about the MTD potential
!**********************************************
    class(mtdpot) :: self
    integer,intent(in) :: iunit

    !write (iunit,'(" --- metadynamics parameter ---")')
    select case (self%mtdtype)
    case (cv_std_mtd)
      write (iunit,'(" MTD/CV type   :",1x,a)') 'standard'
    case (cv_rmsd)
      write (*,'(" MTD/CV type   :",1x,a)') 'RMSD bias'
    case (cv_rmsd_static)
      write (iunit,'(" MTD/CV type   :",1x,a)') 'RMSD bias (static)'
    end select
    write (iunit,'(" kpush /Eh     :",f10.4)') self%kpush
    write (iunit,'(" alpha /bohr⁻² :",f10.4)') self%alpha

    select case (self%mtdtype)
    case (cv_rmsd)
      write (iunit,'(" ramp          :",f10.4,1x,i0)') self%ramp,check_dump_steps_rmsd(self)
      write (iunit,'(" dump/fs       :",f10.4,1x,i0 )') self%cvdump_fs,self%cvdumpstep
      write (iunit,'(" # CVs (max)   :",i10 )') self%maxsave
    case (cv_rmsd_static)
      if (allocated(self%biasfile)) write (iunit,'(" reading from  :",1x,a)') self%biasfile
      write (iunit,'(" ramp (adjust.):",f10.4,1x,i0)') self%ramp,check_dump_steps_rmsd(self)
      write (iunit,'(" # CVs (loaded):",i10 )') self%maxsave
    end select

    return
  end subroutine mtd_info

!========================================================================================!
  subroutine cv_dump(mol,pot,cv,pr)
!*****************************************************
!* subroutine cv_dump
!* update the list of CVs at the current MD timestep
!*****************************************************
!$  use omp_lib
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(in) :: cv
    logical :: pr

    select case (pot%mtdtype)
    case (cv_std_mtd) !>--- CV update for "standard" MTD
      pot%ncur = pot%ncur+1
      pot%cv(pot%ncur) = cv

    case (cv_rmsd) !>--- structure mapping for RMSD MTD
      pot%cvdump = pot%cvdump+1  !> cvdump counts the MD step since the last CV was added
      if (pot%cvdump == pot%cvdumpstep) then !> the MTD tracks when it needs to be updated
        pot%cvdump = 0  !> reset if new CV is added
        pot%ncur = pot%ncur+1
        pot%cvxyz(:,:,pot%ncur) = mol%xyz(:,:)
        if (pot%ncur == 1) then
          !>--- The first one should be sligthly distorted
          call rmsdcv_perturb(mol%nat,pot%cvxyz(:,:,pot%ncur))
        end if
        if (pr) then
          write (*,'(2x,"adding snapshot to metadynamics bias, now at ",i0," CVs")') pot%ncur
        end if
      end if

    case (cv_rmsd_static)
      pot%cvdump = pot%cvdump+1 !> the cvdump is equal to the MD step
      !> no further update necessary

    case default
      return
    end select

    return

  end subroutine cv_dump

!=========================================================================================!
  subroutine rmsdcv_perturb(nat,xyz)
!************************************************
!* Slightly perturb a given geometry for RMSD CV
!* to avoid singularities if the CV is exactly the
!* current structure
!*************************************************
    implicit none
    integer,intent(in) :: nat
    real(wp),intent(inout) :: xyz(3,nat)
    real(wp) :: r(3)
    integer :: i,j
    real(wp),parameter :: tol = 1.0e-8_wp
    real(wp),parameter :: displace = 1.0e-6_wp
    do i = 1,nat
      do
        !> generate a random vector r in [-1,1]
        call random_number(r)
        r = (r-0.5_wp)*2.0_wp
        !> check that displacement is large enough
        if (norm2(r) >= 1e-8_wp) exit
      end do
      !> normalize
      r = r/norm2(r)
      !> displace
      xyz(:,i) = xyz(:,i)+displace*r
    end do
  end subroutine rmsdcv_perturb

!========================================================================================!
  subroutine calc_mtd(mol,pot,emtd,grdmtd)
!*******************************************************
!* subroutine calc_mtd
!* select how the MTD potential is calculated.
!* On Output:
!*             emtd - final MTD energy contribution
!*           mtdgrd - final MTD gradient contribution
!********************************************************
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(out) :: emtd
    real(wp),intent(out) :: grdmtd(3,mol%nat)
    real(wp) :: dum
    emtd = 0.0_wp
    grdmtd = 0.0_wp

    select case (pot%mtdtype)
    case (cv_std_mtd)
      call calc_damp(pot,pot%damptype,0.0_wp)

    case (cv_rmsd)
      dum = float(pot%cvdump)
      call calc_damp(pot,cv_rmsd,dum)
      call calc_rmsd_mtd(mol,pot,emtd,grdmtd)

    case (cv_rmsd_static)
      dum = float(pot%cvdump)
      call calc_damp(pot,cv_rmsd_static,dum)
      call calc_rmsd_mtd(mol,pot,emtd,grdmtd)

    case default
      emtd = 0.0_wp
      grdmtd = 0.0_wp

    end select

    return
  end subroutine calc_mtd

!========================================================================================!
  subroutine calc_damp(pot,dt,x)
!**************************************************
!* subroutine calc_damp
!* calculate a MTD-type-specific damping factor.
!* If/how/where the damping factor is applied
!* depends on the MTD type
!**************************************************
    implicit none
    type(mtdpot) :: pot
    integer :: dt
    real(wp),intent(in) :: x
    real(wp),parameter :: tol = 0.9999_wp
    select case (dt) !>-- select damping parameter calculation
    case (cv_rmsd,cv_rmsd_static)
      pot%damp = (2.0_wp/(1.0_wp+ &
      &       exp(-pot%ramp*x))-1.0_wp)  !> x is pot%cvdump as float

    case (damp_heaviside) !> simple heaviside switch
      pot%damp = sign(0.5_wp,x)+0.5_wp

    case default
      pot%damp = 1.0_wp

    end select

    return
  end subroutine calc_damp

  function check_dump_steps_rmsd(pot) result(steps)
!**************************************************
!* Check the number of MD steps that are affected
!* by the damping
!**************************************************
    implicit none
    type(mtdpot) :: pot
    integer :: steps
    real(wp) :: dum
    real(wp),parameter :: tol = 0.9999_wp
    steps = 0
    do
      steps = steps+1
      dum = float(steps)
      call calc_damp(pot,cv_rmsd,dum)
      if (pot%damp > tol) then
        pot%damp = 0.0_wp
        exit
      end if
    end do
  end function check_dump_steps_rmsd

!========================================================================================!
  subroutine calc_damp2(pot,t,damp)
!*********************************************
!* subroutine calc_damp2
!* damping routine for snapshot-cv-specific
!* damping parameter
!*********************************************
    implicit none
    type(mtdpot) :: pot
    integer :: t !> snapshot
    real(wp) :: damp

    select case (pot%damptype) !>-- select damping parameter calculation
    case (damp_heaviside_cv) !> simple heaviside switch
      damp = sign(0.5_wp,pot%damping(t))+0.5_wp
    case default
      pot%damp = 1.0_wp
    end select

    return
  end subroutine calc_damp2

!========================================================================================!
  subroutine calc_rmsd_mtd(mol,pot,ebias,grdmtd)
!**************************************************************
!* subroutine calc_rmsd_mtd
!* calculate energy and gradient contribution from the RMSD
!* of the current structure (mol) to any structure in a list
!* of documented references.
!* Optionally, atoms for which the RMSD is to be calculated
!* can be specified.
!* Since RMSD calculation can be costly for many structures
!* there is some OMP parallelization going on.
!**************************************************************
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(out) :: ebias
    real(wp),intent(out) :: grdmtd(3,mol%nat)

    real(wp),allocatable :: xyzref(:,:)
    real(wp),allocatable :: xyzcp(:,:)
    real(wp),allocatable :: grad(:,:)
    real(wp) :: U(3,3),x_center(3),y_center(3)
    real(wp) :: rmsdval,E,dEdr

    integer :: i,j,k,l

    ebias = 0.0_wp
    grdmtd = 0.0_wp

    if (pot%ncur < 1) return

    if (.not.allocated(pot%atinclude)) then !>-- include all atoms in RMSD
      allocate (xyzref(3,mol%nat),grad(3,mol%nat),source=0.0_wp)
      !$omp parallel default(none) &
      !$omp shared(pot,mol) &
      !$omp private(grad,xyzref,U,x_center,y_center,rmsdval,E,dEdr) &
      !$omp reduction(+:ebias,grdmtd)
      !$omp do schedule(dynamic)
      do i = 1,pot%ncur
        grad = 0.0_wp
        xyzref = pot%cvxyz(:,:,i)
        call rmsd(mol%nat,mol%xyz,xyzref,1,U,x_center,y_center,rmsdval, &
        &          .true.,grad)
        E = pot%kpush*exp(-pot%alpha*rmsdval**2)
        if (i == pot%ncur.or.pot%mtdtype == cv_rmsd_static) then
          E = E*pot%damp
        end if
        ebias = ebias+E
        dEdr = -2.0_wp*pot%alpha*e*rmsdval
        grdmtd = grdmtd+dEdr*grad
      end do
      !$omp enddo
      !$omp end parallel
      deallocate (grad,xyzref)

    else !>--- use only selected atoms in RMSD
      k = count(pot%atinclude,1)
      if (k < 1) return
      allocate (xyzcp(3,k),xyzref(3,k),grad(3,k),source=0.0_wp)
      !$omp parallel default(none) &
      !$omp shared(pot,mol,k) &
      !$omp private(grad,xyzref,U,x_center,y_center,rmsdval,E,dEdr) &
      !$omp private(xyzcp,j,l) &
      !$omp reduction(+:ebias,grdmtd)
      !$omp do schedule(dynamic)
      do i = 1,pot%ncur
        grad = 0.0_wp
        l = 0
        do j = 1,mol%nat
          if (pot%atinclude(j)) then
            l = l+1
            xyzcp(:,l) = mol%xyz(:,j)
            xyzref(:,l) = pot%cvxyz(:,j,i)
          end if
        end do
        call rmsd(k,xyzcp,xyzref,1,U,x_center,y_center,rmsdval, &
        &          .true.,grad)
        E = pot%kpush*exp(-pot%alpha*rmsdval**2)
        if (i == pot%ncur.or.pot%mtdtype == cv_rmsd_static) then
          E = E*pot%damp
        end if
        ebias = ebias+E
        dEdr = -2.0_wp*pot%alpha*e*rmsdval
        l = 0
        do j = 1,mol%nat
          if (pot%atinclude(j)) then
            l = l+1
            grdmtd(:,j) = grdmtd(:,j)+dEdr*grad(:,l)
          end if
        end do
      end do
      !$omp enddo
      !$omp end parallel
      deallocate (grad,xyzref,xyzcp)
    end if

    return
  end subroutine calc_rmsd_mtd

!========================================================================================!
  subroutine calc_std_mtd(mol,pot,cvt,ebias,grdmtd)
!*******************************************************
!* subroutine calc_std_mtd
!* calculate energy and gradient contribution from the
!* standard MTD formulation (list of CVs)
!*******************************************************
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp) :: cvt  !> value of the CV at the current timestep
    real(wp),intent(out) :: ebias
    real(wp),intent(out) :: grdmtd(3,mol%nat)

    real(wp) :: U(3,3),x_center(3),y_center(3)
    real(wp) :: rmsdval,E,dEdr,dcv,damp2

    integer :: i,j,k,l

    ebias = 0.0_wp
    grdmtd = 0.0_wp

    if (pot%ncur < 1) return

    !$omp parallel default(none) &
    !$omp shared(pot,mol,cvt) &
    !$omp private(U,x_center,y_center,dcv,damp2,E,dEdr) &
    !$omp reduction(+:ebias,grdmtd)
    !$omp do schedule(dynamic)
    do i = 1,pot%ncur
      dcv = cvt-pot%cv(i)
      E = pot%kpush*exp(-pot%alpha*dcv**2)  !> Gaussian shaped potential
      E = E*pot%damp
      call calc_damp2(pot,i,damp2)
      E = E*damp2
      ebias = ebias+E
      dEdr = -2.0_wp*pot%alpha*e*dcv
      grdmtd = grdmtd+dEdr*pot%cvgrd
    end do
    !$omp enddo
    !$omp end parallel

    return

  end subroutine calc_std_mtd

!========================================================================================!
end module metadynamics_module
