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

module metadynamics_module

  use iso_fortran_env,only:wp => real64,error_unit
  use ls_rmsd
  use strucrd

  implicit none

  !======================================================================================!
  !--- private module variables and parameters
  private

  !--- some constants and name mappings
  real(wp),parameter :: bohr = 0.52917726_wp
  real(wp),parameter :: autokcal = 627.509541_wp
  real(wp),parameter :: autoaa = 0.52917726_wp
  real(wp),parameter :: aatoau = 1.0_wp / autoaa
  real(wp),parameter :: amutokg = 1.660539040e-27_wp
  real(wp),parameter :: metokg = 9.10938356e-31_wp
  real(wp),parameter :: kgtome = 1.0_wp / metokg
  real(wp),parameter :: amutoau = amutokg * kgtome
  real(wp),parameter :: fstoau = 41.3413733365614_wp
  real(wp),parameter :: kB = 3.166808578545117e-06_wp !in Eh/K

  integer,parameter :: std_mtd = 1
  integer,parameter :: rmsd_mtd = 2
  integer,parameter :: damp_heaviside = 3
  integer,parameter :: damp_heaviside_cv = 4

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
    logical,allocatable :: atinclude(:) !specify atoms to include in RMSD potentail
    real(wp),allocatable :: cvxyz(:,:,:) !ensemble of CV structures to calculate RMSD from

    !>--- damping of the MTD potential
    integer :: damptype = 0
    real(wp) :: ramp = 0.03_wp
    real(wp) :: damp = 1.0_wp
    real(wp),allocatable :: damping(:)  !input for snapshot-specific damping
    

  contains
    procedure :: deallocate => mtd_deallocate

  end type mtdpot

  public :: mtdpot
  public :: mtd_ini
  public :: cv_dump
  public :: calc_mtd

contains
!========================================================================================!
! subroutine mtd_ini
! initialize mtd settings
!-------------------------!
  subroutine mtd_ini(mol,pot,tstep,mdlength,pr)
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(in) :: tstep !> MD timestep in fs
    real(wp),intent(in) :: mdlength !> MD length in ps
    logical,intent(in) :: pr
    real(wp) :: dum1,dum2
    integer :: idum1,idum2

    if (pr) then
      write (*,'(" --- metadynamics parameter ---")')
      write (*,'(" kpush   :",f9.3)') pot%kpush
      write (*,'(" alpha   :",f9.3)') pot%alpha
    end if

    dum1 = anint((mdlength * 1000.0_wp) / tstep)
    idum1 =nint(dum1)

    select case (pot%mtdtype)
    case (std_mtd) !>--- "standard" MTD
            pot%nmax = idum1
      allocate(pot%cv(idum1),source=0.0_wp)
      allocate(pot%cvgrd(3,mol%nat),source=0.0_wp)
    case (rmsd_mtd) !>--- RMSD MTD
      if (pot%cvdump_fs <= 0.0_wp) return !> structure dumpstep in fs must be given
      dum2 = max(1.0_wp, (pot%cvdump_fs / tstep))
      dum2 = anint(dum2)
      pot%cvdumpstep = nint(dum2) !> bias structure dump step in MD steps
      dum1 = dum1 / dum2
      dum1 = floor(dum1)
      pot%nmax = nint(dum1) !> max number of bias structure dump
      if (pot%maxsave == 0) pot%maxsave = nint(dum1)
      allocate (pot%cvxyz(3,mol%nat,pot%maxsave),source=0.0_wp)
      !>--- automatic ramp parameter
      !> (should yield damp≈0.5 for cvdumpstep/2, but is at least 0.03 as in xtb)
      pot%ramp = log(3.0_wp) / (0.5_wp * dum2)
      pot%ramp = max(pot%ramp,0.03_wp)
      if (pr) then
        write (*,'(" ramp    :",f9.3)') pot%ramp
        write (*,'(" dump/fs :",f9.3,i9 )') pot%cvdump_fs,pot%cvdumpstep
        write (*,'(" # CVs   :",i9 )') pot%maxsave
      end if
    case default
      return
    end select

    return
  end subroutine mtd_ini

!========================================================================================!
! subroutine mtd_deallocate
! type internal procedure to deallocate data
  subroutine mtd_deallocate(self)
    class(mtdpot) :: self
    if (allocated(self%cvxyz)) deallocate (self%cvxyz)
    if (allocated(self%atinclude)) deallocate (self%atinclude)
    if(allocated(self%cv)) deallocate(self%cv)
    if(allocated(self%cvgrd)) deallocate(self%cvgrd)

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
    self%ramp = 0.03_wp
    self%damp = 1.0_wp

    return
  end subroutine mtd_deallocate

!========================================================================================!
! subroutine cv_dump
! update the list of CVs at the current MD timestep
!---------------------------------------------------!
  subroutine cv_dump(mol,pot,cv,pr)
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(in) :: cv
    logical :: pr

    select case (pot%mtdtype)
    case (std_mtd) !>--- CV update for "standard" MTD
      pot%ncur = pot%ncur + 1
      pot%cv(pot%ncur) = cv
    case (rmsd_mtd) !>--- structure mapping for RMSD MTD
      pot%cvdump = pot%cvdump + 1
      if (pot%cvdump == pot%cvdumpstep) then
        pot%cvdump = 0
        pot%ncur = pot%ncur + 1
        pot%cvxyz(:,:,pot%ncur) = mol%xyz(:,:)
        if(pr)then
        write (*,'(2x,"adding snapshot to metadynamics bias, now at ",i0," CVs")')pot%ncur
        endif
      end if
    case default
      return
    end select

    return
  end subroutine cv_dump

!========================================================================================!
! subroutine calc_mtd
! select how the MTD potential is calculated.
! On Output:
!             emtd - final MTD energy contribution
!           mtdgrd - final MTD gradient contribution
!----------------------------------------------------!
  subroutine calc_mtd(mol,pot,emtd,grdmtd)
    implicit none
    type(coord) :: mol
    type(mtdpot) :: pot
    real(wp),intent(out) :: emtd
    real(wp),intent(out) :: grdmtd(3,mol%nat)

    emtd = 0.0_wp
    grdmtd = 0.0_wp

    select case (pot%mtdtype)
    case (std_mtd)
      call calc_damp(pot,pot%damptype,0.0_wp)

    case (rmsd_mtd)
      call calc_damp(pot,rmsd_mtd,0.0_wp)
      call calc_rmsd_mtd(mol,pot,emtd,grdmtd)
    case default
      emtd = 0.0_wp
      grdmtd = 0.0_wp
    end select

    return
  end subroutine calc_mtd

!========================================================================================!
! subroutine calc_damp
! calculate a MTD-type-specific damping factor.
! If/how/where the damping factor is applied
! depends on the MTD type
!-----------------------------------------------!
  subroutine calc_damp(pot,dt,x)
    implicit none
    type(mtdpot) :: pot
    integer :: dt
    real(wp) :: x

    select case (dt) !>-- select damping parameter calculation
    case (rmsd_mtd) 
      pot%damp = (2.0_wp / (1.0_wp + &
      &       exp(-pot%ramp * float(pot%cvdump))) - 1.0_wp)
    case( damp_heaviside ) !> simple heaviside switch
      pot%damp = sign(0.5_wp,x) + 0.5_wp
    case default
      pot%damp = 1.0_wp
    end select

    return
  end subroutine calc_damp

!========================================================================================!
! subroutine calc_damp2
! damping routine for snapshot-cv-specific 
! damping parameter
!-----------------------------------------------!
  subroutine calc_damp2(pot,t,damp)
    implicit none
    type(mtdpot) :: pot
    integer :: t !> snapshot
    real(wp) :: damp


    select case (pot%damptype) !>-- select damping parameter calculation
    case( damp_heaviside_cv ) !> simple heaviside switch
      damp = sign(0.5_wp,pot%damping(t)) + 0.5_wp
    case default
      pot%damp = 1.0_wp
    end select

    return
  end subroutine calc_damp2

!========================================================================================!
! subroutine calc_rmsd_mtd
! calculate energy and gradient contribution from the RMSD
! of the current structure (mol) to any structure in a list
! of documented references.
! Optionally, atoms for which the RMSD is to be calculated
! can be specified.
! Since RMSD calculation can be costly for many structures
! there is some OMP parallelization going on.
!-----------------------------------------------------------!
  subroutine calc_rmsd_mtd(mol,pot,ebias,grdmtd)
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

    if (.not. allocated(pot%atinclude)) then !>-- include all atoms in RMSD
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
        E = pot%kpush * exp(-pot%alpha * rmsdval**2)
        if (i == pot%ncur) then
          E = E * pot%damp
        end if
        ebias = ebias + E
        dEdr = -2.0_wp * pot%alpha * e * rmsdval
        grdmtd = grdmtd + dEdr * grad
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
            l = l + 1
            xyzcp(:,l) = mol%xyz(:,j)
            xyzref(:,l) = pot%cvxyz(:,j,i)
          end if
        end do
        call rmsd(k,xyzcp,xyzref,1,U,x_center,y_center,rmsdval, &
        &          .true.,grad)
        E = pot%kpush * exp(-pot%alpha * rmsdval**2)
        if (i == pot%ncur) then
          E = E * pot%damp
        end if
        ebias = ebias + E
        dEdr = -2.0_wp * pot%alpha * e * rmsdval
        l = 0
        do j = 1,mol%nat
          if (pot%atinclude(j)) then
            l = l + 1
            grdmtd(:,j) = grdmtd(:,j) + dEdr * grad(:,l)
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
! subroutine calc_std_mtd
! calculate energy and gradient contribution from the 
! standard MTD formulation (list of CVs)
!-----------------------------------------------------------!
  subroutine calc_std_mtd(mol,pot,cvt,ebias,grdmtd)
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
        dcv = cvt - pot%cv(i) 
        E = pot%kpush * exp(-pot%alpha * dcv**2)  !> Gaussian shaped potential
        E = E * pot%damp
        call calc_damp2(pot,i,damp2)
        E = E * damp2
        ebias = ebias + E
        dEdr = -2.0_wp * pot%alpha * e * dcv
        grdmtd = grdmtd + dEdr * pot%cvgrd
      end do
      !$omp enddo
      !$omp end parallel

    return

  end subroutine calc_std_mtd

!========================================================================================!
end module metadynamics_module
