!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2026 Philipp Pracht
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

!> Biased single-point Hessian (SPH, xtb --bhess) on top of the CREST calculator.
!> Ported from xtb (src/getkopt.f90, src/hessian.F90); Spicher & Grimme,
!> J. Chem. Theory Comput. 2021, 17, 1701.

module bhess_module
  use crest_parameters
  use crest_calculator
  use strucrd
  use optimize_module
  use irmsd_module,only:rmsd
  use axis_module,only:axis
  use thermochem_module,only:prj_mw_hess,mass_weight_hess
  implicit none
  private

  public :: bhess_frequencies

  real(wp),parameter :: bhess_target_rmsd = 0.1_wp  !> Angstroem
  real(wp),parameter :: bhess_alpha = 1.0_wp        !> Bohr^-2
  real(wp),parameter :: bhess_accuracy = 0.3_wp     !> SCC accuracy of the Hessian

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine bhess_frequencies(mol,calc,molh,freq,kpush,io,pr)
    !********************************************************************
    !* RMSD-biased Hessian of mol and its rescaled frequencies.
    !* A Gaussian RMSD bias -k*exp(-a*RMSD^2) centered on the input
    !* keeps a biased optimization within 0.1 A of the input (k by
    !* bisection); the Hessian is taken there and the frequencies are
    !* scaled back by the ratio of the unbiased to the biased curvature.
    !*
    !* mol   : input geometry (Bohr), also the bias reference
    !* calc  : unbiased calculator (not modified)
    !* molh  : geometry of the Hessian (Bohr)
    !* freq  : scaled frequencies in cm-1, 3*nat entries
    !* kpush : total bias strength (Eh), 0 if no bias was needed
    !********************************************************************
    type(coord),intent(in) :: mol
    type(calcdata),intent(in) :: calc
    type(coord),intent(out) :: molh
    real(wp),intent(out) :: freq(3*mol%nat)
    real(wp),intent(out) :: kpush
    integer,intent(out) :: io
    logical,intent(in) :: pr

    type(calcdata),target :: calcb,calcp
    type(calculation_settings) :: pen
    real(wp),allocatable :: h(:,:),hbias(:,:),htb(:,:),eig(:),scal(:),fctb(:),fcb(:)
    real(wp) :: rabc(3),avmom
    integer :: nat3,ip,i
    logical :: linear

    io = 0
    nat3 = 3*mol%nat
    freq = 0.0_wp

!>--- calculator with the RMSD bias as additional level
    call calcb%copy(calc)
    calcb%pr_energies = .false.
    pen%id = jobtype%penalty
    pen%penalty%alpha = bhess_alpha
    pen%penalty%kpush = 0.0_wp
    allocate (pen%penalty%biastmp(1))
    pen%penalty%biastmp(1) = mol
    call calcb%add(pen)
    ip = calcb%ncalculations
    calcb%calcs(ip)%penalty%biaslist => calcb%calcs(ip)%penalty%biastmp

    call bhess_kopt(mol,calcb,ip,molh,kpush,io,pr)
    if (io /= 0) return

!>--- biased Hessian at the final geometry, and the bias alone
    do i = 1,ip-1
      calcb%calcs(i)%accuracy = bhess_accuracy
    end do
    allocate (h(nat3,nat3),hbias(nat3,nat3),source=0.0_wp)
    call numhess1(molh%nat,molh%at,molh%xyz,calcb,h,io)
    if (io /= 0) return
    call calcp%add(calcb%calcs(ip))
    calcp%calcs(1)%penalty%biaslist => calcp%calcs(1)%penalty%biastmp
    call numhess1(molh%nat,molh%at,molh%xyz,calcp,hbias,io)
    if (io /= 0) return

!>--- project translations/rotations (not for linear molecules) and mass-weight
    call axis(molh%nat,molh%at,molh%xyz*autoaa,rabc,avmom)
    linear = rabc(1) < 1.0e-10_wp
    if (linear) then
      call mass_weight_hess(molh%nat,molh%at,nat3,h)
      call mass_weight_hess(molh%nat,molh%at,nat3,hbias)
    else
      call prj_mw_hess(molh%nat,molh%at,nat3,molh%xyz,h)
      call prj_mw_hess(molh%nat,molh%at,nat3,molh%xyz,hbias)
    end if
    htb = h-hbias

    allocate (eig(nat3),scal(nat3),fctb(nat3),fcb(nat3))
    call bhess_diag(nat3,h,eig,io)
    if (io /= 0) return
    call rescale_freq(nat3,htb,h,hbias,eig,fctb,fcb,scal)
    do i = 1,nat3
      freq(i) = autorcm*sign(sqrt(abs(eig(i))),eig(i))*scal(i)
    end do

    if (pr) then
      write (stdout,'(1x,a,f10.5,a)') 'kpush   :',kpush,' Eh'
      write (stdout,'(1x,a,f10.5,a)') 'alpha   :',bhess_alpha,' Bohr^-2'
    end if
  end subroutine bhess_frequencies

!========================================================================================!

  subroutine bhess_kopt(mol,calcb,ip,molh,kpush,io,pr)
    !********************************************************************
    !* Bisection for the bias strength that keeps the biased optimum
    !* within the target RMSD of the input (xtb get_kopt). molh is the
    !* optimized geometry of the last trial, i.e. at the returned kpush.
    !********************************************************************
    type(coord),intent(in) :: mol
    type(calcdata),intent(inout) :: calcb
    integer,intent(in) :: ip
    type(coord),intent(out) :: molh
    real(wp),intent(out) :: kpush
    integer,intent(out) :: io
    logical,intent(in) :: pr

    real(wp) :: ax,cx,kbias,rmsdval
    integer :: i

    ax = 0.0_wp
    cx = 0.005_wp*(real(mol%nat,wp)/(bhess_target_rmsd+0.01_wp))
    kbias = ax
    call bhess_trial(mol,calcb,ip,kbias,molh,rmsdval,io)
    if (io /= 0) return
    if (pr) then
      write (stdout,'(1x,a,f9.6)') 'target rmsd / A        ',bhess_target_rmsd
      write (stdout,'(1x,a,f9.6)') 'unbiased initial rmsd  ',rmsdval
    end if
    if (rmsdval > bhess_target_rmsd) then
      if (pr) write (stdout,'(1x,a)') 'iter. min.        max.        rmsd       kpush'
      do i = 1,10
        kbias = 0.5_wp*(cx+ax)
        call bhess_trial(mol,calcb,ip,kbias,molh,rmsdval,io)
        if (io /= 0) return
        if (pr) write (stdout,'(1x,i2,4f12.6)') i,ax,cx,rmsdval,-kbias
        if (rmsdval > bhess_target_rmsd) then
          ax = kbias
        else
          cx = kbias
        end if
        if (abs(rmsdval-bhess_target_rmsd) < 1.0e-3_wp) exit
      end do
    end if
    kpush = kbias
    if (pr) write (stdout,'(1x,a,f9.6)') 'final kpush: ',-kpush
  end subroutine bhess_kopt

!========================================================================================!

  subroutine bhess_trial(mol,calcb,ip,kbias,molopt,rmsdval,io)
    !********************************************************************
    !* Optimize mol with the bias -kbias*exp(-a*RMSD^2) and return the
    !* RMSD (A) of the optimum to the input
    !********************************************************************
    type(coord),intent(in) :: mol
    type(calcdata),intent(inout) :: calcb
    integer,intent(in) :: ip
    real(wp),intent(in) :: kbias
    type(coord),intent(out) :: molopt
    real(wp),intent(out) :: rmsdval
    integer,intent(out) :: io
    real(wp),allocatable :: grd(:,:)
    real(wp) :: etot

!>--- the penalty level scales kpush with nat
    calcb%calcs(ip)%penalty%kpush = -kbias/real(mol%nat,wp)
    allocate (grd(3,mol%nat),source=0.0_wp)
    call optimize_geometry(mol,molopt,calcb,etot,grd,.false.,.false.,io)
    if (io /= 0) return
    rmsdval = rmsd(molopt,mol)*autoaa
  end subroutine bhess_trial

!========================================================================================!

  subroutine bhess_diag(n,h,eig,io)
    !**********************************************
    !* eigenvalues and eigenvectors (in h) of h
    !**********************************************
    integer,intent(in) :: n
    real(wp),intent(inout) :: h(n,n)
    real(wp),intent(out) :: eig(n)
    integer,intent(out) :: io
    real(wp),allocatable :: work(:)
    integer,allocatable :: iwork(:)
    external :: dsyevd
    allocate (work(1+6*n+2*n**2),iwork(3+5*n))
    call dsyevd('V','U',n,h,n,eig,work,size(work),iwork,size(iwork),io)
  end subroutine bhess_diag

!========================================================================================!

  subroutine rescale_freq(n3,htb,hess,hbias,eig,fc_tb,fc_bias,freq_scal)
    !********************************************************************
    !* Frequency scaling of the biased Hessian (xtb rescale_freq):
    !* per mode, the curvatures of the unbiased (htb) and the bias
    !* (hbias) parts along the eigenvector give the scaling factor.
    !* Modes with |eigenvalue| <= 1e-6 au stay unscaled.
    !********************************************************************
    integer,intent(in) :: n3
    real(wp),intent(in) :: htb(n3,n3),hess(n3,n3),hbias(n3,n3),eig(n3)
    real(wp),intent(out) :: fc_tb(n3),fc_bias(n3),freq_scal(n3)
    real(wp),parameter :: alp1 = 1.27_wp,alp2 = 1.5e-4_wp
    integer :: j

    do j = 1,n3
      fc_tb(j) = dot_product(hess(:,j),matmul(htb,hess(:,j)))
      fc_bias(j) = dot_product(hess(:,j),matmul(hbias,hess(:,j)))
      if (abs(eig(j)) > 1.0e-6_wp) then
        freq_scal(j) = sqrt((fc_tb(j)+alp2)/((fc_tb(j)+alp2)+alp1*fc_bias(j)))
        if (fc_tb(j) < 0.0_wp.and.fc_bias(j) /= 0.0_wp) then
          freq_scal(j) = -sqrt((abs(fc_tb(j))+alp2)/((abs(fc_tb(j))+alp2)+alp1*fc_bias(j)))
        end if
      else
        freq_scal(j) = 1.0_wp
      end if
    end do
  end subroutine rescale_freq

!========================================================================================!
!========================================================================================!
end module bhess_module
