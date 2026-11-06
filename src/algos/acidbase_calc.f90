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

!> Calculator replacements for the xtb calls of the pKa and acid/base tools

module acidbase_calc
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use bhess_module,only:bhess_frequencies
  use thermochem_module,only:calcthermo
  implicit none
  private

  public :: pka_set_charge
  public :: pka_opt_calc
  public :: ab_calc

!>--- thermo settings of the xtb G(RRHO) the pKa FER was fitted with
  real(wp),parameter :: ab_temp = 298.15_wp
  real(wp),parameter :: ab_ithr = -20.0_wp
  real(wp),parameter :: ab_sthr = 50.0_wp

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine pka_set_charge(env,chrg)
    !*****************************************************************
    !* Set the global charge and rebuild env%calc with it, dropping
    !* any model setup of the previous system
    !*****************************************************************
    type(systemdata),intent(inout) :: env
    integer,intent(in) :: chrg
    type(calcdata) :: tmp
    integer :: i
    env%chrg = chrg
    call tmp%copy(env%calc)
    do i = 1,tmp%ncalculations
      tmp%calcs(i)%chrg = chrg
    end do
    call env%calc%copy(tmp)
  end subroutine pka_set_charge

!========================================================================================!

  subroutine pka_opt_calc(env,fname,io)
    !*****************************************************************
    !* Tight (vtight) optimization of the structure on the xyz file
    !* fname with env%calc; the file is overwritten with the result
    !*****************************************************************
    type(systemdata),intent(inout) :: env
    character(len=*),intent(in) :: fname
    integer,intent(out) :: io
    type(calcdata) :: calc
    type(coord) :: mol,molnew
    real(wp),allocatable :: grd(:,:)
    real(wp) :: etot

    call mol%open(fname)
    call calc%copy(env%calc)
    calc%optlev = 2
    calc%maxcycle = 0
    allocate (grd(3,mol%nat),source=0.0_wp)
    call optimize_geometry(mol,molnew,calc,etot,grd,.false.,.false.,io)
    if (io /= 0) return
    molnew%energy = etot
    molnew%wrextxyz = .false.
    call molnew%write(fname)
  end subroutine pka_opt_calc

!========================================================================================!

  subroutine ab_calc(env,mol,chrg,bhess,etot,gsolv,grrho,wbo,q,io)
    !*****************************************************************
    !* Energies and electronic data of an acid or base with env%calc
    !* (xtb --wbo [--bhess] replacement).
    !*
    !* etot  : total energy (incl. solvation) at the input geometry,
    !*         for bhess at the geometry of the Hessian, without bias
    !* gsolv : E(solution) - E(gas) at the input geometry
    !* grrho : G(RRHO) at 298.15 K from the biased Hessian (bhess)
    !* wbo,q : Wiberg bond orders and atomic charges, taken at the
    !*         geometry of the Hessian for bhess (as in xtb)
    !*****************************************************************
    type(systemdata),intent(inout) :: env
    type(coord),intent(in) :: mol
    integer,intent(in) :: chrg
    logical,intent(in) :: bhess
    real(wp),intent(out) :: etot,gsolv,grrho
    real(wp),intent(out) :: wbo(mol%nat,mol%nat),q(mol%nat)
    integer,intent(out) :: io
    type(calcdata) :: calc,gcalc
    type(coord) :: molh,mold
    real(wp),allocatable :: grd(:,:),freq(:)
    real(wp) :: egas,kpush,temps(1),et(1),ht(1),gt(1),st(1)
    integer :: i

    etot = 0.0_wp
    gsolv = 0.0_wp
    grrho = 0.0_wp
    wbo = 0.0_wp
    q = 0.0_wp
    call calc%copy(env%calc)
    calc%pr_energies = .false.
    do i = 1,calc%ncalculations
      calc%calcs(i)%chrg = chrg
      calc%calcs(i)%rdwbo = .true.
      calc%calcs(i)%rdqat = .true.
    end do
    call gcalc%copy(calc)
    do i = 1,gcalc%ncalculations
      if (allocated(gcalc%calcs(i)%solvmodel)) deallocate (gcalc%calcs(i)%solvmodel)
      if (allocated(gcalc%calcs(i)%solvent)) deallocate (gcalc%calcs(i)%solvent)
    end do

    allocate (grd(3,mol%nat),source=0.0_wp)
    mold = mol
    call engrad(mold,calc,etot,grd,io)
    if (io /= 0) return
    call ab_getwboq(calc,mol%nat,wbo,q)
    call engrad(mold,gcalc,egas,grd,io)
    if (io /= 0) return
    gsolv = etot-egas

    if (bhess) then
      allocate (freq(3*mol%nat))
      call bhess_frequencies(mol,calc,molh,freq,kpush,io,.false.)
      if (io /= 0) return
      temps(1) = ab_temp
      call calcthermo(molh%nat,molh%at,molh%xyz,freq,.false.,ab_ithr,1.0_wp,ab_sthr, &
      &               1,temps,et,ht,gt,st)
      grrho = gt(1)
!>--- xtb reports the energy at the Hessian geometry (there including the bias)
      call engrad(molh,calc,etot,grd,io)
      if (io /= 0) return
      call ab_getwboq(calc,mol%nat,wbo,q)
    end if
  end subroutine ab_calc

!========================================================================================!

  subroutine ab_getwboq(calc,nat,wbo,q)
    !*********************************************
    !* WBO and charges of the first level that
    !* provides them
    !*********************************************
    type(calcdata),intent(in) :: calc
    integer,intent(in) :: nat
    real(wp),intent(out) :: wbo(nat,nat),q(nat)
    integer :: i
    wbo = 0.0_wp
    q = 0.0_wp
    do i = 1,calc%ncalculations
      if (allocated(calc%calcs(i)%wbo)) then
        wbo = calc%calcs(i)%wbo
        exit
      end if
    end do
!>--- xtb's wbo file (read by the legacy route) lists only bond orders > 0.1
    where (wbo <= 0.1_wp) wbo = 0.0_wp
    do i = 1,calc%ncalculations
      if (allocated(calc%calcs(i)%qat)) then
        q = calc%calcs(i)%qat
        exit
      end if
    end do
  end subroutine ab_getwboq

!========================================================================================!
!========================================================================================!
end module acidbase_calc
