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

!> Calculator replacements for the xtb calls of MSREACT

module msreact_calc
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use msmod
  use aiss_iff_prepare,only:lmo_basic_atoms
  implicit none
  private

  !> one distortion: optimize xyz (Angstroem), with a bond constraint
  !> A-B at distance D (Angstroem) and xtb force constant fc if constr
  type :: msjob
    real(wp),allocatable :: xyz(:,:)
    logical  :: constr = .false.
    integer  :: A = 0
    integer  :: B = 0
    real(wp) :: D = 0.0_wp
    real(wp) :: fc = 0.0_wp
  end type msjob
  public :: msjob

  public :: msjob_add
  public :: msreact_sp_calc
  public :: msreact_wbo_calc
  public :: msreact_basic_atoms
  public :: msreact_optjobs

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine msjob_add(jobs,njob,mol,constr,A,B,D,fc)
    !*******************************************
    !* append a distortion job to the list
    !*******************************************
    type(msjob),allocatable,intent(inout) :: jobs(:)
    integer,intent(inout) :: njob
    type(msmol),intent(in) :: mol
    logical,intent(in) :: constr
    integer,intent(in),optional :: A,B
    real(wp),intent(in),optional :: D,fc
    type(msjob),allocatable :: tmp(:)

    if (.not.allocated(jobs)) allocate (jobs(16))
    if (njob >= size(jobs)) then
      allocate (tmp(2*size(jobs)))
      tmp(1:njob) = jobs(1:njob)
      call move_alloc(tmp,jobs)
    end if
    njob = njob+1
    jobs(njob)%xyz = mol%xyz
    jobs(njob)%constr = constr
    if (constr) then
      jobs(njob)%A = A
      jobs(njob)%B = B
      jobs(njob)%D = D
      jobs(njob)%fc = fc
    end if
  end subroutine msjob_add

!========================================================================================!

  subroutine msreact_basecalc(env,etemp,chrg,uhf,calc)
    !*******************************************************
    !* env%calc at the given electronic temperature, charge
    !* and uhf, without constraints
    !*******************************************************
    type(systemdata),intent(inout) :: env
    real(wp),intent(in) :: etemp
    integer,intent(in) :: chrg,uhf
    type(calcdata),intent(out) :: calc
    integer :: i
    call calc%copy(env%calc,ignore_constraints=.true.)
    calc%pr_energies = .false.
    do i = 1,calc%ncalculations
      calc%calcs(i)%etemp = etemp
      calc%calcs(i)%etemp_user_set = .true.
      calc%calcs(i)%chrg = chrg
      calc%calcs(i)%uhf = uhf
      call calc%calcs(i)%sync_multiplicity()
    end do
  end subroutine msreact_basecalc

!========================================================================================!

  subroutine msreact_sp_calc(env,etemp,etot)
    !*******************************************************
    !* singlepoint energy of the input structure
    !*******************************************************
    type(systemdata),intent(inout) :: env
    real(wp),intent(in) :: etemp
    real(wp),intent(out) :: etot
    type(calcdata) :: calc
    type(coord) :: mol
    real(wp),allocatable :: grd(:,:)
    integer :: io

    write (stdout,*)
    call smallhead('Singlepoint Calculation')
    call mol%open(env%inputcoords)
    call msreact_basecalc(env,etemp,env%chrg,env%uhf,calc)
    allocate (grd(3,mol%nat),source=0.0_wp)
    call engrad(mol,calc,etot,grd,io)
    if (io /= 0) then
      write (stdout,*)
      write (stdout,*) ' Initial singlepoint calculation failed!'
      write (stdout,*) ' Please check your input.'
      error stop
    end if
    write (stdout,*) 'Input energy successfully calculated.'
  end subroutine msreact_sp_calc

!========================================================================================!

  subroutine msreact_wbo_calc(env,etemp)
    !*******************************************************
    !* Wiberg bond orders of the input structure, written to
    !* the file "wbo" (entries > 0.1, as xtb) for the topology.
    !* EI: neutral, low-spin molecule; CID: the ion itself
    !*******************************************************
    type(systemdata),intent(inout) :: env
    real(wp),intent(in) :: etemp
    type(calcdata) :: calc
    type(coord) :: mol
    real(wp),allocatable :: grd(:,:),wbo(:,:)
    real(wp) :: etot
    integer :: io,i,j,ich,chrg,uhf

    write (stdout,*)
    call smallhead('WBO Calculation')
    call mol%open(env%inputcoords)
    chrg = env%chrg
    uhf = env%uhf
    if (env%msei) then
      chrg = 0
      call getspin(mol%nat,0,env%inputcoords,uhf)
    end if
    call msreact_basecalc(env,etemp,chrg,uhf,calc)
    do i = 1,calc%ncalculations
      calc%calcs(i)%rdwbo = .true.
    end do
    allocate (grd(3,mol%nat),source=0.0_wp)
    call engrad(mol,calc,etot,grd,io)
    if (io /= 0) then
      write (stdout,*)
      write (stdout,*) ' Initial singlepoint calculation failed!'
      write (stdout,*) ' Please check your input.'
      error stop
    end if
    do i = 1,calc%ncalculations
      if (allocated(calc%calcs(i)%wbo)) then
        wbo = calc%calcs(i)%wbo
        exit
      end if
    end do
    if (.not.allocated(wbo)) then
      write (stdout,*) 'No WBOs from this level, the topology falls back to CNs.'
      return
    end if
    open (newunit=ich,file='wbo')
    do i = 1,mol%nat
      do j = 1,i-1
        if (wbo(j,i) > 0.1_wp) write (ich,*) j,i,wbo(j,i)
      end do
    end do
    close (ich)
    write (stdout,*) 'WBO successfully calculated.'
  end subroutine msreact_wbo_calc

!========================================================================================!

  subroutine msreact_basic_atoms(env,basicatlist)
    !*******************************************************
    !* atoms of the input carrying a LP or pi LMO
    !* (targets of the attractive H-shift potential)
    !*******************************************************
    type(systemdata),intent(inout) :: env
    integer,allocatable,intent(out) :: basicatlist(:)
    type(coord) :: mol
    logical,allocatable :: basat(:)
    character(len=:),allocatable :: lvl
    integer :: i,io

    call mol%open('coord')
    mol%chrg = env%chrg
    mol%uhf = env%uhf
    select case (trim(env%gfnver))
    case ('--gfn1')
      lvl = 'gfn1'
    case ('--gfn0')
      lvl = 'gfn0'
    case default
      lvl = 'gfn2'
    end select
    write (stdout,*)
    write (stdout,'(a)',advance='no') ' LMO calculation ... '
    allocate (basat(mol%nat))
    call lmo_basic_atoms(mol,lvl,300.0_wp,basat,io)
    if (io /= 0) error stop 'LMO calculation for the H-shift sites failed'
    write (stdout,'(a)') 'done.'
    basicatlist = pack([(i,i=1,mol%nat)],basat)
    write (stdout,"('Protonation sites at atom positions:', *(i5))") basicatlist
  end subroutine msreact_basic_atoms

!========================================================================================!

  subroutine msreact_optjobs(env,mso,njob,jobs)
    !**********************************************************
    !* Optimize all distortions and append the results to the
    !* product list, in job order. Constrained jobs: loose
    !* optimization with the bond constraint, then at most
    !* maxc crude cycles without it (partial results are kept,
    !* the cycle limit keeps fragments from recombining).
    !* Other jobs: crude optimization.
    !**********************************************************
    type(systemdata),intent(inout) :: env
    type(msobj),intent(inout) :: mso
    integer,intent(in) :: njob
    type(msjob),intent(in) :: jobs(:)
    type(calcdata) :: base
    type(coord) :: ref
    real(wp),allocatable :: xyzout(:,:,:),eout(:)
    logical,allocatable :: ok(:)
    integer :: i,T,Tn,nat,chrg
    logical :: nested
    external :: ompmklset

    if (njob < 1) return
    ref%nat = mso%pl%mol(1)%nat
    ref%at = mso%pl%mol(1)%at
    chrg = mso%pl%mol(1)%chrg
    nat = ref%nat
    call msreact_basecalc(env,mso%T,env%chrg,env%uhf,base)
    allocate (xyzout(3,nat,njob),source=0.0_wp)
    allocate (eout(njob),source=0.0_wp)
    allocate (ok(njob),source=.false.)

    call new_ompautoset(env,'auto_nested',njob,T,Tn)
    nested = env%omp_allow_nested
    !$omp parallel do schedule(dynamic) default(shared) private(i)
    do i = 1,njob
      if (nested) call ompmklset(Tn)
      call msreact_runjob(base,ref,mso%maxc,jobs(i),xyzout(:,:,i),eout(i),ok(i))
    end do
    !$omp end parallel do

    do i = 1,njob
      if (.not.ok(i)) cycle
      call mso%pl%append(nat,ref%at,xyzout(:,:,i),eout(i),chrg,0,mso%pl%mol(1)%molmass)
    end do
    write (stdout,'(1x,i0,a,i0,a)') count(ok),' of ',njob,' optimizations successful'
  end subroutine msreact_optjobs

!========================================================================================!

  subroutine msreact_runjob(base,ref,maxc,job,xyz,etot,ok)
    !*******************************************************
    !* one distortion job, xyz in and out in Angstroem
    !*******************************************************
    type(calcdata),intent(in) :: base
    type(coord),intent(in) :: ref
    integer,intent(in) :: maxc
    type(msjob),intent(in) :: job
    real(wp),intent(out) :: xyz(3,ref%nat),etot
    logical,intent(out) :: ok
    type(calcdata) :: calc
    type(constraint) :: cons
    type(coord) :: mol,molnew
    real(wp),allocatable :: grd(:,:)
    integer :: io

    ok = .false.
    etot = 0.0_wp
    xyz = 0.0_wp
    mol%nat = ref%nat
    mol%at = ref%at
    mol%xyz = job%xyz*aatoau
    allocate (grd(3,mol%nat),source=0.0_wp)

    !$omp critical
    call calc%copy(base)
    !$omp end critical
    if (job%constr) then
!>--- xtb's fc*(r-r0)**2 is 1/2*k*(r-r0)**2 with k = 2*fc
      call cons%bondconstraint(job%A,job%B,job%D*aatoau,2.0_wp*job%fc)
      call calc%add(cons)
      calc%optlev = -1
      calc%maxcycle = 0
      call optimize_geometry(mol,molnew,calc,etot,grd,.false.,.false.,io)
      if (io /= 0) return
      mol%xyz = molnew%xyz
      !$omp critical
      call calc%copy(base)
      !$omp end critical
      calc%optlev = -3
      calc%maxcycle = maxc
      call optimize_geometry(mol,molnew,calc,etot,grd,.false.,.false.,io)
      if (io /= 0.and.io /= maxc) return
    else
      calc%optlev = -3
      calc%maxcycle = 0
      call optimize_geometry(mol,molnew,calc,etot,grd,.false.,.false.,io)
      if (io /= 0) return
    end if
    xyz = molnew%xyz*autoaa
    ok = .true.
  end subroutine msreact_runjob

!========================================================================================!
!========================================================================================!
end module msreact_calc
