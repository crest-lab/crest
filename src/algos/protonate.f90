!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022 Philipp Pracht
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

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!> New implementation of protonation routines
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
subroutine crest_new_protonate(env,tim)
!********************************************************************
!* Standalone runtype for protonating a molecule
!*
!* Input/Output:
!*  env  -  crest's systemdata object
!*  tim  -  timer object
!********************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use parallel_interface
  use cregen_interface
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn,np
  logical :: pr,wr
  character(len=80) :: atmp
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp

  real(wp) :: energy,dip
  real(wp),allocatable :: grad(:,:)
  type(calcdata),allocatable :: tmpcalc
  type(calculation_settings) :: tmpset
  real(wp),allocatable :: protxyz(:,:)
  integer :: natp,pstep
  integer,allocatable :: atp(:)
  real(wp),allocatable :: xyzp(:,:,:)
  real(wp),allocatable :: ep(:)
  character(len=*),parameter :: partial = '∂E/∂'
!========================================================================================!
  write (stdout,*)
  !call system('figlet singlepoint')
  write (stdout,*) "                 _                    _        "
  write (stdout,*) " _ __  _ __ ___ | |_ ___  _ __   __ _| |_ ___  "
  write (stdout,*) "| '_ \| '__/ _ \| __/ _ \| '_ \ / _` | __/ _ \ "
  write (stdout,*) "| |_) | | | (_) | || (_) | | | | (_| | ||  __/ "
  write (stdout,*) "| .__/|_|  \___/ \__\___/|_| |_|\__,_|\__\___| "
  write (stdout,*) "|_|                                            "
  write (stdout,*) "-----------------------------------------------"
  write (stdout,*) "  automated protonation site screening script  "
  write (stdout,*) 'Cite as:'
  write (stdout,*) 'P.Pracht, C.A.Bauer, S.Grimme'
  write (stdout,*) 'JCC, 2017, 38, 2618–2631.'
  write (stdout,*)

!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()
!========================================================================================!
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!
!>--- The first step is to perfom a LMO calculation to identify suitable protonation sites
  call tim%start(14,'LMO center calculation')
  write (stdout,'(a)') repeat('-',80)
  write (stdout,'(a)')
  write (stdout,'(a)',advance='no') '> Setting up GFN0-xTB for LMO center calculation ... '
  flush (stdout)

  allocate (grad(3,mol%nat),source=0.0_wp)
!>--- GFN0 job adapted from global settings
  allocate (tmpcalc)
  call env2calc(env,tmpcalc,mol)
  tmpcalc%calcs(1)%id = jobtype%gfn0
  tmpcalc%calcs(1)%rdwbo = .true.
  tmpcalc%calcs(1)%getlmocent = .true.
  tmpcalc%calcs(1)%chrg = env%chrg
  tmpcalc%ncalculations = 1
  write (stdout,'(a)') 'done.'
!>--- and then start it
  write (stdout,'(a)',advance='yes') '> Performing singlepoint calculation ... '
  flush (stdout)
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
  call engrad(mol,tmpcalc,energy,grad,io)
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<!
  call tim%stop(14)
  write (stdout,'(a)') '> done.'
  write (atmp,'(a)') '> Total wall time for calculations'
  call tim%write_timing(stdout,14,trim(atmp),.true.)
  write (stdout,'(a)') repeat('-',80)
  write (stdout,*)
  if (io /= 0) then
    write (stdout,*)
    write (stdout,*) 'WARNING: Calculation exited with error!'
  end if
!>--- check LMO center
  np = tmpcalc%calcs(1)%nprot
  if (np > 0) then
    write (stdout,'(a,i0,a)') '> ',np,' π- or LP-centers identified as protonation candidates:'
    call move_alloc(tmpcalc%calcs(1)%protxyz,protxyz)
    write (stdout,'(1x,a5,1x,a,5x,a)') 'LMO','type','center(xyz/Ang)'
    do i = 1,np
      select case (nint(protxyz(4,i)))
      case (2)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'LP   ',protxyz(1:3,i)*autoaa
      case (3)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'π    ',protxyz(1:3,i)*autoaa
      case (4)
        write (stdout,'(1x,i5,1x,a,3F12.5)') i,'del.π',protxyz(1:3,i)*autoaa
      end select
    end do
  else
    write (stdout,*)
    write (stdout,*) 'WARNING: No suitable protonation sites found!'
    write (stdout,*) '         Check if you expect π- or LP-centers for your molecule!'
    write (stdout,*)
    return
  end if
  deallocate (tmpcalc)
  deallocate (grad)

!========================================================================================!
!>--- If we reached this point, we have candidate positions for our protonation!
  pstep = 0
  write (stdout,'(a)',advance='yes') '> Generating candidate structures ... '
  flush (stdout)
  natp = mol%nat+1
  allocate (atp(natp),source=0)
  allocate (xyzp(3,natp,np),ep(np),source=0.0_wp)

  call protonation_candidates(env,mol,natp,np,protxyz,atp,xyzp)
!>--- NOTE: after this the global charge (env%chrg) and all charges saved in the calc levels have been increased

  write (atmp,'(a,i0,a)') 'protonate_',pstep,'.xyz'
  write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with candidates ... '

  call wrensemble('protonate_0.xyz',natp,np,atp,xyzp*autoaa,ep)

  write (stdout,'(a)') '> done.'
  write (stdout,*)

!>--- Enforce further constraints, conditions, etc.
!     (TODO)

!>--- Optimize candidates, optional FF pre-step
  if (env%protb%ffopt) then
    call smallhead('Protomer Ensemble FF Pre-Optimization')
    allocate (tmpcalc)
    tmpcalc%optnewinit = .true.
    call tmpset%create('gfnff')
    tmpset%chrg = env%chrg
    call tmpcalc%add(tmpset)
    tmpcalc%maxcycle = 10000

    call tim%start(15,'Ensemble optimization (FF)')
    call print_opt_data(tmpcalc,stdout)
    write (stdout,'(a,i0,a)') '> ',np,' structures to optimize ...'
    call crest_oloop(env,natp,np,atp,xyzp,ep,.true.,tmpcalc)
    call tim%stop(15)

    deallocate (tmpcalc)

    pstep = pstep+1
    write (atmp,'(a,i0,a)') 'protonate_',pstep,'.xyz'
    write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
    call rename(ensemblefile,trim(atmp))

    !>--- sorting
    write (stdout,'(a)') '> Sorting structures by energy to remove failed opts. ...'
    call newcregen(env,6,trim(atmp))
    call rename(trim(atmp)//'.sorted',trim(atmp))
    write (stdout,*)

    !>--- re-read sorted ensemble
    deallocate(xyzp,atp) !> clear this space to re-use it
    call rdensemble(trim(atmp),natp,np,atp,xyzp)
  end if

!>--- H-position-only optimization
!    (TODO)

!>--- Optimize with global settings
  call smallhead('Protomer Ensemble Optimization')
  call tim%start(16,'Ensemble optimization')
  call print_opt_data(env%calc,stdout)
  write (stdout,'(a,i0,a)') '> ',np,' structures to optimize ...'
  call crest_oloop(env,natp,np,atp,xyzp,ep,.true.)
  call tim%stop(16)

  pstep = pstep+1
  write (atmp,'(a,i0,a)') 'protonate_',pstep,'.xyz'
  write (stdout,'(a,a,a)') '> Write ',trim(atmp),' with optimized structures ... '
  call rename(ensemblefile,trim(atmp))

!>--- sorting
  write (stdout,'(a)') '> Sorting structures by energy ...'
  call newcregen(env,6,trim(atmp))
  write (stdout,'(a)') '> sorted file was renamed to protonated.xyz'
  call rename(trim(atmp)//'.sorted','protonated.xyz')

!========================================================================================!
  return
end subroutine crest_new_protonate

!========================================================================================!
!========================================================================================!

subroutine protonation_candidates(env,mol,natp,np,protxyz,at,xyz)
!********************************************************
!* generate protonation/ionization candidate structures
!* The outputs are at and xyz, the latter being in Bohr
!********************************************************
  use crest_data
  use crest_parameters
  use strucrd,only:coord
  implicit none
  !> INPUT
  type(systemdata),intent(inout) :: env
  type(coord),intent(in) :: mol
  integer,intent(in) :: natp
  integer,intent(in) :: np
  real(wp),intent(in) :: protxyz(4,np)
  !> OUTPUT
  integer,intent(out)  :: at(natp)
  real(wp),intent(out) :: xyz(3,natp,np)
  !> LOCAL
  integer :: i,j,k,l
  integer :: ati,ichrg

  if (natp .ne. mol%nat+1) then
    write (stdout,'(a)') 'WARNING: Inconsistent number of atoms in protonation routine'
    error stop
  end if

  if (env%protb%swelem) then
!>--- User-defined monoatomic ion
    ichrg = env%protb%swchrg
    write (stdout,'(a,i0)') '> Increasing the molecular charge by ',ichrg
    call env%calc%increase_charge(ichrg)
    env%chrg = env%chrg+ichrg
    ati = env%protb%swat
  else
!>--- DEFAULT: H⁺
    write (stdout,'(a)') '> Increasing the molecular charge by 1'
    call env%calc%increase_charge()
    env%chrg = env%chrg+1
    ati = 1
  end if

!>--- Populate
  do i = 1,np
    do j = 1,mol%nat
      xyz(1:3,j,i) = mol%xyz(1:3,j)
      at(j) = mol%at(j)
    end do
    xyz(1:3,natp,i) = protxyz(1:3,i)
    at(natp) = ati
  end do
end subroutine protonation_candidates

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
