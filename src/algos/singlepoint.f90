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
!> Implementation of a singlepoint calculation
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!========================================================================================!
subroutine crest_singlepoint(env,tim)
!********************************************************************
!* Standalone runtype for a singlepoint calculation
!*
!* Input/Output:
!*  env  -  crest's systemdata object
!*  tim  -  timer object
!********************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use gradreader_module,only:write_engrad
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,T,Tn
  logical :: pr,wr
  character(len=80) :: atmp
!========================================================================================!
  type(calcdata) :: calc
  real(wp) :: accuracy,etemp

  real(wp) :: energy,dip
  real(wp),allocatable :: grad(:,:)

  character(len=*),parameter :: partial = '∂E/∂'
!========================================================================================!
  write (stdout,*)
  !call system('figlet singlepoint')
  write (stdout,*) "     _             _                  _       _   "
  write (stdout,*) " ___(_)_ __   __ _| | ___ _ __   ___ (_)_ __ | |_ "
  write (stdout,*) "/ __| | '_ \ / _` | |/ _ \ '_ \ / _ \| | '_ \| __|"
  write (stdout,*) "\__ \ | | | | (_| | |  __/ |_) | (_) | | | | | |_ "
  write (stdout,*) "|___/_|_| |_|\__, |_|\___| .__/ \___/|_|_| |_|\__|"
  write (stdout,*) "             |___/       |_|                      "
  write (stdout,*)
!========================================================================================!
  call new_ompautoset(env,'max',0,T,Tn)
  call ompprint_intern()
  call tim%start(14,'Singlepoint calculation')
!========================================================================================!
  call env%ref%to(mol)
  write (stdout,*)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)
!========================================================================================!

  write (stdout,'(a)') repeat('-',80)
  write (stdout,'(a)')

  allocate (grad(3,mol%nat),source=0.0_wp)
  calc = env%calc

!>--- print some info about the calculation
  call calc%info(stdout)

!>--- and then start it
  write (stdout,'(a)') repeat('-',80)
  write (stdout,'(a)',advance='no') '> Performing singlepoint calculations ... '
  flush (stdout)
  call engrad(mol,calc,energy,grad,io)
  call tim%stop(14)
  write (stdout,*) 'done.'
  write (atmp,'(a)') '> Total wall time for calculations'
  call tim%write_timing(stdout,14,trim(atmp),.true.)
  write (stdout,'(a)') repeat('-',80)
  if (io /= 0) then
    write (stdout,*)
    write (stdout,*) 'WARNING: Calculation exited with error!'
    return
  end if

!>--- print out the results
  if (any(calc%calcs(:)%rdwbo)) then
    write (stdout,*)
    write (stdout,*) 'Connectivity information (bond order):'
    do k = 1,calc%ncalculations
      if (calc%calcs(k)%rdwbo) then
        write (stdout,'("> ",a,i0)') 'Calculation level ',k
        write (stdout,'(a12,a12,a10)') 'Atom A','Atom B','BO(A-B)'
        do i = 1,mol%nat
          do j = 1,i-1
            if (calc%calcs(k)%wbo(i,j) > 0.01_wp) then
              write (stdout,*) i,j,calc%calcs(k)%wbo(i,j)
            end if
          end do
        end do
      end if
    end do
    write (stdout,*)
    write (stdout,'(a)') repeat('-',80)
  else
    write (stdout,*)
  end if

  if (any(calc%calcs(:)%rddip)) then
    write (stdout,*)
    write (stdout,*) 'Molecular dipole moments (a.u.):'
    do k = 1,calc%ncalculations
      if (calc%calcs(k)%rddip) then
        dip = norm2(calc%calcs(k)%dipole)
        write (stdout,'("> ",a,i0)') 'Calculation level ',k
        write (stdout,'(a10,a10,a10,a12)') 'x','y','z','tot (Debye)'
        write (stdout,'(4f10.3)') calc%calcs(k)%dipole(:),dip*autod
      end if
    end do
    write (stdout,*)
    write (stdout,'(a)') repeat('-',80)
  else
    write (stdout,*)
  end if

  if (all(calc%calcs(:)%rdgrad.eqv..false.)) then
    write (stdout,'(a)') '> No gradients calculated'
  else
    write (stdout,'(a)') '> Final molecular gradient ( Eh/a0 ):'
    write (stdout,'(13x,a,13x,a,13x,a)') partial//'x',partial//'y',partial//'z'
    do i = 1,mol%nat
      write (stdout,'(3f18.8)') grad(1:3,i)
    end do
    write (stdout,'(a,f18.8,a)') '> Gradient norm:',norm2(grad),' Eh/a0'
  end if

  if (calc%ncalculations > 1) then
    write (stdout,*)
    write (stdout,'(a)') '> Individual energies and gradient norms:'
    do k = 1,calc%ncalculations
      write (stdout,'(1x,a,i3,2f18.8)') 'calculation ',k,calc%etmp(k),norm2(calc%grdtmp(:,:,k))
    end do
    if (calc%nconstraints > 0) then
      write (stdout,'(1x,a)') '(+ constraints contribution)'
    end if
  end if

  write (stdout,*)
  write (stdout,'(a)') repeat('=',40)
  write (stdout,'(1x,a,f20.10,a)') 'TOTAL ENERGY ',energy,' Eh'
  write (stdout,'(1x,a,f20.10,a)') 'GRADIENT NORM',norm2(grad),' Eh/a0'
  write (stdout,'(a)') repeat('=',40)

  write (stdout,'(1x,a)') 'Writing crest.engrad ...'
  call write_engrad('crest.engrad',energy,grad)

  if (env%testnumgrad) then
    call numgrad(mol,calc,grad)
  end if

  deallocate (grad)
!========================================================================================!
  return
end subroutine crest_singlepoint

!========================================================================================!
!========================================================================================!

subroutine crest_xtbsp(env,xtblevel,molin)
!********************************************************************
!* Replacement for the legacy xtbsp routine, makes use of gfn0 or tblite.
!* The purpose of this routine is usually to generate WBOs
!*
!* Input/Output:
!*  env      - crest's systemdata object
!*  xtblevel - quick selection of calc. level
!*  molin    - molecule data
!********************************************************************
  use crest_parameters
  use crest_data
  use crest_calculator
  use strucrd
  use wiberg_mayer,only:write_wbo
  implicit none
  !> INPUT
  type(systemdata) :: env
  integer,intent(in),optional :: xtblevel
  type(coord),intent(in),optional :: molin
  !> LOCAL
  integer :: lv,io
  type(calcdata) :: tmpcalc
  type(calculation_settings) :: cal
  type(coord) :: mol
  real(wp) :: energy
  real(wp),allocatable :: grad(:,:)

  !>--- transfer settings from env to tmpcalc
  call env2calc(env,tmpcalc,molin)

  tmpcalc%calcs(1)%rdwbo = .true. !> obtain WBOs
  if (present(xtblevel)) then
    select case (xtblevel) !> no default
    case (0)
      tmpcalc%calcs(1)%id = jobtype%gfn0
    case (1)
      tmpcalc%calcs(1)%id = jobtype%tblite
      tmpcalc%calcs(1)%tblitelvl = 1
    case (2)
      tmpcalc%calcs(1)%id = jobtype%tblite
      tmpcalc%calcs(1)%tblitelvl = 2
    end select
  end if
  if (present(molin)) then
    mol = molin
  else
    call mol%open('coord')
  end if
  allocate (grad(3,mol%nat),source=0.0_wp)

  call engrad(mol,tmpcalc,energy,grad,io)
  if (io .ne. 0) then
    error stop 'crest_xtbsp failed'
  end if

  !>--- write wbo file
  if (tmpcalc%calcs(1)%rdwbo) then
    call write_wbo(tmpcalc%calcs(1)%wbo,0.002_wp)

  end if

  deallocate (grad)
  call tmpcalc%reset()
  call mol%deallocate()
end subroutine crest_xtbsp

!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
subroutine crest_ensemble_singlepoints(env,tim)
!***********************************************
!* subroutine crest_ensemble_singlepoints
!* This routine implements a standalone runtype
!* to perform singlepoint evaluations along an
!* ensemble or trajectory file.
!***********************************************
  use crest_parameters,only:wp,stdout,bohr,angstrom
  use crest_data
  use crest_calculator
  use strucrd
  use optimize_module
  use utilities,only:dumpenergies
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,c
  logical :: pr,wr,ex
!========================================================================================!
  type(calcdata) :: calc

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)

  character(len=:),allocatable :: ensnam
  integer :: nat,nall,T,Tn
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  character(len=80) :: atmp
  real(wp) :: percent
  character(len=52) :: bar
!========================================================================================!
  write (*,*)
!>--- check for the ensemble file
  inquire (file=env%ensemblename,exist=ex)
  if (ex) then
    ensnam = env%ensemblename
  else
    write (stdout,*) 'no ensemble file provided.'
    return
  end if

!>--- start the timer
  call tim%start(14,'Ensemble singlepoints')

!>---- read the input ensemble
  call rdensembleparam(ensnam,nat,nall)
  if (nall .lt. 1) return
  allocate (xyz(3,nat,nall),at(nat),eread(nall))
  call rdensemble(ensnam,nat,nall,at,xyz,eread)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: crest_oloop requires coordinates in Bohrs
  xyz = xyz/bohr
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

!>--- set OMP parallelization
  call new_ompautoset(env,'auto',nall,T,Tn)

!========================================================================================!
  !>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,14x,"│")') "ENSEMBLE SINGLEPOINTS"
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)
  write (stdout,'(1x,a,i0,a,1x,a)') 'Evaluationg all ',nall,' structures of file',trim(ensnam)
  !>--- call the loop
  call crest_sploop(env,nat,nall,at,xyz,eread,.true.)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- Important: ensemble file must be written in AA
  xyz = xyz/angstrom
!>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!>--- write output ensemble
  call wrensemble(ensemblefile,nat,nall,at,xyz,eread)
  write (stdout,'(/,a,a,a)') 'Ensemble with updated energies written to <',ensemblefile,'>'

  call dumpenergies('crest.energies',eread)
  write (stdout,'(/,a,a,a)') 'List of energies written to <','crest.energies','>'

  deallocate (eread,at,xyz)
!========================================================================================!
  call tim%stop(14)
  return
end subroutine crest_ensemble_singlepoints

