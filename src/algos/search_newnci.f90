!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
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

subroutine crest_search_newnci(env,tim)
!*******************************************************************
!* This is a mofidied implementation of CREST's iMTD workflow
!* aimed at non-covalent systems.
!*******************************************************************
  use crest_parameters,only:wp,stdout
  use crest_data
  use crest_calculator
  use strucrd
  use dynamics_module
  use shake_module
  use iomod
  use utilities
  use cregen_interface
  implicit none
  type(systemdata),intent(inout) :: env
  type(timer),intent(inout)      :: tim
  type(coord) :: mol,molnew
  integer :: i,j,k,l,io,ich,m
  logical :: pr,wr
!===========================================================!
  type(calcdata) :: calc
  type(mddata) :: mddat
  type(shakedata) :: shk

  type(mddata),allocatable :: mddats(:)
  integer :: nsim,nallout

  real(wp) :: energy,gnorm
  real(wp),allocatable :: grad(:,:)
  character(len=:),allocatable :: ensnam
  integer :: nat,nall
  real(wp),allocatable :: eread(:)
  real(wp),allocatable :: xyz(:,:,:)
  integer,allocatable  :: at(:)
  logical :: dump,ex
  character(len=80) :: atmp,btmp,str
  logical :: multilevel(6)
  logical :: start,lower
!===========================================================!
!>--- printout header
  write (stdout,*)
  write (stdout,'(10x,"┍",49("━"),"┑")')
  write (stdout,'(10x,"│",14x,a,13x,"│")') " CREST NCI SAMPLING   "
  write (stdout,'(10x,"┕",49("━"),"┙")')
  write (stdout,*)

!===========================================================!
!>--- setup
  call env%ref%to(mol)
  write (stdout,*) 'Input structure:'
  call mol%append(stdout)
  write (stdout,*)

!>--- sets the MD length according to a flexibility measure
  call md_length_setup(env)
!>--- create the MD calculator saved to env
  call env_to_mddat(env)

  if (env%performMTD) then
!>--- (optional) calculate a short 1ps test MTD to check settings
    call tim%start(1,'Trial metadynamics (MTD)')
    call trialmd(env)
    call tim%stop(1)
  end if

!===========================================================!
!>--- Start mainloop
  env%nreset = 0
  start = .true.
  MAINLOOP: do
    call printiter
    if (.not.start) then
!>--- clean Dir for new iterations, but leave iteration backup files
      call clean_V2i
      env%nreset = env%nreset+1
    else
!>--- at the beginning, wipe directory clean
      call V2cleanup(.false.)
    end if
!===========================================================!
!>--- Meta-dynamics loop
    mtdloop: do i = 1,env%Maxrestart

      write (stdout,*)
      write (stdout,'(1x,a)') '------------------------------'
      write (stdout,'(1x,a,i0)') 'Meta-Dynamics Iteration ',i
      write (stdout,'(1x,a)') '------------------------------'

      nsim = -1 !>--- enambles automatic MTD setup in init routines
      call crest_search_multimd_init(env,mol,mddat,nsim)
      allocate (mddats(nsim),source=mddat)
      call crest_search_multimd_init2(env,mddats,nsim)

      call tim%start(2,'Metadynamics (MTD)')
      call crest_search_multimd(env,mol,mddats,nsim)
      call tim%stop(2)
!>--- a file called crest_dynamics.trj should have been written
      ensnam = 'crest_dynamics.trj'
!>--- deallocate for next iteration
      if (allocated(mddats)) deallocate (mddats)

!==========================================================!
!>--- Reoptimization of trajectories
      call tim%start(3,'Geometry optimization')
      call optlev_to_multilev(env%optlev,multilevel)
      call crest_multilevel_oloop(env,ensnam,multilevel)
      call tim%stop(3)

!>--- save the CRE under a backup name
      call checkname_xyz(crefile,atmp,str)
      call checkname_xyz('.cre',str,btmp)
      call rename(atmp,btmp)
!>--- save cregen output
      call checkname_tmp('cregen',atmp,btmp)
      call rename('cregen.out.tmp',btmp)

!=========================================================!
!>--- cleanup after first iteration and prepare next
      if (i .eq. 1.and.start) then
        start = .false.
!>-- obtain a first lowest energy as reference
        env%eprivious = env%elowest
!>-- remove the two extreme-value MTDs
        if (.not.env%readbias.and.env%runver .ne. 33.and. &
        &   env%runver .ne. 787878) then
          env%nmetadyn = env%nmetadyn-2
        end if
!>-- the cleanup
        call clean_V2i
!>-- and always do two cycles of MTDs
        cycle mtdloop
      end if
!=========================================================!
!>--- Check for lowest energy
      call elowcheck(lower,env)
      if (.not.lower) then
        exit mtdloop
      end if
    end do mtdloop
!=========================================================!
!>--- collect all ensembles from mtdloop and merge
    write (stdout,*)
    write (stdout,'(''========================================'')')
    write (stdout,'(''           MTD Simulations done         '')')
    write (stdout,'(''========================================'')')
    write (stdout,'(1x,''Collecting ensmbles.'')')
!>-- collecting all ensembles saved as ".cre_*.xyz"
    call collectcre(env)
    call newcregen(env,0)
    call checkname_xyz(crefile,atmp,btmp)
!>--- remaining number of structures
    call remaining_in(atmp,env%ewin,nallout)

!==========================================================!
!>--- exit mainloop
    exit MAINLOOP
  end do MAINLOOP

!==========================================================!
!>--- final ensemble optimization
  write (stdout,'(/)')
  write (stdout,'(3x,''================================================'')')
  write (stdout,'(3x,''|           Final Geometry Optimization        |'')')
  write (stdout,'(3x,''================================================'')')
  call tim%start(3,'Geometry optimization')
  call checkname_xyz(crefile,atmp,str)
  call crest_multilevel_wrap(env,trim(atmp),0)
  call tim%stop(3)

!==========================================================!
!>--- print CREGEN results and clean up Directory a bit
  write (stdout,'(/)')
  call smallhead('Final Ensemble Information')
  call V2terminating()

!========================================================================================!
!========================================================================================!
end subroutine crest_search_imtdgc
!========================================================================================!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!========================================================================================!
