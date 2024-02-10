!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Philipp Pracht
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

!=========================================================================================!
! iMTD(RMSD)-GC Algorithm (also Entropy mode and iMTD-sMTD Algo --v4)
! This is the algo for the conformational search.
!=========================================================================================!
subroutine confscript2i_legacy(env,tim)
  use iso_fortran_env,only:wp => real64
  use crest_data
  use iomod
  use strucrd,only:coord2xyz,xyz2coord
  use utilities
  use cregen_interface
  implicit none
  type(systemdata) :: env
  type(timer)   :: tim
  integer :: i,k,m
  integer :: eit,eit2
  real(wp) :: ewin
  integer :: dum,bref
  character(len=256) :: str,atmp,btmp
  integer :: nallout
  logical :: lower
  logical :: start
  logical :: ex
  logical :: stopiter,fail

  real(wp) :: autokcal
  parameter(autokcal=627.509541d0)

  settingLogs:associate (performMTD => env%performMTD,restartopt => env%restartopt, &
  &         rotamermds => env%rotamermds,performCross => env%performCross,           &
  &         doNMR => env%doNMR)

  settingData:associate (nat => env%nat,cgf => env%cgf,thresholds => env%thresholds)

!--- some defaults
  ewin = env%ewin    !EWIN

!--- prepare the directory
  if (performMTD) then
    call V2cleanup(restartopt)
  end if
!--- if we only have 2 atoms, do not do anything, except copying the coords
  if (nat .le. 2) then
    call catchdiatomic(env)
    return
  end if

!---- Get do a single trial MTD to test the settings
  call md_length_setup(env)  !> sets the MD length according to a flexibility measure
!---- set number of regular MDs
  call adjustnormmd(env)

  if (performMTD) then
    call tim%start(1,'test MD')
    call trialmd(env)    !calculate a short 1ps test MTD to check settings
    call tim%stop(1)
  end if

!---- copy the original coord
  call copy('coord','coord.original')
  call coord2xyz('coord','.history.xyz')

!===========================================================================================!
!===========================================================================================!

  if (env%performMTD) then
    write (*,*)
    write (*,*) 'list of Vbias parameters applied:'
    do m = 1,env%nmetadyn
      write (*,'(''$metadyn '',f10.5,f8.3,i5)') env%metadfac(m) / env%rednat,env%metadexp(m)
    end do
  end if

  env%nreset = 0
  start = .true.
  MAINLOOP: do
    if (env%iterativeV2) then
      call printiter
    end if
    if (.not. start) then
      call clean_V2i  !--- clean Dir for new iterations
      env%nreset = env%nreset + 1
    else !--at the beginning clean existing backup-ensembles
      call rmrfw('.cre_')
    end if
!-------- iterative loop over MTDs
    mtdloop: do i = 1,env%Maxrestart
      !---- Small Header
      write (*,*)
      write (*,'(''========================================'')')
      if (env%Maxrestart .gt. 1) then
        write (*,'(''            MTD Iteration '',i2,13x)') i
      else
        write (*,'(''            MTD Simulations '',15x)')
      end if
      write (*,'(''========================================'')')
      !---- do the MTDs
      call tim%start(2,'MTD')
      call MetaMD_para_OMP(env)
      call tim%stop(2)

      MULTILEVELSKIP: if (env%multilevelopt) then
        !---- Optimize using confopt (parallel)
        call tim%start(3,'multilevel OPT')
        write (*,*)
        write (*,'(''-----------------------'')')
        write (*,'(''Multilevel Optimization'')')
        write (*,'(''-----------------------'')')
        write (*,*)
        if (env%optlev >= -2.0d0) then
          call multilevel_opt(env,4)
        end if
        call append_INPUT_to('coord','input')   !include the input structure into the last optimization
        if (.not. env%superquick .and. env%optlev >= 1.0d0) then
          if (.not. env%entropic) then
            call multilevel_opt(env,5)
          else
            call multilevel_opt(env,99)   !--- the last CREGEN is done within this subroutine
          end if
        end if
        call tim%stop(3)
        !---- save the CRE under a backup name
        call checkname_xyz(crefile,atmp,str)
        call checkname_xyz('.cre',str,btmp)
        call rename(atmp,btmp)
        !---- save cregen output
        call checkname_tmp('cregen',atmp,btmp)
        call rename('cregen.out.tmp',btmp)
      else
        exit mtdloop
      end if MULTILEVELSKIP

      !---- in the first cycle just save the energy and cycle
      if (i .eq. 1 .and. start) then
        !---- clean the dir
        start = .false.       !--- only in the first cycle of MAINLOOP the MTDs are done at least 2 times
        if (.not. env%readbias .and. &
        &   env%runver .ne. 33 .and. &
        &   env%runver .ne. 787878 ) then
          dum = env%nmetadyn - 2
          env%nmetadyn = dum !--- only in the first cycle two MTDs are done additionally with extreme Vbias settings
        end if
        call clean_V2i
        env%eprivious = env%elowest !in the first cycle only save the energy
        !---- save the new best conformer
        inquire (file='crest_best.xyz',exist=ex)
        if (ex) then
          call XYZappendto('crest_best.xyz','.history.xyz')
          call xyz2coord('crest_best.xyz','coord')             !new reference coord to start the MTDs with
        end if

        cycle  !--- always do at least 2 cycles
      end if

      !---- check elowest
      call elowcheck(lower,env)
      if (.not. lower) then
        exit mtdloop
      end if
    end do mtdloop

!--------
    write (*,'(''========================================'')')
    if (env%Maxrestart .gt. 1) then
      write (*,'(''            MTD Iterations done         '')')
    else
      write (*,'(''           MTD Simulations done         '')')
    end if
    write (*,'(''========================================'')')
    write (*,'(1x,''Collecting ensembles.'')')
    call collectcre(env)                      !--- collecting all ensembles saved as ".cre_*.xyz"
      if (.not. env%entropic .and. env%crestver .ne. 22) then
        call newcregen(env,0)
      else
        call newcregen(env,2)
      end if
    call checkname_xyz(crefile,atmp,btmp)
    call remaining_in(atmp,ewin,nallout) !--- remaining number of structures
    write (*,*)
!====================================================================!
!---- (Optional) sampling of additional XH positions
    if (env%doOHflip) then
      call XHorient(env,conformerfile)
      inquire (file='oh_ensemble.xyz',exist=ex)
      if (ex) then
        call checkname_xyz(crefile,atmp,btmp)
        call appendto('oh_ensemble.xyz',atmp)
        call remove('oh_ensemble.xyz')
        if (.not. env%entropic .and. env%crestver .ne. 22) then
          call newcregen(env,0)
        else
          call newcregen(env,2)
        end if
        call remaining_in(btmp,ewin,nallout)
        write (*,*)
      end if
    end if
!====================================================================!
!---- Perform additional MDs on the lowest conformers
    if (env%rotamermds) then
      call tim%start(4,'MD ')
      call normalMD_para_OMP(env,env%nrotammds,env%temps)
      if (env%multilevelopt) then
        if (env%optlev >= 1.0d0) then
          call multilevel_opt(env,6)
        else
          call multilevel_opt(env,99)
        end if
        call elowcheck(lower,env)
      else
        lower = .false.
      end if
      call tim%stop(4)
      if (lower) then
        call checkname_xyz(crefile,atmp,str)
        call checkname_xyz('.cre',str,btmp)
        call rename(atmp,btmp)
        if (env%iterativeV2) cycle MAINLOOP
      end if
    end if
!====================================================================!
!---- Genetic crossing
    if (env%performCross) then
      call tim%start(5,'GC')
      call cross3(env)
      call tim%stop(5)
      call confg_chk3(env)
      call elowcheck(lower,env)
      if (lower) then
        call checkname_xyz(crefile,atmp,str)
        call checkname_xyz('.cre',str,btmp)
        call rename(atmp,btmp)
        if (env%iterativeV2) cycle MAINLOOP
      end if
    end if
!====================================================================!
!---- Entropy mode iterative statically biased MDs
    if (env%entropymd) then
      call mtdatoms(env)
      call tim%start(6,'static MTD')
      call emtdcopy(env,0,stopiter,fail)
      bref = env%emtd%nbias

      ENTROPYITER: do eit = 1,env%emtd%iter
        dum = nint(float(env%emtd%nbias) * env%emtd%nbiasgrow)
        !env%emtd%nbias = nint(float(env%emtd%nbias) * env%emtd%nbiasgrow)
        env%emtd%nbias = max(env%emtd%nbias + 1,dum)
        fail = .false.
        EFALLBACK: do k = 1,env%emtd%maxfallback
          call printiter2(eit)
          call tim%start(6,'static MTD')
          call entropyMD_para_OMP(env)
          call tim%stop(6)
          call emtdcheckempty(env,fail,env%emtd%nbias)
          if (fail) then
            if (k == env%emtd%maxfallback) then
              stopiter = .true.
            else
              cycle EFALLBACK
            end if
          else
            call tim%start(3,'multilevel OPT')
            if (env%optlev >= -1.0d0) then
              call multilevel_opt(env,2)
            end if
            call multilevel_opt(env,99)

            call tim%stop(3)
            !--- if in the entropy mode a lower structure was found
            !    --> cycle, required for extrapolation
            call elowcheck(lower,env)
            if (lower .and. env%entropic) then
              env%emtd%nbias = bref  !IMPORTANT, reset for restart
              cycle MAINLOOP
            end if
            !--- file handling
            eit2 = eit
            call emtdcopy(env,eit2,stopiter,fail)
            env%emtd%iterlast = eit2
          end if
          if (.not. lower .and. fail .and. .not. stopiter) then
            cycle EFALLBACK
          end if
          exit EFALLBACK  !fallback loop is exited on first opportuinity
        end do EFALLBACK
        if (stopiter) then
          exit ENTROPYITER
        end if
      end do ENTROPYITER
    end if

!-------
    exit MAINLOOP !--- if this point is reached, i.e., there weren't any further restarts, exit the loop
  end do MAINLOOP

!==========================================================================================!
!==========================================================================================!

  !if(.not.env%entropic .and. .not.(env%crestver == 22))then
  if (.not. env%entropymd) then
    !------ last optimization (with user set optlevel)
    write (*,*)
    write (*,*)
    write (*,'(3x,''================================================'')')
    write (*,'(3x,''|           Final Geometry Optimization        |'')')
    write (*,'(3x,''================================================'')')

    if (doNMR) cgf(3) = .true.  !--- if NMR equivalencies are requested, turn them on here
    call tim%start(7,'')
    call multilevel_opt(env,99) !--- the last CREGEN is done within this subroutine
    call tim%stop(7)                                 !--- optlevel is userset
  else
    !------ or just sort the last ensemble for entropy mode
    !if(doNMR) cgf(3)=.true.
    !call confg_chk3(env)
  end if

!>---- print CREGEN results and clean up Directory a bit
  if (env%crestver .ne. crest_solv) then
    call V2terminating()
  end if

  end associate settingData
  end associate settingLogs

end subroutine confscript2i_legacy

