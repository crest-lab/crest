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

!----------------------------------------------------------------------------------------
! CREST modus 3: MDOPT
!----------------------------------------------------------------------------------------
subroutine mdopt_legacy(env,tim)
  use crest_data
  use iomod
  use strucrd,only:rdensembleparam,rdensemble
  use utilities
  use cregen_interface
  implicit none

  type(systemdata) :: env
  type(timer)     :: tim
  character(len=256) :: tmpname
  character(len=128) :: newname

  integer :: iz2,nall

  associate (ensemblename => env%ensemblename)

!---- start timer
    call tim%start(1,'MDOPT')

!---- output file name
    newname = ensemblefile !'crest_ensemble.xyz'
    if (env%optpurge) newname = 'crest_optpurge.xyz'

!---- which file to optimize?
    tmpname = trim(ensemblename)
    if (index(tmpname,'none selected') .ne. 0) then
      write (*,*) 'No file selected with "--mdopt" option!'
      error stop
    end if
    call rdensembleparam(tmpname,iz2,nall)

!---- Small Header
    if (any((/-9224,-9225/) == env%properties)) then
      continue
    else
      write (*,*)
      if (.not.env%optpurge) then
        write (*,'(5x,''======================================='')')
        write (*,'(5x,''|               M D O P T             |'')')
        write (*,'(5x,''======================================='')')
        write (*,'(3x,''Optimization along trajectory (or ensemble).'')')
      else
        write (*,'(5x,''======================================='')')
        write (*,'(5x,''|          O P T - P U R G E          |'')')
        write (*,'(5x,''======================================='')')
        write (*,'(4x,''Optimization along ensemble to crosscheck'')')
        write (*,'(4x,''convergence into the same PES minima.'')')
      end if
      write (*,*)
    end if
    write (*,'(1x,a,a)') 'Input file: ','<'//trim(tmpname)//'>'
    write (*,'(1x,a,i0,a)') 'Containing ',nall,' structures.'
    if (env%optpurge) then
      call copy(trim(tmpname),trim(tmpname)//'.backup')
      write (*,'(1x,a)') 'Saved backup of input file.'
    end if

!---- Do the optimization
    call MDopt_para(env,trim(tmpname),0)

!---- printout and copy
    call rename('OPTIM'//'/'//'opt.xyz',trim(newname))
    write (*,'(1x,a,1x,a)') 'Optimized ensemble on file','<'//trim(newname)//'>'
!---- clean up
    call rmrf('OPTIM')

!--- purge identicals
    if (env%optpurge) then
      env%ensemblename = newname
      env%confgo = .true.
      call newcregen(env,9)
      call cleanpurge(tmpname,env%ethrpurge) !applied to original file!
    end if

!---- get the timing
    call tim%stop(1)
  end associate
contains
  subroutine cleanpurge(fname,ethr)
    use crest_parameters
    use strucrd
    implicit none
    character(len=*) :: fname
    integer :: nall,nat
    integer,allocatable :: groups(:)
    integer :: ngrps
    integer :: ich,i,j,ich2,jref
    real(wp),allocatable :: xyz(:,:,:)
    integer,allocatable :: at(:)
    real(wp),allocatable :: er(:)
    real(wp) :: elast,ediff,ethr
    character(len=*),parameter :: newfile = 'crest_ensemble.xyz'
    logical,allocatable :: track(:)
    logical :: paste

    call rdensembleparam(fname,nat,nall)
    allocate (xyz(3,nat,nall),at(nat),er(nall))
    call rdensemble(fname,nat,nall,at,xyz,er)

    open (newunit=ich,file='.groups') !written by newcregen(env,9)
    read (ich,*) nall,ngrps
    allocate (groups(nall))
    do i = 1,nall
      read (ich,*) j,groups(i)
    end do
    close (ich,status='delete')

    write (*,'(/,1x,a,1x,a)') 'Check duplicate info and energy windows (ETHR):'
    allocate (track(nall),source=.false.)
    !-- write new ensemble
    open (newunit=ich,file=newfile)
    open (newunit=ich2,file='.duplicates')
    do i = 1,ngrps
      paste = .true.
      elast = 0.0_wp
      do j = 1,nall
        if (groups(j) == i) then
          if (paste) then
            !call wrxyz(ich,nat,at,xyz(:,:,j),er(j))
            track(j) = .true.
            paste = .false.
            jref = j
            elast = er(j)
          else
            ediff = (er(j)-elast)*autokcal
            if (ediff .lt. ethr) then
              write (*,'(1x,a,i0,a)') 'Structure ',j,' was discarded'
              track(j) = .false.
              write (ich2,*) j,jref
            else
              !call wrxyz(ich,nat,at,xyz(:,:,j),er(j))
              track(j) = .true.
            end if
            elast = er(j)
          end if
        end if
      end do
    end do
    !-- write the ensemble
    do j = 1,nall
      if (track(j)) then
        call wrxyz(ich,nat,at,xyz(:,:,j),er(j))
      end if
    end do
    close (ich2)
    close (ich)
    write (*,'(/,1x,a,1x,a)') 'Unique structures written to file','<'//newfile//'>'
    write (*,'(1x,a,1x,a)') 'Duplicative structures in file','<'//trim(fname)//'>'
    write (*,'(1x,a)') ' have been documented in hidden file <.duplicates>'

    deallocate (er,at,xyz)
    return
  end subroutine cleanpurge
end subroutine mdopt_legacy

!----------------------------------------------------------------------------------------
! CREST modus 4: SCREEN
!----------------------------------------------------------------------------------------
subroutine screen_legacy(env,tim)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:rdensembleparam,rdensemble
  use utilities
  implicit none

  type(systemdata) :: env
  type(timer)      :: tim
  character(len=256) :: tmpname,newname

  integer :: iz2,nall,nremain

  character(len=128) :: inpnam,outnam
  character(len=512) :: thispath,filename

  real(wp) :: ewinbackup,rthrbackup

  call getcwd(thispath)

  ewinbackup = env%ewin
  rthrbackup = env%rthr

!---- start timer
  call tim%start(1,'SCREEN')

!---- output file name
  newname = ensemblefile !'crest_ensemble.xyz'

!---- which file to optimize?

  tmpname = trim(env%ensemblename)
  if (index(tmpname,'none selected') .ne. 0) then
    write (*,*) 'No file slected with "--screen" option!'
    error stop
  end if
  call rdensembleparam(tmpname,iz2,nall)

!---- Small Header
  write (*,*)
  write (*,'(5x,''======================================='')')
  write (*,'(5x,''|              S C R E E N            |'')')
  write (*,'(5x,''======================================='')')
  write (*,'(1x,''Multilevel optimization and structure screening.'')')
  write (*,*)
  write (*,'(1x,a,a)') 'Input file: ','<'//trim(tmpname)//'>'
  write (*,'(1x,a,i0,a)') 'Containing ',nall,' structures.'
  write (*,*)
!---- Multilevel job calls start here

  call smallhead('1. crude pre-optimization')
  call checkname_xyz(crefile,inpnam,outnam)
  call MDopt_para(env,trim(tmpname),1)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  env%ewin = ewinbackup*4.0d0
  env%rthr = rthrbackup*2.0d0
  !---using cregen to sort the optimized structures
  call checkname_xyz(crefile,inpnam,outnam)
  call confg_chk3(env)
  !--- only in the first ste  the input file is overwritten with the sorted file
  call remove(inpnam)
  call rename(outnam,inpnam)
  !-----
  call remaining_in(trim(inpnam),env%ewin,nremain)
  env%ewin = ewinbackup
  env%rthr = rthrbackup
  call rmrf('OPTIM')
  write (*,*)
!--- then vloose optimization
  call smallhead('2. optimization with loose thresholds')
  call checkname_xyz(crefile,inpnam,outnam)
  call MDopt_para(env,trim(inpnam),2)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  env%ewin = ewinbackup*2.0d0
  call confg_chk3(env)
  call remaining_in(trim(inpnam),env%ewin,nremain)
  env%ewin = ewinbackup
  env%rthr = rthrbackup
  write (*,*)
!--- lastly vtight optimization
  call smallhead('3. optimization with very tight thresholds')
  call checkname_xyz(crefile,inpnam,outnam)
  call MDopt_para(env,trim(inpnam),0)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call confg_chk3(env)
  call checkname_xyz(crefile,inpnam,outnam)
  call remaining_in(trim(inpnam),env%ewin,nremain)
  call rmrf('OPTIM')
  write (*,*)

!---- printout
  call catdel('cregen.out.tmp')
  write (*,*)
  write (*,'(1x,a,1x,a)') 'Final ensemble on file','<'//trim(newname)//'>'

  call rename(conformerfile,trim(newname))

!---- clean up
  call screen_cleanup

!---- get the timing
  call tim%stop(1)
end subroutine screen_legacy

