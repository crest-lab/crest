!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht
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

!================================================================================!
! AUTOMATED PROTONATION AND ENERGETIC RANKING SCRIPT
! To use run:
!    crest <input> --protonate
!================================================================================!

subroutine prothead
  implicit none
  write (*,*) '       __________________________________________'
  write (*,*) '      |                                          |'
  write (*,*) '      |       automated protonation script       |'
  write (*,*) '      |__________________________________________|'
  write (*,*) ' Universitaet Bonn, MCTC'
  write (*,*) ' P.Pracht, Wed 28. Nov 13:11:52 CEST 2018'
  write (*,*)
  write (*,*) ' Cite as:'
  write (*,*) ' P.Pracht, C.A.Bauer, S.Grimme'
  write (*,*) ' JCC, 2017, 38, 2618–2631.'
  write (*,*)
end subroutine prothead

!--------------------------------------------------------------------------------------------
! Protonation workflow with GFNn-xTB
!--------------------------------------------------------------------------------------------
subroutine protonate(env,tim)
  use crest_parameters
  use crest_data
  use iomod
  use strucrd,only:coord2xyz
  use utilities
  implicit none
  type(systemdata) :: env
  type(timer)      :: tim
  type(protobj)    :: prot

  character(len=32)  :: dirn
  character(len=64)  :: protname
  character(len=256) :: thispath
  character(len=256) :: filename
  character(len=128) :: inpnam,outnam

  integer :: ich
  integer :: natp,nallout,refchrg

  logical :: ex

  interface
    subroutine xtblmo(env, print)
      import :: systemdata
      type(systemdata) :: env
      logical, optional :: print
    end subroutine xtblmo
  end interface

!--- printout & clean directory
  call protclean
  call prothead

  if (.not.allocated(env%ptb%atmap)) allocate (env%ptb%atmap(env%nat))
  if (.not.env%ptb%strictPDT.and..not.env%ptb%fixPDT) then
!--- sort the input file (H atoms to the bottom)
    call htothebottom('coord',env%chrg,env%nat,env%ptb%atmap)
  else
!--- or sort AND apply heavy atom bond constraints
    call PDT_constraints(env)
  end if

!--- get some settings
  call getcwd(thispath)
  dirn = 'PROT'
  protname = 'protonate_0.xyz'
  prot = env%ptb
  refchrg = env%chrg
  prot%newchrg = env%chrg+1  !increase chrg by one
  natp = env%nat+1 !additional proton, Nat is increased by one

!--- do the xTB calculation for the LMOs
  call tim%start(1,'LMO calc.')
  call xtblmo(env,.true.)
  call tim%stop(1)
  inquire (file='coordprot.0',exist=ex)
  if (.not.ex) then
    write (*,*)
    write (*,*) '***Warning***'
    write (*,*) 'No "coordprot.0" file was written, it is possible that'
    write (*,*) 'there are no suitable LP- or π-centers in the molecule.'
    write (*,*) 'Hence the procedure could not be automatized. (sorry)'
    write (*,*) '***Warning***'
    return
  end if

!--- get the new charge and set up the calculations
  call tim%start(2,'multilevel OPT')
  open (newunit=ich,file='.CHRG')
  write (ich,*) prot%newchrg     !new charge written here
  close (ich)
  write (*,*)
  write (*,'(''-----------------------'')')
  write (*,'(''Multilevel Optimization'')')
  write (*,'(''-----------------------'')')

  call coord2xyz('coordprot.0',trim(protname))
  call appendto('xtbscreen.xyz',protname)
  env%nat = natp

  !write(*,*) 'switching stuff:'
  !write(*,*) prot%swelem
  !write(*,*) prot%swat
  !write(*,*) prot%swchrg
  if (prot%swelem) then
    call swelem(protname,env)
  end if
  env%chrg = prot%newchrg !!all optimizations access env%chrg!!!

  call smallhead('1. crude pre-optimization')
  call checkname_xyz('protonate',inpnam,outnam)
  call MDopt_para(env,protname,1)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  !write(*,*) trim(filename)
  !call copy(trim(filename),'ensemble-test.xyz')
  if (prot%ABcorrection) then
    call prot_correction(env,trim(filename))
  end if
  prot%ewin = prot%ewin*3.0d0
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures
  !call prot_correction(env,trim(outnam))
  write (*,*)

  call smallhead('2. loose optimization')
  call checkname_xyz('protonate',inpnam,outnam)
  call MDopt_para(env,inpnam,2)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  if (prot%ABcorrection) then
    call prot_correction(env,trim(filename))
  end if
  prot%ewin = prot%ewin*(2.0d0/3.0d0)
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures
  write (*,*)

  call smallhead('3. optimization with user-defined thresholds')
  call checkname_xyz('protonate',inpnam,outnam)
  call MDopt_para(env,inpnam,0)
  filename = trim(thispath)//'/'//trim(outnam)
  call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
  call rmrf('OPTIM')
  if (prot%ABcorrection) then
    call prot_correction(env,trim(filename))
  end if
  prot%ewin = prot%ewin/2.0d0
  call sort_ens(prot,outnam,.false.)
  call remaining_in(outnam,prot%ewin,nallout) !--- remaining number of structures

  !call rename(outnam,'protonated.xyz')
  call cosort(outnam,'protonated.xyz',.false.,.false.)
  call sort_ens(prot,'protonated.xyz',.true.)
  call tim%stop(2)

!>--- (optional) post-processing
  if (env%relax) then
    env%rednat = env%rednat+1
    call relaxensemble('protonated.xyz',env,tim)
  end if

  if (env%outputsdf) then
    call new_wrsdfens(env,'protonated.xyz','protonated.sdf',.true.)
  end if

!--- reset data for main dir
  env%chrg = refchrg
  if (env%chrg .eq. 0) then
    call remove('.CHRG')
  else
    open (newunit=ich,file='.CHRG')
    write (ich,*) env%chrg
    close (ich)
  end if
  env%nat = natp-1 !reset nat
end subroutine protonate

!--------------------------------------------------------------------------------------------
! A quick single point xtb calculation and calculate LMOs
!--------------------------------------------------------------------------------------------
subroutine xtblmo(env,print)
  use crest_parameters
  use iomod
  use crest_data
  implicit none
  type(systemdata) :: env
  character(len=80) :: fname
  character(len=:),allocatable :: jobcall
  integer :: io
  character(len=*),parameter :: pipe = ' > xtb.out 2>/dev/null'
  logical, optional :: print ! leave the xtb.out file (e.g. for msreact mode)

!---- setting threads
  if (env%autothreads) then
    call ompautoset(env%threads,7,env%omp,env%MAXRUN,1) !set the global OMP/MKL variables for the xtb jobs
  end if
!---- new plain coord file
  fname = 'tmpcoord'
  call copy('coord',fname)
!---- jobcall
  write (*,*)
  write (*,'('' LMO calculation ... '')',advance='no')
  jobcall = trim(env%ProgName)
  jobcall = trim(jobcall)//' '//trim(fname)
  jobcall = trim(jobcall)//' '//trim(env%gfnver)
  jobcall = trim(jobcall)//' --sp --lmo '//trim(env%solv) 
  jobcall = trim(jobcall)//trim(pipe)
  call command(trim(jobcall),io)
  write (*,'(''done.'')')

!---- cleanup
  call remove(fname)
  if (.not. print) call remove('xtb.out')
  call remove('energy')
  call remove('charges')
  call remove('xtbrestart')
end subroutine xtblmo

!--------------------------------------------------------------------------------------------
! swithc the added proton to a nother element
!--------------------------------------------------------------------------------------------
subroutine swelem(iname,env)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd,only:rdensembleparam,rdensemble,wrxyz
  implicit none
  type(systemdata) :: env
  type(protobj) :: prot
  character(len=*) :: iname

  integer :: i,ich
  integer :: nat,nall
  integer :: nchrg
  real(wp),allocatable :: xyz(:,:,:)
  real(wp),allocatable :: eread(:)
  integer,allocatable  :: at(:)

  prot = env%ptb
  nchrg = env%chrg+prot%swchrg
  prot%newchrg = nchrg

  call rdensembleparam(iname,nat,nall)
  allocate (xyz(3,nat,nall),eread(nall),at(nat))
  call rdensemble(iname,nat,nall,at,xyz,eread)

  !---- write updated .CHRG file
  open (newunit=ich,file='.CHRG')
  write (ich,'(i6)') nchrg    !new charge written here
  close (ich)
  call remove(iname)

  open (newunit=ich,file=iname)
  at(nat) = prot%swat
  do i = 1,nall
    call wrxyz(ich,nat,at,xyz(:,:,i))
  end do
  close (ich)
  deallocate (at,eread,xyz)

  env%ptb = prot
  return
end subroutine swelem

subroutine swparse(iname,prot)
  use crest_parameters
  use iomod
  use crest_data
  use strucrd,only:i2e,e2i
  implicit none
  type(protobj) :: prot
  character(len=*) :: iname

  integer :: i,slen
  character(len=1) :: sig
  character(len=10) :: el
  character(len=10) :: numbers
  character(len=10) :: elchrg
  character(len=2)  :: chrg
  character(len=1)  :: chrg2

  numbers = '0123456789'
  chrg = '+-'
  elchrg = ''
  chrg2 = ''
  el = ''

  !write(*,*)iname

  slen = len_trim(iname)
  do i = 1,slen
    sig = iname(i:i)
    !write(*,*) sig
    if (sig == '') cycle
    if (index(numbers,sig) .ne. 0) then
      elchrg = trim(elchrg)//sig
    elseif (index(chrg,sig) .ne. 0) then
      chrg2 = sig
    else
      el = trim(el)//sig
      !write(*,*) el
    end if
  end do

  prot%swat = e2i(el)

  if (elchrg .ne. '') then
    read (elchrg,*) prot%swchrg
  elseif (chrg2 == '+') then
    prot%swchrg = 1
  elseif (chrg2 == '-') then
    prot%swchrg = -1
  else
    prot%swchrg = 0
  end if
  if (chrg2 == '-'.and.prot%swchrg .gt. 0) prot%swchrg = prot%swchrg*(-1)

  if (prot%swat .ne. 0.and.prot%swat .le. 86) then
    write (*,'(2x,a,1x,a,1x,a,a,1x,i0,a)') '-swel :','using',trim(i2e(prot%swat,'nc')), &
    & '-atom with charge',prot%swchrg,' instead of H⁺'

    prot%swelem = .true.

  end if

end subroutine swparse

!----------------------------------------------------!
! for every structure calculate an correction
! to the acid/base reaction
!----------------------------------------------------!
subroutine prot_correction(env,iname)
  use crest_parameters
  use crest_data
  use strucrd
  implicit none
  type(systemdata) :: env
  character(len=*) :: iname
  integer :: nat,nall
  integer,allocatable :: at(:)
  real(wp),allocatable :: xyz(:,:,:)
  real(wp),allocatable :: eread(:)
  integer :: i
  real(wp) :: dE
  real(wp) :: acidchrg
  real(wp) :: d1,d2,d3,d4,d5,d6

  write (*,'(1x,a)') 'Calculate acid/base correction ...'
  call rdensembleparam(iname,nat,nall)
  allocate (xyz(3,nat,nall),eread(nall),at(nat))
  call rdensemble(iname,nat,nall,at,xyz,eread)

  acidchrg = env%chrg+1
  do i = 1,nall
    call wrxyz('acid.xyz',nat,at,xyz(:,:,i))
    call acidbase(env,'acid.xyz','coord',acidchrg,.true.,.false.,dE, &
        & .false.,d1,d2,d3,d4,d5,d6)
    !eread(i) = eread(i) - dE
    eread(i) = d1+d3-dE
  end do

  call wrensemble(iname,nat,nall,at,xyz,eread)

  return
end subroutine prot_correction
