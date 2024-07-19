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

!=========================================================================!
! AUTOMATED DEPROTONATION AND ENERGETIC RANKING SCRIPT
! To use run:
!    crest <input> --deprotonate
!========================================================================!

subroutine deprothead
      implicit none
      write(*,*)'       __________________________________________'
      write(*,*)'      |                                          |'
      write(*,*)'      |      automated deprotonation script      |'
      write(*,*)'      |__________________________________________|'
      write(*,*)' Universitaet Bonn, MCTC'
      write(*,*)' P.Pracht, Wed 28. Nov 13:11:52 CEST 2018'
      write(*,*)
end subroutine deprothead

!--------------------------------------------------------------------------------------------
! Protonation workflow with GFNn-xTB
!--------------------------------------------------------------------------------------------
subroutine deprotonate(env,tim)
      use crest_parameters
      use crest_data
      use iomod
      use strucrd, only: rdnat,rdcoord,i2e
      use utilities
      implicit none
      type(systemdata) :: env
      type(timer)      :: tim
      type(protobj)    :: deprot

      real(wp),allocatable :: xyz(:,:)
      integer,allocatable  :: at(:)

      integer :: i,j
      character(len=64)  :: deprotname
      character(len=256) :: thispath
      character(len=256) :: filename
      character(len=128) :: inpnam,outnam

      integer :: ich
      integer :: natp,nallout,refchrg

!--- printout & clean directory
      call deprotclean
      call deprothead

      if(.not.allocated(env%protb%atmap))allocate(env%protb%atmap(env%nat))
      if(.not.env%protb%strictPDT .and. .not.env%protb%fixPDT)then
!--- sort the input file (H atoms to the bottom)
        call htothebottom('coord',env%chrg,env%nat,env%protb%atmap)
      else
!--- or sort AND apply heavy atom bond constraints
         call PDT_constraints(env)
      endif

!--- get some settings
      call getcwd(thispath)
      deprotname='deprotonate_0.xyz'
      deprot=env%protb
      refchrg = env%chrg
      deprot%newchrg = env%chrg - 1  !increase chrg by one
      env%chrg = env%chrg - 1  !in the new version all calculations access env%chrg!!!
      natp=env%nat - 1 !additional proton, Nat is increased by one
     
!--- write a file with all the possible deprotonated structures
      call tim%start(1,'INPUT generation')
      allocate(xyz(3,env%nat),at(env%nat))
      call rdcoord('coord',env%nat,at,xyz)
      open(newunit=ich,file=deprotname)

      do i=1,env%nat
         if(at(i).ne.1)cycle
         write(ich,'(1x,i6)') natp
         write(ich,*)
         do j=1,env%nat
             if(i.eq.j) then
                 cycle
             else
                 write(ich,'(a2,2x,3F16.10)')i2e(at(j),'nc'),xyz(1:3,j)*bohr
             endif
         enddo
      enddo

      close(ich)
      deallocate(at,xyz)
      call tim%stop(1)
      
!--- get the new charge and set up the calculations
     call tim%start(2,'multilevel OPT')
     open(newunit=ich,file='.CHRG')
     write(ich,*) deprot%newchrg     !new charge written here
     close(ich)
     env%nat=natp  !required so that MDopt_para works

     write(*,'(''-----------------------'')')
     write(*,'(''Multilevel Optimization'')')
     write(*,'(''-----------------------'')')

     call smallhead('1. crude pre-optimization')
     call checkname_xyz('deprotonate',inpnam,outnam)
     call MDopt_para(env,deprotname,1)
        filename=trim(thispath)//'/'//trim(outnam)
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
        call rmrf('OPTIM')
        if(deprot%ABcorrection)then
          call deprot_correction(env,trim(filename))
        endif
        deprot%ewin=deprot%ewin*3.0d0
        call sort_ens(deprot,outnam,.false.)
        call remaining_in(outnam,deprot%ewin,nallout) !--- remaining number of structures
        write(*,*)

     call smallhead('2. loose optimization')
     call checkname_xyz('deprotonate',inpnam,outnam)
     call MDopt_para(env,inpnam,2)
        filename=trim(thispath)//'/'//trim(outnam)
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
        call rmrf('OPTIM')
        if(deprot%ABcorrection)then
          call deprot_correction(env,trim(filename))
        endif
        deprot%ewin=deprot%ewin*(2.0d0/3.0d0)
        call sort_ens(deprot,outnam,.false.)
        call remaining_in(outnam,deprot%ewin,nallout) !--- remaining number of structures
        write(*,*)

     call smallhead('3. optimization with user-defined thresholds')
     call checkname_xyz('deprotonate',inpnam,outnam)
     call MDopt_para(env,inpnam,0)
        filename=trim(thispath)//'/'//trim(outnam)
        call rename('OPTIM'//'/'//'opt.xyz',trim(filename))
        call rmrf('OPTIM')
        if(deprot%ABcorrection)then
          call deprot_correction(env,trim(filename))
        endif
        deprot%ewin=deprot%ewin/2.0d0
        call sort_ens(deprot,outnam,.false.)
        call remaining_in(outnam,deprot%ewin,nallout) !--- remaining number of structures


     call cosort(outnam,'deprotonated.xyz',.false.,.false.)
     call sort_ens(deprot,'deprotonated.xyz',.true.)
     call tim%stop(2)

!>--- (optional) post-processing
     if(env%relax)then
       env%rednat = env%rednat - 1
       call relaxensemble('deprotonated.xyz',env,tim)
     endif

     if(env%outputsdf)then
     call new_wrsdfens(env,'deprotonated.xyz','deprotonated.sdf',.true.)
     endif

!--- reset data for main dir
     env%chrg = refchrg
     if(env%chrg .eq. 0) then
       call remove('.CHRG')
     else
       open(newunit=ich,file='.CHRG')
       write(ich,*) env%chrg
       close(ich)
     endif
     env%nat=natp + 1 !reset nat
     !call chdir(thispath)

     return
end subroutine deprotonate

!----------------------------------------------------!
! for every structure calculate an correction
! to the acid/base reaction
!----------------------------------------------------!
subroutine deprot_correction(env,iname)
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
    write(*,'(1x,a)') 'Calculate acid/base correction ...'
    call rdensembleparam(iname,nat,nall)
    allocate(xyz(3,nat,nall),eread(nall),at(nat))
    call rdensemble(iname,nat,nall,at,xyz,eread)

    acidchrg = env%chrg
    do i=1,nall
     call wrxyz('base.xyz',nat,at,xyz(:,:,i))
     call acidbase(env,'coord','base.xyz',acidchrg,.true.,.false.,dE, &
         &  .false.,d1,d2,d3,d4,d5,d6)
     !eread(i) = eread(i) + dE
     eread(i) = d4+d6+dE
    enddo

    call wrensemble(iname,nat,nall,at,xyz,eread)

    return
end subroutine deprot_correction

