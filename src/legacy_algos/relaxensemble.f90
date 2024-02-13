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

!==================================================================!
! A function to be used in conjunction with the
! protonation/deprotonation/tautomerization tools
! --> relax a each generated structure by performing a (reduced)
!     conformational search
!==================================================================!
subroutine relaxensemble(fname,env,tim)
      use crest_parameters
      use crest_data
      use strucrd
      use iomod
      use cregen_interface
      implicit none

      character(len=*) :: fname
      type(systemdata) :: env
      type(timer)   :: tim

      integer :: nall
      integer :: nat
      real(wp),allocatable :: xyz(:,:,:)
      real(wp),allocatable :: c0(:,:)
      real(wp) :: energy
      real(wp),allocatable :: eread(:)
      integer,allocatable :: at(:)
      character(len=40),allocatable :: origin(:)

      integer :: i,r,ich

      character(len=:),allocatable :: relaxdir
      character(len=:),allocatable :: oname
      character(len=512) :: thisdir
      character(len=128) :: atmp

!--- deactivate possible constraints of the strict modi
      if(env%ptb%strictPDT .or. env%ptb%fixPDT)then
         env%cts%used=.false.
      endif

      call getcwd(thisdir)
      relaxdir='RELAX'

      !--- read the ensemble
      call rdensembleparam(fname,nat,nall)
      if(nall.lt.1) return
      allocate(xyz(3,nat,nall),at(nat),eread(nall))
      if(.not.env%trackorigin)then
        call rdensemble(fname,nat,nall,at,xyz,eread)
      else
        allocate(origin(nall))
        call rdensemble_origin(fname,nat,nall,at,xyz,eread,origin)
      endif

      write(*,'(/)')
      call largehead('Relaxation of generated structures')
      !--- make the new subdir and change to it
      r = makedir(relaxdir)
      call chdir(relaxdir)

      !---  write .CHRG and .UHR files
      if(env%chrg.ne.0)then
         open(newunit=ich,file='.CHRG')
         write(ich,*)env%chrg
         close(ich)
      endif
      if(env%uhf.ne.0)then
         open(newunit=ich,file='.UHF')
         write(ich,*)env%uhf
         close(ich)
      endif

      allocate(c0(3,nat), source=0.0d0)

      !--- do the reduced conf. search for all the structures
      do i=1,nall
         !--- NOTE: all conformational searches are done in series and the Dir is cleared in between
         call V2cleanup(.false.)
         write(atmp,'(a,i0)') 'Structure ',i
         call largehead(trim(atmp))
         call wrc0('coord',nat,at,xyz(:,:,i)/bohr)
         call confscript2i(env,tim) !MTD-GC algo
         !--- read the new best conformer
         call rdcoord('crest_best.xyz',nat,at,c0,energy)
         c0=c0*bohr
         xyz(:,:,i) = c0(:,:)
         eread(i) = energy
         !stop
      enddo

      call chdir(thisdir)

      !--- write the new ensemble file
      oname='relax.'//trim(fname)
      call wrensemble(oname,nat,nall,at,xyz,eread)

      deallocate(c0)
      if(allocated(origin))deallocate(origin)
      deallocate(eread,at,xyz)
      call rmrf(relaxdir)

      !energy sorting
      call largehead('Sorting relaxed ensemble')
      env%ensemblename=trim(oname)
      env%confgo = .true.
      call newcregen(env,6)
      env%confgo = .false.
      env%ensemblename='none selected'  !RESET, IMPORTANT!
      call rename(trim(oname)//'.sorted',trim(oname))


      return
end subroutine relaxensemble
