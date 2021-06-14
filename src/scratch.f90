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

!===================================================================================================!
!c scratch dir handling
!===================================================================================================!
subroutine scrdir(env)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod

      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
      !type(options) :: opt       ! MAIN STORAGE OF BOOLEAN SETTINGS

      integer :: ich,io

      if(len_trim(env%scratchdir).lt.1)then
         !call system('mktemp -d > tmpconf 2>/dev/null')
         call execute_command_line('mktemp -d > tmpconf 2>/dev/null', exitstat=io)
         open(newunit=ich,file='tmpconf')
         read(ich,'(a)',iostat=io) env%scratchdir
         if(io < 0 ) then   ! if mktemp failed and tmpconf is empty
            env%scratchdir=''
            return
         endif
         close(ich,status='delete')
      endif

      write(*,'(a,a)') 'Using scratch directory: ',trim((env%scratchdir))

      io = makedir(trim(env%scratchdir))

      if (env%crestver .eq. crest_solv) then
         call copy('solute',trim(env%scratchdir)//'/'//'solute')
         call copy('solvent',trim(env%scratchdir)//'/'//'solvent')
      else
         call copy('coord',trim(env%scratchdir)//'/'//'coord')
      end if
      call copy('.CHRG',trim(env%scratchdir)//'/'//'.CHRG')
      call copy('.UHF',trim(env%scratchdir)//'/'//'.UHF')
      call copy(env%fixfile,trim(env%scratchdir)//'/'//trim(env%fixfile))
      call copy(env%ensemblename,trim(env%scratchdir)//'/'//trim(env%ensemblename))
      call copy(env%constraints,trim(env%scratchdir)//'/'//trim(env%constraints))
      call copy(trim(env%fixfile),trim(env%scratchdir)//'/'//trim(env%fixfile))

!      io = sylnk('./scratch',trim(env%scratchdir))
      io = sylnk(trim(env%scratchdir),'./scratch')

      call chdir(trim(env%scratchdir))

end subroutine scrdir

subroutine scrend(env)
      use iso_fortran_env, wp => real64
      use crest_data
      use iomod

      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA
      !type(options) :: opt       ! MAIN STORAGE OF BOOLEAN SETTINGS

      character(len=1024) :: crefi,crefi2
      integer :: io
      logical :: ex

      if(len_trim(env%scratchdir).lt.1)then

         return
      endif

      call copy(trim(env%scratchdir)//'/'//'coord','coord')
      call copy(trim(env%scratchdir)//'/'//conformerfile,conformerfile)
      call checkname_xyz(trim(env%scratchdir)//'/'//crefile,crefi,crefi2)
      call copy(trim(crefi),crefile//'.xyz')

      inquire(file=trim(env%scratchdir)//'/'//conformerfilebase//'.sdf',exist=ex)
      if(ex)then
        call copy(trim(env%scratchdir)//'/'//conformerfilebase//'.sdf',conformerfilebase//'.sdf')
      endif

      inquire(file=trim(env%scratchdir)//'/'//'crest_ensemble.xyz',exist=ex)
      if(ex)then
        call copy(trim(env%scratchdir)//'/'//'crest_ensemble.xyz','crest_ensemble.xyz')
      endif

      call system('cp -r '//trim(env%scratchdir)//'/* ./')

      if(.not.env%keepScratch)then
      call rmrf(env%scratchdir)
      endif
 
      return 
end subroutine scrend
