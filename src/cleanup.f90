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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C All subroutines related to removing stuff
!C 
!C  system calls are the easiest way to remove files, although they introduce
!C  an additional dependence on the operating system.
!C  But since there are no simple kernel alternatives to do the same thing
!C  I will keep them for now.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!-------------------------------------------------------------------------
! Basic subroutine for removing files or directories
!-------------------------------------------------------------------------
subroutine rmrf(dataname)
     implicit none
     character(len=*) :: dataname
     integer :: io
     !call system('rm -rf '//trim(dataname)//' 2>/dev/null')
     call execute_command_line('rm -rf '//trim(dataname)//' 2>/dev/null', exitstat=io)
     return
end subroutine rmrf
!-------------------------------------------------------------------------
! Basic subroutine for removing files or directories (with wild card ending)
!-------------------------------------------------------------------------
subroutine rmrfw(dataname)
     implicit none
     character(len=*) :: dataname
     call rmrf(trim(dataname)//'*')
     return
end subroutine rmrfw


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C Specific cleanup routines for different parts of the CREST code
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!-------------------------------------------------------------------------
! Cleanup function for MTD (confscript2_misc.f90)
!-------------------------------------------------------------------------
subroutine V2cleanup(restartopt)
      use iomod
      use crest_data
      implicit none
      logical :: restartopt
      character(len=12),parameter :: pipe =' 2>/dev/null'

      call rmrf('.tmpxtbmodef hessian .xtboptok')
      call rmrf('xtbrestart xtbmdok')
      if(.not.restartopt)then
       call rmrfw('METADYN')
       call rmrfw('NORMMD')
       call rmrf('MRMSD')
       call rmrf('OPTIM')
       call rmrf('PROP')
       call rmrf('crest_best.xyz')
       call rmrf(conformerfile)
       call rmrfw(crefile)
       call rmrfw('.cre_')
       call rmrf('cregen_*.tmp')
       call rmrf('gfnff_topo')
      else
       call rmrfw('NORMMD')
       call rename(crefile//'_0.xyz','.tmp_full.xyz')
       call rmrfw(crefile)
       call rename('.tmp_full.xyz',crefile//'_0.xyz')
      endif
      return
end subroutine V2cleanup

!-------------------------------------------------------------------------
! Cleanup function for MTD and MD dirs (confscript2_misc.f90)
!-------------------------------------------------------------------------
subroutine cleanMTD
      implicit none
      call rmrfw('METADYN')
      call rmrfw('NORM')
      call rmrfw('STATICMTD')
      call rmrf('MRMSD')
      return
end subroutine cleanMTD

!------------------------------------------------------------------------
! last bit of file removing and renaming (confscript2_misc.f90)
!------------------------------------------------------------------------
subroutine V2terminating
       use iomod
       use crest_data
       implicit none
       character(len=256) :: inpnam,outnam
       call remove('scoord.1')
       call catdel('cregen.out.tmp') !print to screen
       call checkname_xyz(crefile,inpnam,outnam)
       call rename(inpnam,"crest_rotamers.xyz")
       call rmrfw('crest_rotamers_')
       call rmrfw('crest_smtd_')
       return
end subroutine V2terminating

!------------------------------------------------------------------------
! clean Dir between iterations (confscript2_misc.f90)
!------------------------------------------------------------------------
subroutine clean_V2i
      use iomod
      use crest_data
      implicit none
      call rmrfw('METADYN')
      call rmrfw('NORM')
      call rmrfw('STATICMTD')
      call rmrfw(crefile)
      return
end subroutine clean_V2i

!-----------------------------------------------------------------------
! Remove all crest_conformers and crest_rotamers files (confscript2_misc.f90)
!-----------------------------------------------------------------------
subroutine rmcres()
      use crest_data
      implicit none
      call rmrfw(crefile)
      call rmrf(conformerfile)
      return
end subroutine rmcres

!-----------------------------------------------------------------------
! Remove the OPTIM dir and scoord files (confscript2_misc.f90)
!-----------------------------------------------------------------------
subroutine rmoptim()
      implicit none
      call rmrf('OPTIM')
      call rmrfw('scoord.')
      return
end subroutine rmoptim

!-----------------------------------------------------------------------
! Cleanup options for SCREEN mode (confscript3.f90)
!-----------------------------------------------------------------------
subroutine screen_cleanup
      use crest_data
      implicit none
      call rmrf('OPTIM')
      call rmrf('scoord.1')
      call rmrf(conformerfile)
      call rmrfw(crefile)
      return
end subroutine screen_cleanup

!-----------------------------------------------------------------------
! Cleanup routine for protonation tool (protonate.f90)
!-----------------------------------------------------------------------
subroutine protclean
      implicit none
      call rmrf('xtblmoinfo coordprot.0 lmocent.coord')
      call rmrf('protonate_*.xyz xtbscreen.xyz protonated.xyz')
      return
end subroutine protclean

!-----------------------------------------------------------------------
! Cleanup routine for deprotonation tool (deprotonate.f90)
!-----------------------------------------------------------------------
subroutine deprotclean
      implicit none
      call rmrf('xtblmoinfo coordprot.0 lmocent.coord')
      call rmrf('deprotonate_*.xyz xtbscreen.xyz')
      return
end subroutine deprotclean

!-----------------------------------------------------------------------
! Cleanup routine for tautomerization (tautomerize.f90)
!-----------------------------------------------------------------------
subroutine tautclean
      implicit none
      call rmrf('xtblmoinfo coordprot.0 lmocent.coord')
      call rmrf('protonate_*.xyz deprotonate_*.xyz')
      call rmrf('tautomerize_*.xyz tautomers.xyz xtbscreen.xyz')
      return
end subroutine tautclean
subroutine tautclean2
      implicit none
      call rmrf('xtblmoinfo coordprot.0 lmocent.coord')
      call rmrf('protonate_*.xyz deprotonate_*.xyz')
      call rmrf('protonated.xyz deprotonated.xyz xtbscreen.xyz')
      return
end subroutine tautclean2

!-----------------------------------------------------------------------
! Cleanup routine for MF-MD-GC (confscript1.f90)
!-----------------------------------------------------------------------
subroutine clean
      use crest_data                      
      implicit none
      call rmrfw('MODEF')
      call rmrfw('NORMMD')
      call rmrfw('TMPCONF')
      call rmrfw(crefile)
      return
end subroutine clean



!-----------------------------------------------------------------------
! change the dir and then remove mos and dh (dft_propcalc.f90)
!-----------------------------------------------------------------------
subroutine cleanDFT(TMPCONF)
      implicit none
      integer :: TMPCONF
      character(len=:),allocatable :: str
      character(len=:),allocatable :: pipe
      character(len=10) :: nmmr
      integer :: i,io
      pipe=' 2>/dev/null'
      do i=1,TMPCONF
         write(nmmr,'(i0)')i
         str='cd TMPCONF'//trim(nmmr)//' && rm -rf '
         !call system(str//'mos'//pipe)
         call execute_command_line(str//'mos'//pipe, exitstat=io)
         !call system(str//'dh'//pipe)
         call execute_command_line(str//'dh'//pipe, exitstat=io)
      enddo
      return
end subroutine cleanDFT




