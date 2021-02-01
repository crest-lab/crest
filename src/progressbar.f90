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
subroutine progbar(percent,bar)
      implicit none
      real*8,intent(in) :: percent
      character(len=52),intent(inout) :: bar
      integer :: i,j,k,l
      integer :: done,notdone

      bar='['

      done=nint(percent/2)
      notdone=50-done

      do i=1,done
         bar=trim(bar)//'#'
      enddo


      do i=1,notdone
         bar=trim(bar)//'-'
      enddo

      bar=trim(bar)//']'

end subroutine progbar

subroutine printprogbar(percent,bar)
      implicit none
      real*8,intent(in) :: percent
      character(len=52),intent(in) :: bar

      write(0,FMT="(A1,A52,2x,F6.2,A)",ADVANCE="NO") achar(13), &
      & bar, percent, '% finished.'

      flush(0)

end subroutine printprogbar

subroutine printemptybar()
     implicit none
     real*8 :: percent
     character(len=52) :: bar
 
     percent=0.00d0
     call progbar(percent,bar)
     call printprogbar(percent,bar)

end subroutine printemptybar
