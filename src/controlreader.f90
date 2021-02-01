!================================================================================!
! This file is part of crest.
!
! Copyright (C)2020 Philipp Pracht
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

!=================================================================!
! This is the new control file reader of crest
!=================================================================!
subroutine crest_readcontrolfile(fname,env)
    use crest_data
    use filemod
    use iomod
    implicit none
    character(len=*) :: fname
    type(systemdata) :: env
    !type(options)    :: opt

    type(filetype) :: ctrl
    logical :: ex
    character(len=:),allocatable :: atmp,btmp,ctmp

    !--- FIRST check if we have a file given and if it is
    !    is a crest control file.
    inquire(file=fname,exist=ex)
    if(.not.ex)return
    call ctrl%open(fname)
    call ctrl%clearblanks()
    atmp=trim(ctrl%line(1))
    call to_lower(atmp)
    if(atmp=='$crest')then
        write(*,'(1x,a,1x,a)') 'Reading crest control file:',trim(fname)
    else
        return
    endif

    !--- do some formatting of the control file
    do i=1,ctrl%nlines
       atmp = trim(ctrl%line(i))
       call to_lower(atmp)
       if(atmp=='$end')then
           call ctrl%replace(i,'')
       endif
    enddo
    call ctrl%write('$end')
    call ctrl%clearblanks()




    return
end subroutine crest_readcontrolfile



