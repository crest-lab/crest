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

subroutine topofunctional(env,fname)
    use crest_data
    use zdata
    implicit none
    type(systemdata) :: env
    !type(options)    :: opt 

    character(len=*) :: fname
    type(zmolecule) :: zmol
    call  simpletopo_file(fname,zmol,.true.,.true.,'none')

    !--- do a xtb SP calculation to get a wbo file
    call xtbsp(env)


    stop
end subroutine topofunctional



function has_benzene(zmol) result(does)
    use zdata
    implicit none
    type(zmolecule) :: zmol
    logical :: does
    logical :: lcount
    integer :: i,j,k,l
    does = .false.
    lcount=.false.


    return
end function has_benzene


