!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Stefan Grimme
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

!========================================================================================!
subroutine readl(a1,x,n)
  use iso_fortran_env,only:wp => real64
  implicit real(wp) (a-h,o-z)
  character(*) a1
  dimension x(*)
  i = 0
  is = 1
10 i = i+1
  x(i) = readaa(a1,is,ib,ie)
  if (ib .gt. 0.and.ie .gt. 0) then
    is = ie
    goto 10
  end if
  n = i-1
  return
end subroutine readl

!========================================================================================!
function readaa(a,istart,iend,iend2)
  use iso_fortran_env,only:wp => real64
  implicit real(wp) (a-h,o-z)
  real(wp) readaa
  character(*) a

  NINE = ICHAR('9')
  IZERO = ICHAR('0')
  MINUS = ICHAR('-')
  IDOT = ICHAR('.')
  ND = ICHAR('D')
  NE = ICHAR('E')
  IBL = ICHAR(' ')

  iend = 0
  iend2 = 0
  idig = 0
  c1 = 0
  c2 = 0
  one = 1.d0
  x = 1.d0
  nl = len(a)
  do j = istart,nl-1
    n = ichar(a(j:j))
    m = ichar(a(j+1:j+1))
    if (n .le. nine.and.n .ge. izero.or.n .eq. idot) goto 20
    if (n .eq. minus.and.(m .le. nine.and.m .ge. izero &
       & .or.m .eq. idot)) goto 20
  end do
  readaa = 0.d0
  return
20 continue
  iend = j
  do i = j,nl
    n = ichar(a(i:i))
    if (n .le. nine.and.n .ge. izero) then
      idig = idig+1
      if (idig .gt. 10) goto 60
      c1 = c1*10+n-izero
    elseif (n .eq. minus.and.i .eq. j) then
      one = -1.d0
    elseif (n .eq. idot) then
      goto 40
    else
      goto 60
    end if
  end do
40 continue
  idig = 0
  do ii = i+1,nl
    n = ichar(a(ii:ii))
    if (n .le. nine.and.n .ge. izero) then
      idig = idig+1
      if (idig .gt. 10) goto 60
      c2 = c2*10+n-izero
      x = x/10
    elseif (n .eq. minus.and.ii .eq. i) then
      x = -x
    else
      goto 60
    end if
  end do
  !
  ! put the pieces together
  !
60 continue
  readaa = one*(c1+c2*x)
  do j = iend,nl
    n = ichar(a(j:j))
    iend2 = j
    if (n .eq. ibl) return
    if (n .eq. nd.or.n .eq. ne) goto 57
  end do
  return

57 c1 = 0.0d0
  one = 1.0d0
  do i = j+1,nl
    n = ichar(a(i:i))
    iend2 = i
    if (n .eq. ibl) goto 70
    if (n .le. nine.and.n .ge. izero) c1 = c1*10.0d0+n-izero
    if (n .eq. minus) one = -1.0d0
  end do
  continue
70 readaa = readaa*10**(one*c1)
  return
end function readaa

!========================================================================================!
!cuts the at blanks and tabstops and returns all floats and strings in order of occurence
subroutine cutline(line,floats,strings)
  use iso_fortran_env,only:wp => real64
  implicit none
  real(wp) floats(*),num
  character(len=128) line,str,stmp
  character(len=80) strings(3)
  character(len=1) digit
  integer i,ty,cs,cf

  stmp = ''
  cs = 1
  cf = 1
  strings(:) = ''
  do i = 1,len(trim(line))
    digit = line(i:i)
    if (digit .ne. ' '.and.digit .ne. char(9)) then  !should exclude tabstops and blanks, 9 is ascii code for tab
      stmp = trim(stmp)//trim(digit)
    elseif (stmp .ne. '') then
      call checktype(stmp,num,str,ty)      !get type of string, 0=number, 1=character
      if (ty .eq. 0) then
        floats(cf) = num
        cf = cf+1
      elseif (ty .eq. 1) then
        strings(cs) = trim(str)
        cs = cs+1
      else
        write (*,*) 'Problem in checktype, must abort'
        exit
      end if
      stmp = ''
    end if
    if (i .eq. len(trim(line))) then  !special case: end of line
      call checktype(stmp,num,str,ty)
      if (ty .eq. 0) then
        floats(cf) = num
        cf = cf+1
      elseif (ty .eq. 1) then
        strings(cs) = trim(str)
        cs = cs+1
      else
        write (*,*) 'Problem in checktype, must abort'
        exit
      end if
      stmp = ''
    end if
  end do
end subroutine cutline

!========================================================================================!
!this checks the type of the string and returns it cast to real or as string.
subroutine checktype(field,num,str,ty)
  use iso_fortran_env,only:wp => real64
  implicit none
  character(len=*) field,str
  real(wp) :: num
  integer :: e,ty
  logical :: is_num

  ty = 99
  str = ''
  is_num = .false.
  read (field,'(F10.5)',IOSTAT=e) num !cast string on real and get error code; 0 means success.
  if (e .eq. 0) is_num = .true.
  if (is_num) then
    if (index(field,'.') .ne. 0) then  !check for integer/real
      read (field,'(F30.16)') num
      ty = 0
    else                       !if integer, add .0 to string; otherwise cast to real does not work
      str = trim(field)//'.0'
      read (str,'(F30.16)') num
      str = ''
      ty = 0
    end if
  else
    str = field
    ty = 1
  end if
end subroutine checktype

