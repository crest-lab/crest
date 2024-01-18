! This file is part of xtb, modified for crest
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

Subroutine get_schoenflies(n,iat,xyz,sfsym,paramar)
  Use iso_c_binding
  Implicit None
  integer,parameter :: wp = selected_real_kind(15,307)

  Interface c_interface
    !Interface to c routine for symmetry recognition
    !attypes are Atom types as integers (e.g 6 for Carbon etc...)
    !coord must be ``one dimensional'' sequential(!) arrays of doubles
    !symbol is the recognized schoenflies symbol
    Subroutine schoenflies(natoms,attypes,coord,symbol,paramar) &
          &    bind(C,name="schoenflies")
      Use iso_c_binding
      import
      Implicit None
      Integer(c_int),Intent(In),value :: natoms
      integer(c_int),intent(in) :: attypes(*)
      real(c_double),intent(in) :: coord(3,*)
      Character(kind=c_char),Intent(out)  :: symbol(*)
      real(c_double),intent(in) :: paramar(*)
    End Subroutine schoenflies
  End Interface c_interface

  !Dummy Arguments
  Character(Len=*) :: sfsym
  Integer,Intent(In)  :: n
  Integer,Intent(In)  :: iat(n)
  Real(wp),Intent(In) :: xyz(3,n)
  Real(wp),Intent(In) :: paramar(11)

  !local variables for passing to c routine:
  Integer(c_int) :: natoms
  Integer(c_int),Allocatable,Dimension(:) :: attypes
  Real(c_double),Allocatable,Dimension(:,:) :: coord
  Real(c_double),Allocatable,Dimension(:) :: c_paramar
  character(kind=c_char) :: symbol(6)

  !local stack:
  Integer :: i

  Allocate (attypes(n))
  Allocate (coord(3,n))
  Allocate (c_paramar(11))

  !now, copy contents
  natoms = n
  attypes = iat
  coord = xyz
  c_paramar = paramar
  symbol = C_NULL_CHAR

  Call schoenflies(natoms,attypes,coord,symbol,c_paramar)

  sfsym = ""
  do i = 1,size(symbol)
    if (symbol(i) .eq. c_null_char) exit
    sfsym(i:i) = symbol(i)
  end do

  !deallocate arrays:
  Deallocate (attypes,coord,c_paramar)
End Subroutine get_schoenflies

subroutine getsymmetry2(pr,iunit,n,iat,xyz,symthr,maxatdesy,sfsym)
  use iso_c_binding,only:c_char,c_null_char
  use iso_fortran_env,only:wp => real64
  implicit none
  integer,intent(in) :: iunit
  integer n,iat(n),maxatdesy
  real(wp) xyz(3,n)
  real(wp) symthr
  Character(len=*) sfsym
  logical pr
  Character(len=4) atmp

  Real(wp) :: paramar(11)  !parameter array for get_schoenflies_

  if (n .gt. maxatdesy) then
    if (pr) write (iunit,*) 'symmetry recognition skipped because # atoms >',maxatdesy
    sfsym = 'none'
    return
  end if

  if (pr) write (iunit,'(a)')
  !parameters for symmetry recognition:
  paramar(1) = -1         ! verbose, increase for more detailed output (to stdout)
  paramar(2) = 10          ! MaxAxisOrder
  paramar(3) = 100         ! MaxOptCycles
  paramar(4) = 0.001d0     ! ToleranceSame
  paramar(5) = 0.5d0       ! TolerancePrimary
  paramar(6) = symthr      ! ToleranceFinal, THIS IS THE IMPORTANT VALUE
  paramar(7) = 0.5d0       ! MaxOptStep
  paramar(8) = 1.0D-7      ! MinOptStep
  paramar(9) = 1.0D-7      ! GradientStep
  paramar(10) = 1.0D-8     ! OptChangeThreshold
  paramar(11) = 5          ! OptChangeHits

  atmp = '    '
  Call get_schoenflies(n,iat,xyz,atmp,paramar)
  !call flush(iunit)

  !TM stuff (trafo table)
  sfsym(1:3) = atmp(1:3)
  if (sfsym(1:1) .eq. 'D') sfsym(1:1) = 'd'
  if (sfsym(1:1) .eq. 'C') sfsym(1:1) = 'c'
  if (sfsym(1:1) .eq. 'T') sfsym(1:1) = 't'
  if (sfsym(1:1) .eq. 'O') sfsym(1:1) = 'o'
  if (sfsym(1:1) .eq. 'I') sfsym(1:1) = 'i'
  if (sfsym(1:1) .eq. 'S') sfsym(1:1) = 's'
  if (sfsym .eq. 'dih') sfsym = 'd6h'
  if (sfsym .eq. 'civ') sfsym = 'c6v'
  if (sfsym(3:3) .gt. 'v'.or.sfsym(3:3) .lt. 'a') sfsym(3:3) = ' '

  if (pr) then
    write (iunit,'(a3,'' symmetry found (for desy threshold: '',e9.2,'')'')') sfsym,symthr
  end if

End subroutine getsymmetry2

