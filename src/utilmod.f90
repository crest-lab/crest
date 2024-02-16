!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2023 Stefan Grimme, Philipp Pracht
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

module utilities
  use crest_parameters
  use ls_rmsd
  implicit none
  private

  !> subroutines
  public :: boltz2
  public :: boltz
  public :: checkname_tmp
  public :: checkname_xyz
  public :: distance
  public :: getanmrrc
  public :: getname1
  public :: heavyrmsd
  public :: rdarg
  public :: revlin_count
  public :: revlin
  public :: TRJappendto_skipfirst
  public :: XYZappendto
  public :: dumpenergies

  !> functions
  public :: lin
  public :: lina
  public :: linr
  public :: ohbonded2
  public :: ohbonded
  public :: distcma

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

!*****************************************
!* "lin" routines address in packed array
!* formerly in "lin.f"
!*****************************************
  integer function lin(i1,i2)
    integer ::  i1,i2,idum1,idum2
    idum1 = max(i1,i2)
    idum2 = min(i1,i2)
    lin = idum2+idum1*(idum1-1)/2
    return
  end function lin
!==========================================!
  subroutine revlin(k,i,j)
!******************************************
!* determine i and j for a given index k
!* (i.e., the reverse of the lin function)
!******************************************
    implicit none
    integer(int64) :: k,b,x
    integer :: i,j
    real(wp) :: kf
    real(wp) :: idum
    kf = real(k,8)
    idum = 1.0+8.0*kf
    !idum = sqrt(idum)
    idum = idum**0.5
    idum = (1.0+idum)/2.0
    i = floor(idum)
    x = i-1
    x = x*i
    b = x/2
    b = k-b
    j = int(b)
    return
  end subroutine revlin
!==========================================!
  subroutine revlin_count(k,i,j,dimij)
!********************************************
!* determine i and j for a given index k
!* by counting until k is reached
!* (only for testing, because too expensive)
!********************************************
    implicit none
    integer(int64) :: k
    integer :: i,j
    integer :: a,b
    integer :: dimij
    integer(int64) :: kdum
    OUTER: do a = 1,dimij
      do b = 1,a
        kdum = lina(a,b)
        if (kdum == k) then
          i = a
          j = b
          exit OUTER
        end if
      end do
    end do OUTER
    return
  end subroutine revlin_count
!==========================================!
  integer(int64) function lina(i1,i2)
!*********************************
!* int64 version of lin function
!*********************************
    integer :: i1,i2,idum1,idum2
    idum1 = max(i1,i2)
    idum2 = min(i1,i2)
    lina = idum1
    lina = lina*(idum1-1)
    lina = lina/2
    lina = lina+idum2
    return
  end function lina
!==========================================!
  integer(int64) function linr(o1,o2,i)
    integer(int64) :: o1
    integer   :: o2,i
    linr = (o1+1)
    linr = linr+i
    linr = linr-o2
    return
  end function linr

!========================================================================================!

  subroutine boltz(n,t,e,p)
!*********************************
!* Boltzmann weighting routines
!* formerly in "boltz.f"
!*********************************
    integer,intent(in) :: n
    real(wp) :: e(*),p(*)
    real(wp),allocatable :: e2(:)
    real(wp) :: t,f,hsum,esum
    integer :: i
    allocate (e2(n))
    !f=8.314*t/4.184d+3
    f = 0.593d0/298.15d0
    f = f*t
    esum = 0
    do i = 1,n
      !e2(i) = (e(i) -e(1))* 627.5095
      e2(i) = e(i)
    end do

    do i = 1,n
      esum = esum+exp(-e2(i)/f)
    end do
    hsum = 0
    do i = 1,n
      p(i) = exp(-e2(i)/f)/esum
    end do
    deallocate (e2)
  end subroutine boltz

  subroutine boltz2(n,e,p)
    implicit none

    integer,intent(in)   :: n     ! Number of molecules
    real(wp),intent(in)  :: e(n)  ! Molecule energies
    real(wp),intent(out) :: p(n)  ! Population

    !> Boltzman constant k in Eh/K
    real(wp),parameter :: kh = 3.1668114d-6
    !> Room temperature
    real(wp),parameter :: T = 298.15_wp

    real(wp) :: val
    real(wp) :: denom
    real(wp) :: emin
    integer  :: i

    p = 0.0_wp
    emin = minval(e)

    do i = 1,n
      val = -(e(i)-emin)/(kh*T)
      p(i) = exp(val)
    end do
    denom = sum(p)
    p = p/denom
  end subroutine boltz2

!=========================================================================================!

  subroutine heavyrmsd(n,nall,k,l,xyz,at,rmsdval)
!*********************************
!* Calculate heavy atom (+OH) RMSD
!**********************************
    use ls_rmsd
    implicit none
    integer n,at(n),j,nall,k,l,nn
    real(wp) xyz(3,n,nall),rmsdval

    logical ::  oh

    real(wp),allocatable :: xyz1(:,:),xyz2(:,:)
    !Dummys:
    real(wp) g(3,3),U(3,3),x_center(3),y_center(3)
    integer i

    nn = 0
    do j = 1,n
      oh = ohbonded2(n,j,xyz(1,1,k),at)
      if (at(j) .gt. 2.or.oh) nn = nn+1
    end do
    allocate (xyz1(3,nn),xyz2(3,nn))

    i = 0
    do j = 1,n
      oh = ohbonded2(n,j,xyz(1,1,k),at)
      if (at(j) .gt. 2.or.oh) then
        i = i+1
        xyz1(1:3,i) = xyz(1:3,j,k)*bohr
        xyz2(1:3,i) = xyz(1:3,j,l)*bohr
      end if
    end do

    call rmsd(i,xyz1,xyz2,0,U,x_center,y_center,rmsdval,.false.,g)

    deallocate (xyz1,xyz2)
  end subroutine heavyrmsd

  logical function ohbonded2(n,m,xyz,at)
    integer :: n,at(n),m,i
    real(wp) :: xyz(3,n)
    real(wp) :: r
    ohbonded2 = .false.
    if (at(m) .ne. 1) return

    do i = 1,n
      if (i .eq. m.or.at(i) .ne. 8) cycle
      r = sqrt((xyz(1,i)-xyz(1,m))**2  &
      &       +(xyz(2,i)-xyz(2,m))**2  &
      &       +(xyz(3,i)-xyz(3,m))**2)
      if (r*0.52917726d0 .lt. 1.1) then
        ohbonded2 = .true.
        exit
      end if
    end do

  end function ohbonded2

!--- formerly in "ohbonded.f"
  logical function ohbonded(n,m,xyz,at,acid)
    integer :: n,at(n),m,acid(86),i
    real(wp) :: xyz(3,n)
    real(wp) :: r

    ohbonded = .false.
    if (at(m) .ne. 1) return

    do i = 1,n
      if (i .eq. m) cycle
      if (acid(at(i)) .eq. 0) cycle
      r = sqrt((xyz(1,i)-xyz(1,m))**2  &
      &       +(xyz(2,i)-xyz(2,m))**2  &
      &       +(xyz(3,i)-xyz(3,m))**2)
      if (r*0.52917726d0 .lt. 1.2) then
        ohbonded = .true.
        exit
      end if
    end do

  end function ohbonded

!========================================================================================!

  subroutine distance(n,xyz,r)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(inout) :: r(:,:)
    real(wp),intent(in) :: xyz(3,n)
    real(wp) :: dx,dy,dz
    integer :: i,j
    do i = 1,n
      do j = 1,n
        dx = xyz(1,j)-xyz(1,i)
        dy = xyz(2,j)-xyz(2,i)
        dz = xyz(3,j)-xyz(3,i)
        r(j,i) = sqrt(dx*dx+dy*dy+dz*dz)
      end do
      r(i,i) = 0
    end do
  end subroutine distance

  real(wp) function distcma(n,j,xyz)
    implicit none
    integer,intent(in) :: n,j
    real(wp),intent(in) :: xyz(3,n)
    real(wp) :: dx,dy,dz
    dx = xyz(1,j)
    dy = xyz(2,j)
    dz = xyz(3,j)
    distcma = sqrt(dx*dx+dy*dy+dz*dz)
  end function distcma

!========================================================================================!

  subroutine checkname_tmp(base,fname,checkname)
!************************************************
!* Iterate X through files called <base>_<X>.tmp
!* And returns the file names
!* <base>_<X>.tmp as <fname>
!* and <base>_<X+1>.tmp as <checkname>
!***********************************************
    character(len=*) :: base,fname,checkname
    integer :: i,j
    logical :: ex
    i = 0
    do
      write (checkname,'(a,''_'',i0,''.tmp'')') trim(base),i
      inquire (file=trim(checkname),exist=ex)
      if (ex) then
        i = i+1
      else
        exit
      end if
    end do
    j = max(0,i-1)
    write (fname,'(a,''_'',i0,''.tmp'')') trim(base),j
  end subroutine checkname_tmp

!===============================================================!
  subroutine checkname_xyz(base,fname,checkname)
!************************************************
!* Iterate X through files called <base>_<X>.xyz
!* And returns the file names
!* <base>_<X>.xyz as <fname>
!* and <base>_<X+1>.xyz as <checkname>
!************************************************
    character(len=*) :: base,fname,checkname
    integer :: i,j
    logical :: ex
    i = 0
    do
      write (checkname,'(a,''_'',i0,''.xyz'')') trim(base),i
      inquire (file=trim(checkname),exist=ex)
      if (ex) then
        i = i+1
      else
        exit
      end if
    end do
    j = max(0,i-1)
    write (fname,'(a,''_'',i0,''.xyz'')') trim(base),j
  end subroutine checkname_xyz

!========================================================================================!

  subroutine getname1(i,atmp)
!*****************************************************
!* generate file name scoord.<i> and return as <atmp>
!*****************************************************
    integer :: i
    character(len=*) :: atmp
    write (atmp,'(''scoord.'',i0)') i
  end subroutine getname1

!========================================================================================!

  subroutine getanmrrc(atmp,fail)
!*************************************
!* check for existance of file <atmp>
!* fail==.true. if it does not exist
!*************************************
    implicit none
    character(len=*),intent(inout) :: atmp
    logical,intent(out) :: fail
    logical :: ex
    inquire (file=trim(atmp),exist=ex)
    fail = .not.ex
    return
  end subroutine getanmrrc

!========================================================================================!

! copy a coord file until an $set-block is encountered
  subroutine clear_setblock(fname)
    implicit none
    character(len=*) :: fname
    character(len=512) :: atmp
    integer :: iost
    integer :: ich,ich2

    open (newunit=ich,file=fname)
    open (newunit=ich2,file='.setdgtmp')

    do
      read (ich,'(a)',iostat=iost) atmp
      if (iost < 0) exit
      if ((index(atmp,'$set') .ne. 0).or.  &
      &  (index(atmp,'$end') .ne. 0)) then
        write (ich2,'(a)') '$end'
        exit
      else
        write (ich2,'(a)') trim(atmp)
      end if
    end do
    close (ich,status='delete')
    close (ich2)
    call rename('.setdgtmp',fname)
  end subroutine clear_setblock

!========================================================================================!

  subroutine rdarg(str,arg,val)
!*************************************************************************
!* read a string <str> and get the value <val> (returned as string, also)
!* that was assignet to the argument <arg>
!* Example:
!* <str>= "reference=foo"
!* <arg>= "reference="  --> <val> will return "foo"
!*************************************************************************
    implicit none
    character(len=*) :: str
    character(len=*) :: arg
    character(len=*) :: val
    character(len=512) :: tmp
    integer :: io
    val = ''
    io = index(str,arg,.true.)
    if (io .ne. 0) then
      io = io+len(arg)
      tmp = str(io:)
      val = trim(tmp)
    else
      val = ''
    end if
    return
  end subroutine rdarg

!========================================================================================!

  subroutine TRJappendto_skipfirst(from,to)
!*********************************************************************************
!* append content of test file "from" into text file "to", similar to "cat A >> B"
!* but specifically for TRJ files, but leaves out the first structure on "from"
!*********************************************************************************
    implicit none
    integer :: io,ich,och
    integer :: nat,i
    character(len=*) :: from
    character(len=*) :: to
    character(len=1024) :: str
    open (newunit=ich,file=to)
    open (newunit=och,file=from)
    do
      read (ich,*,iostat=io)
      if (io < 0) then
        backspace (ich)
        exit
      end if
    end do
    !---first structure is read, but not copied
    read (och,'(a)',iostat=io) str
    if (io == 0) then
      read (str,*) nat
      read (och,'(a)',iostat=io) str
      do i = 1,nat
        read (och,'(a)',iostat=io) str
      end do
    end if
    !---------------------------
    do
      read (och,'(a)',iostat=io) str
      if (io < 0) exit
      write (ich,'(a)') trim(str)
    end do
    close (och)
    close (ich)
  end subroutine TRJappendto_skipfirst

!===============================================================!

  subroutine XYZappendto(from,to)
!**********************************************************************************
!* append content of test file "from" into text file "to", similar to "cat A >> B"
!* but specifically for XYZ files
!**********************************************************************************
    implicit none
    integer :: io
    character(len=*) :: from
    character(len=*) :: to
    integer :: i,nat,tunit,funit
    character(len=1024) :: str
    open (newunit=tunit,file=to)
    open (newunit=funit,file=from)
    do
      read (tunit,*,iostat=io)
      if (io < 0) then
        backspace (tunit)
        exit
      end if
    end do
    read (funit,*) nat
    write (tunit,*) nat
    do i = 1,nat+1
      read (funit,'(a)',iostat=io) str
      if (io < 0) exit
      write (tunit,'(a)') trim(str)
    end do
    close (funit)
    close (tunit)
  end subroutine XYZappendto

!============================================================!

  subroutine dumpenergies(filename,eread)
!****************************
!* write energies to a file
!****************************
    implicit none
    character(len=*),intent(in) :: filename
    real(wp),intent(in) :: eread(:)
    integer :: ich,io,l,i

    open (newunit=ich,file=filename)
    l = size(eread,1)
    do i = 1,l
      write (ich,'(f25.15)') eread(i)
    end do
    close (ich)
  end subroutine dumpenergies

!========================================================================================!
!========================================================================================!
end module utilities
