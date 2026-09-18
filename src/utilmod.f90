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
  public :: get_combinations

  !> functions
  public :: lin
  public :: lina
  public :: linr
  public :: ohbonded2
  public :: ohbonded
  public :: distcma
  public :: binomial
  public :: factorial
  public :: nth_prime

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
    character(len=*),intent(in) :: base
    character(len=*),intent(inout) :: fname,checkname
    character(len=256) :: atmp,btmp
    integer :: i,j
    logical :: ex
    i = 0
    !write(*,*) trim(fname)
    !write(*,*) trim(base)
    do
      write (atmp,'(a,''_'',i0,''.xyz'')') trim(base),i
      inquire (file=trim(atmp),exist=ex)
      checkname = trim(atmp)
      if (ex) then
        i = i+1
      else
        exit
      end if
    end do
    j = max(0,i-1)
    write (btmp,'(a,''_'',i0,''.xyz'')') trim(base),j
    fname = trim(btmp)
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

  function binomial(n,k) result(res)
!**************************************************
!* Function to calculate the binomial coefficient
!*
!*  ⎛ n ⎞
!*  ⎝ k ⎠  = n! / (k! * (n - k)!)
!*
!**************************************************
    implicit none
    integer,intent(in) :: n,k
    real(wp) :: reswp
    integer :: res
    reswp = factorial(n)/(factorial(k)*factorial(n-k))
    res = nint(reswp)
  end function binomial

  function factorial(x) result(fact)
!***************************************************
!* Function to calculate the factorial of a number
!* factorial(x) = x! = x * (x-1) * (x-2) * ... * 1
!***************************************************
    implicit none
    integer,intent(in) :: x
    integer :: i
    real(wp) :: fact
    fact = 1.0_wp
    do i = 2,x
      fact = fact*real(i,wp)
    end do
  end function factorial

  recursive subroutine get_combinations(n,k,ntot,c,combinations,tmp,depth)
    implicit none
    integer,intent(in) :: n,k,ntot,depth !> depth should start out as 0
    integer,intent(inout) :: c,tmp(k)
    integer,intent(inout) :: combinations(k,ntot)
    integer :: i
    if (depth >= k) then
      c = c+1
      combinations(:,c) = tmp(:)
      return
    else if (depth == 0) then
      do i = 1,n
        tmp(depth+1) = i
        call get_combinations(n,k,ntot,c,combinations,tmp,depth+1)
      end do
    else
      do i = 1,tmp(depth)
        if (i == tmp(depth)) cycle
        tmp(depth+1) = i
        call get_combinations(n,k,ntot,c,combinations,tmp,depth+1)
      end do
    end if
  end subroutine get_combinations

!========================================================================================!

  function nth_prime(n) result(prime)
!********************************************
!* get the n-th prime number.
!* The first thousand are saved as a
!* parameter to reduce computational effort.
!********************************************
    implicit none
    integer,intent(in) :: n
    integer :: prime
    integer :: c,num,i
    logical :: is_prime
!&<
    integer,parameter :: maxprimepar = 1000
    integer, parameter :: prime_numbers(maxprimepar) = (/          &
    &  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,                         &
    &  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,                     &
    &  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,                &
    &  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,           &
    &  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,           &
    &  233, 239, 241, 251, 257, 263, 269, 271, 277, 281,           &
    &  283, 293, 307, 311, 313, 317, 331, 337, 347, 349,           &
    &  353, 359, 367, 373, 379, 383, 389, 397, 401, 409,           &
    &  419, 421, 431, 433, 439, 443, 449, 457, 461, 463,           &
    &  467, 479, 487, 491, 499, 503, 509, 521, 523, 541,           &
    &  547, 557, 563, 569, 571, 577, 587, 593, 599, 601,           &
    &  607, 613, 617, 619, 631, 641, 643, 647, 653, 659,           &
    &  661, 673, 677, 683, 691, 701, 709, 719, 727, 733,           &
    &  739, 743, 751, 757, 761, 769, 773, 787, 797, 809,           &
    &  811, 821, 823, 827, 829, 839, 853, 857, 859, 863,           &
    &  877, 881, 883, 887, 907, 911, 919, 929, 937, 941,           &
    &  947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,         &
    &  1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
    &  1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
    &  1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, &
    &  1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
    &  1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
    &  1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
    &  1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
    &  1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
    &  1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
    &  1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
    &  1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
    &  1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
    &  1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, &
    &  1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
    &  2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
    &  2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
    &  2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
    &  2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
    &  2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
    &  2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
    &  2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
    &  2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
    &  2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, &
    &  2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
    &  2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
    &  2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
    &  3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
    &  3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
    &  3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
    &  3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
    &  3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
    &  3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
    &  3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, &
    &  3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
    &  3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
    &  3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
    &  3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
    &  3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
    &  4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
    &  4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
    &  4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
    &  4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
    &  4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, &
    &  4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
    &  4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
    &  4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
    &  4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
    &  4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
    &  4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
    &  4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
    &  5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
    &  5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
    &  5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, &
    &  5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
    &  5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
    &  5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
    &  5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
    &  5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
    &  5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
    &  5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
    &  5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
    &  5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
    &  6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, &
    &  6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
    &  6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
    &  6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
    &  6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
    &  6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
    &  6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
    &  6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
    &  6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
    &  6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
    &  6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, &
    &  7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
    &  7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
    &  7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
    &  7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
    &  7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
    &  7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
    &  7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
    &  7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
    &  7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
    &  7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)
!&>
    if (n <= maxprimepar) then
      prime = prime_numbers(n)
      return
    end if
    c = maxprimepar
    num = prime_numbers(maxprimepar)
    do while (c < n)
      num = num+1
      is_prime = .true.
      do i = 2,int(sqrt(real(num)))
        if (mod(num,i) == 0) then
          is_prime = .false.
          exit
        end if
      end do
      if (is_prime) then
        c = c+1
      end if
    end do
    prime = num
  end function nth_prime

!========================================================================================!
!========================================================================================!
end module utilities
