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
!> non-linear Fit by Levenberg-Marquart Algorithm
!> based on a Pascal program by S.G from the early 90s
!> S.G., 10/08
!========================================================================================!
! the function f(x):
real*8 function f_marq(p,x)
  real*8 x,p(*)
  f_marq = p(1) * (1.0d0 - exp(-p(2) * x**p(3)))
end function f_marq
!========================================================================================!
subroutine Marqfit(pr,pr2,n_in,n_use,y_in,nconf_exp,val,rmsd,ok)
  implicit none
  logical :: pr,pr2      !> print flag
  integer :: n_in        !> # points
  integer :: n_use       !> # points used in fit
  real*8 :: y_in(*)      !> data
  real*8 :: nconf_exp    !> fixed parameter
  real*8 :: val          !> extrapolated value
  real*8 :: rmsd         !> of fit
  logical :: ok          !> fit ok ?

  real*8,allocatable :: Nmat(:,:),fx(:),x(:),y(:)
  real*8,allocatable :: rr(:),ll(:),t(:,:),a(:,:),dp(:),p(:)
  real*8  :: s0,s1,dsq,inkre,dp1,eps,sumx,det,f_marq,r2
  integer :: np
  integer :: n
  integer :: i,j,kk,L,iter

!> # parameters
  np = 2
  if (n_in - 1 .gt. n_use) then
    n = n_use
    j = n_in - n_use + 1    !> first value used for fit
  else
    n = n_in - 1
    j = 2
  end if

  allocate (Nmat(np,np),fx(n),p(np + 1),x(n),y(n), &
  &        rr(n),ll(n),t(np,n),a(n,np),dp(np))

  do i = 1,n
    x(i) = float(i)   !> taker iter blocks
  end do
  y(1:n) = y_in(j:n_in) - y_in(1) !> subtract base (MTD=iter0) value

!> initial set of parameters
  p(1) = y(n) + y(n) * 0.1 !> just start values
  p(2) = 0.3               !>  "     "     "
  p(3) = nconf_exp         !> nconf_exp is a fixed parameter

  if (pr) then
    write (*,*) 'initial P:',p(1:np + 1)
    do i = 1,n
      write (*,*) x(i),y(i),f_marq(p,x(i))
    end do
  end if

  eps = 1.d-4  !> exit thr
  inkre = 0.5  !> damping

  s1 = 1.d+42
  do iter = 1,100

    do i = 1,n
      fx(i) = f_marq(p,x(i))
    end do
    s0 = s1
    s1 = 0
    do i = 1,n
      s1 = s1 + (fx(i) - y(i))**2
    end do
    dsq = s1 - s0

!        if(pr)write(*,*) iter,s1,dsq/s1
    if (abs(dsq / s1) .lt. eps .or. s1 .lt. eps) exit

    do kk = 1,np
      dp1 = abs(0.001 * p(kk))
      if (dp1 .lt. 1.d-8) dp1 = 1.d-8
      p(kk) = p(kk) + dp1
      do i = 1,n
        rr(i) = f_marq(p,x(i))
      end do
      p(kk) = p(kk) - 2.0 * dp1
      do i = 1,n
        ll(i) = f_marq(p,x(i))
      end do
      p(kk) = p(kk) + dp1
      do i = 1,n
        a(i,kk) = (rr(i) - ll(i)) / (2.*dp1)
      end do
    end do

    do i = 1,n
      do j = 1,np
        t(j,i) = a(i,j)
      end do
    end do
!> (*Matrizenmultiplikation *)
!> (* T^*A^ = Nmat          *)
    do i = 1,np
      do j = 1,np
        sumx = 0
        do L = 1,n
          sumx = sumx + t(i,L) * a(L,j)
        end do
        Nmat(i,j) = sumx
      end do
    end do

!> damping of diagonal
    do i = 1,np
      Nmat(i,i) = Nmat(i,i) + inkre * Nmat(i,i)
    end do
    inkre = inkre / 1.1
!> (*  T^ * residuen = dp*)
    do i = 1,np
      sumx = 0
      do L = 1,n
        sumx = sumx + T(i,L) * (y(L) - fx(L))
      end do
      dp(i) = sumx
    end do
!> (* Matrix invertieren *)
    call minv(Nmat,np,det)
    if (det .lt. 1.d-12) then
      if (pr) write (*,*) 'no matrix inversion possible. exit.'
      goto 999
    end if
!> (* N^-1 * (transp.*residuen) *)

    do i = 1,np
      sumx = 0
      do L = 1,np
        sumx = sumx + Nmat(i,L) * dp(L)
      end do
!> (* neue Parameter *)
      p(i) = p(i) + sumx
    end do
  end do

999 rmsd = sqrt(s1 / float(n))
  s0 = 0
  sumx = 0
  dsq = 0
  do i = 1,n
    sumx = sumx + (y(i) - f_marq(p,x(i)))**2
    s0 = s0 + abs(y(i) - f_marq(p,x(i))) / (y(i) + 1.d-6)
    dsq = dsq + (y(i) - s0)**2
  end do
  r2 = 1.0d0 - sumx / (dsq + 1.d-9)
  s0 = 100.*s0 / float(n)
  s1 = 100.*abs(p(1) - y(n)) / y(n)

!     if(pr)then
!     do i=1,n
!        write(42,*) i,y(i)
!     enddo
!     write(42,*)
!     do i=1,n
!        write(42,*) i,f_marq(p,x(i))
!     enddo
!     write(42,*)
!     write(42,*) '0 ',p(1)
!     write(42,*) i+1,p(1)
!     endif

  ok = .true.
  if (s0 .gt. 4.00) ok = .false.! result trustworthy if av. rel. dev of fit < 4 %
  if (s1 .gt. 100.0) ok = .false.! extrapolated value very different from last
  if (r2 .lt. 0.95) ok = .false.! bad correlation

  if (pr) then
    write (*,'(''P         :'',6F12.6)') (p(kk),kk=1,np),nconf_exp
    write (*,'(''RMSD      :'', F12.6)') rmsd
    write (*,'(''av % dev  :'', F12.6)') s0
    write (*,'(''% dev last:'', F12.6)') s1
    write (*,'(''ok fit    :'', L    )') ok
  end if

  if (pr2) then
    write (*,'(''   block      Siter   Siter (fit) '')')
    do i = 1,n
      write (*,'(i9,2F9.3)') int(x(i)),y(i),f_marq(p,x(i))
    end do
    write (*,'(''R²       :'', f12.6)') r2
  end if

  val = p(1) + y_in(1) ! add base value

end subroutine marqfit

!========================================================================================!
SUBROUTINE MINV(A,N,D)
  IMPLICIT doUBLE PRECISION(A - H,O - Z)
  DIMENSION A(*)
!**********************************************************************
!*
!*     INVERT A MATRIX USING GAUSS-JORDAN METHOD.  PART OF DIIS
!*     A - INPUT MATRIX (MUST BE A GENERAL MATRIX), DESTROYED IN
!*        COMPUTATION AND REPLACED BY RESULTANT INVERSE.
!*     N - ORDER OF MATRIX A
!*     D - RESULTANT DETERMINANT
!*
!**********************************************************************
  integer M(N),L(N)
!>
!>     SEARCH FOR LARGEST ELEMENT
!>
  D = 1.0D0
  NK = -N
  do K = 1,N
    NK = NK + N
    L(K) = K
    M(K) = K
    KK = NK + K
    BIGA = A(KK)
    do J = K,N
      IZ = N * (J - 1)
      do I = K,N
        IJ = IZ + I
        IF (ABS(BIGA) .LT. ABS(A(IJ))) THEN
          BIGA = A(IJ)
          L(K) = I
          M(K) = J
        END IF
      end do
    end do
!>
!>     INTERCHANGE ROWS
!>
    J = L(K)
    if (J > K) then
      KI = K - N
      do I = 1,N
        KI = KI + N
        HOLD = -A(KI)
        JI = KI - K + J
        A(KI) = A(JI)
        A(JI) = HOLD
      end do
    end if
!>
!>     INTERCHANGE COLUMNS
!>
    I = M(K)
    !IF (I-K) 80,80,60
    if (I > K) then
      JP = N * (I - 1)
      do J = 1,N
        JK = NK + J
        JI = JP + J
        HOLD = -A(JK)
        A(JK) = A(JI)
        A(JI) = HOLD
      end do
    end if
!>
!>     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS
!>     CONTAINED IN BIGA)
!>
    !IF (BIGA) 100,90,100
    if (BIGA == 0.0d0) then
      D = 0.0d0
      RETURN
    end if
    do I = 1,N
      !IF (I-K) 110,120,110
      if (I == K) cycle
      IK = NK + I
      A(IK) = A(IK) / (-BIGA)
    end do
!>  REDUCE MATRIX
    do I = 1,N
      IK = NK + I
      HOLD = A(IK)
      IJ = I - N
      do J = 1,N
        IJ = IJ + N
        !IF (I-K) 130,150,130
        if (I == K) cycle
        !IF (J-K) 140,150,140
        if (J == K) cycle
        KJ = IJ - I + K
        A(IJ) = HOLD * A(KJ) + A(IJ)
      end do
    end do
!>
!>     DIVIDE ROW BY PIVOT
!>
    KJ = K - N
    do J = 1,N
      KJ = KJ + N
      !IF (J-K) 160,170,160
      if (J == K) cycle
      A(KJ) = A(KJ) / BIGA
    end do
!>
!>     PRODUCT OF PIVOTS
!>
    D = MAX(-1.D25,MIN(1.D25,D))
    D = D * BIGA
!>
!>     REPLACE PIVOT BY RECIPROCAL
!>
    A(KK) = 1.0 / BIGA
  end do
!>
!>     FINAL ROW AND COLUMN INTERCHANGE
!>
  K = N
  do while (K > 1)
    K = (K - 1)
    !IF (K) 260,260,200
    I = L(K)
    !IF (I-K) 230,230,210
    if (I > K) then
      JQ = N * (K - 1)
      JR = N * (I - 1)
      do J = 1,N
        JK = JQ + J
        HOLD = A(JK)
        JI = JR + J
        A(JK) = -A(JI)
        A(JI) = HOLD
      end do
    end if
    J = M(K)
    !IF (J-K) 190,190,240
    if (J > K) then
      KI = K - N
      do I = 1,N
        KI = KI + N
        HOLD = A(KI)
        JI = KI - K + J
        A(KI) = -A(JI)
        A(JI) = HOLD
      end do
    end if
  end do
  RETURN
end subroutine MINV
