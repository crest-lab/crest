!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2022 Stefan Grimme, Philipp Pracht
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
!=========================================================================================!
subroutine preig(io,occ,f,e,istart,norbs)
  use iso_fortran_env,only:wp => real64
  implicit none
  integer,intent(in) :: io
  integer,intent(in) :: istart
  real(wp),intent(in) :: f
  integer,intent(in) :: norbs
  real(wp),intent(in) :: occ(norbs)
  real(wp),intent(in) :: e(norbs)
  integer :: n,ntimes,nrest,j,n2,k,i
  character(len=*),parameter :: fmt1 = '(" #    : ",2x,8(3x,i6,2x))'
  character(len=*),parameter :: fmt2 = &
  &       '(" occ. : ",2x,8(4x,f6.3,1x))'
  character(len=*),parameter :: fmt3 = '(" eps  : ",2x,8f11.3)'

  write (io,'(/,10x,''eigenvalues'')')
  n = 8
  ntimes = (norbs - istart + 1) / n
  nrest = mod(norbs - istart + 1,n)
  if (ntimes .eq. 0) nrest = norbs - istart + 1

  j = istart
  n2 = n + istart - 1

  do k = 1,ntimes
    write (io,fmt1) (i,i=j,n2)
    write (io,fmt2) (occ(i),i=j,n2)
    write (io,fmt3) (f * e(i),i=j,n2)
    j = j + n
    n2 = n2 + n
  end do

  if (nrest .gt. 0 .or. ntimes .eq. 0) then
    write (io,fmt1) (i,i=j,j + nrest - 1)
    write (io,fmt2) (occ(i),i=j,j + nrest - 1)
    write (io,fmt3) (f * e(i),i=j,j + nrest - 1)
  end if

  return
end subroutine preig
!=========================================================================================!
subroutine PREIG2(IO,OCC,E,NORBS)
  use iso_fortran_env,only:wp => real64
  implicit none
  integer,intent(in) :: io
  integer,intent(in) :: norbs
  real(wp),intent(in) :: occ(norbs)
  real(wp),intent(in) :: e(norbs)
  integer :: n,ntimes,nrest,j,n2,k,i
  character(len=*),parameter :: fmt1 = &
  & '("#       :",2x,6(4x,i4,2x))'
  character(len=*),parameter :: fmt2 = &
  & '("# atoms :",2x,6(4x,f5.3,1x))'
  character(len=*),parameter :: fmt3 = &
  & '("shift ev:",2x,6f10.5)'

  n = 6
  ntimes = norbs / n
  nrest = mod(norbs,n)
  if (ntimes .eq. 0) nrest = norbs

  j = 1
  n2 = n

  do k = 1,ntimes
    write (io,fmt1) (i,i=j,n2)
    write (io,fmt2) (occ(i),i=j,n2)
    write (io,fmt3) (e(i),i=j,n2)
    j = j + n
    n2 = n2 + n
  end do

  if (nrest .gt. 0 .or. ntimes .eq. 0) then
    write (io,fmt1) (i,i=j,j + nrest - 1)
    write (io,fmt2) (occ(i),i=j,j + nrest - 1)
    write (io,fmt3) (e(i),i=j,j + nrest - 1)
  end if

  return
end subroutine preig2
!=========================================================================================!
subroutine PREIG3(IO,E,NORBS)
  IMPLICIT REAL * 8(A - H,O - Z)
  DIMENSION E(*)
  character(len=*),parameter :: fmt3 = &
  &     '(" value",i5,"-",i5,":",2X,12F6.2)'
  N = 10
  NTIMES = NORBS / N
  NREST = MOD(NORBS,N)
  if (NTIMES .EQ. 0) NREST = NORBS
  J = 1
  N2 = N
  do K = 1,NTIMES
    write (IO,fmt3) J,N2, (E(i),I=J,N2)
    J = J + N
    N2 = N2 + N
  end do
  if (NREST .GT. 0 .OR. NTIMES .EQ. 0) THEN
    write (IO,fmt3) J,J + NREST - 1, (E(i),I=J,J + NREST - 1)
  end if
end subroutine PREIG3

!=========================================================================================!
! SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
! TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
! ((N+1)*N)/2 WHEN M IS ZERO
subroutine PRMAT(IUOUT,R,N,M,HEAD)
  real*8,dimension(*) :: R
  character(len=*) :: HEAD
  character(len=*),parameter :: fmt1 = &
  & '(/," MATRIX PRINTED:",2x,a)'
  character(len=*),parameter :: fmt2 = '(/,5x,6(3x,i4,3x),/)'
  character(len=*),parameter :: fmt3 = '(1x,i4,6f10.5)'

  write (IUOUT,fmt1) HEAD
  NKPB = 6
  if (M <= 0) then
    IBL = N / NKPB
    IR = N - IBL * NKPB
    J1 = 1
    K1S = 1
    KD = 0
    if (IBL .ne. 0) then
      J2 = NKPB
      do I = 1,IBL
        write (IUOUT,fmt2) (J,J=J1,J2)
        K1 = K1S
        K2 = K1
        KK = 0
        do J = J1,J2
          write (IUOUT,fmt3) J, (R(K),K=K1,K2)
          KK = KK + 1
          K1 = K1 + KD + KK
          K2 = K1 + KK
        end do
        J1 = J1 + NKPB
        if (J1 .GT. N) return
        J2 = J2 + NKPB
        K2 = K1 - 1
        K1 = K2 + 1
        K2 = K1 + (NKPB - 1)
        K1S = K2 + 1
        KK = KD + NKPB
        do J = J1,N
          write (IUOUT,fmt3) J, (R(K),K=K1,K2)
          KK = KK + 1
          K1 = K1 + KK
          K2 = K2 + KK
        end do
        KD = KD + NKPB
      end do
    end if
    if (IR .ne. 0) then
      K1 = K1S
      J2 = J1 + IR - 1
      KK = 0
      K2 = K1
      write (IUOUT,fmt2) (J,J=J1,J2)
      write (IUOUT,*)
      do J = J1,J2
        write (IUOUT,fmt3) J, (R(K),K=K1,K2)
        KK = KK + 1
        K1 = K1 + KD + KK
        K2 = K1 + KK
      end do
    end if
    return
  end if

  IBL = M / NKPB
  IR = M - IBL * NKPB
  I2 = 0
  K2 = 0
  if (IBL .ne. 0) then
    do I = 1,IBL
      I1 = (I - 1) * N * NKPB + 1
      I2 = I1 + (NKPB - 1) * N
      K1 = K2 + 1
      K2 = K1 + (NKPB - 1)
      write (IUOUT,fmt2) (K,K=K1,K2)
      do J = 1,N
        write (IUOUT,fmt3) J, (R(IJ),IJ=I1,I2,N)
        I1 = I1 + 1
        I2 = I1 + (NKPB - 1) * N
      end do
    end do
  end if
  if (IR .ne. 0) then
    I1 = IBL * N * NKPB + 1
    I2 = I1 + (IR - 1) * N
    K1 = K2 + 1
    K2 = M
    write (IUOUT,fmt2) (K,K=K1,K2)
    write (IUOUT,*)
    do J = 1,N
      write (IUOUT,fmt3) J, (R(IJ),IJ=I1,I2,N)
      I1 = I1 + 1
      I2 = I1 + (IR - 1) * N
    end do
  end if
  write (IUOUT,*)
  return
end subroutine prmat

!=========================================================================================!
! SUBROUTINE PRINTS MATRIX RR,WHICH IS SUPPOSED
! TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
! ((N+1)*N)/2 WHEN M IS ZERO
subroutine PRMATI(IUOUT,RR,N,M,HEAD)
  character(len=*) :: HEAD
  integer RR(*)
  DIMENSION R(n * n)
  character(len=*),parameter :: fmt1 = &
  &       '(/," MATRIX PRINTED:",2x,a)'
  character(len=*),parameter :: fmt2 = '(/,5x,6(3x,i4,3x),/)'
  character(len=*),parameter :: fmt3 = '(1x,i4,6f10.5)'

  R(1:n * n) = float(RR(1:n * n))
  write (IUOUT,fmt1) HEAD
  NKPB = 6
  if (M <= 0) then
    IBL = N / NKPB
    IR = N - IBL * NKPB
    J1 = 1
    K1S = 1
    KD = 0
    if (IBL .ne. 0) then
      J2 = NKPB
      do I = 1,IBL
        write (IUOUT,fmt2) (J,J=J1,J2)
        K1 = K1S
        K2 = K1
        KK = 0
        do J = J1,J2
          write (IUOUT,fmt3) J, (R(K),K=K1,K2)
          KK = KK + 1
          K1 = K1 + KD + KK
          K2 = K1 + KK
        end do
        J1 = J1 + NKPB
        if (J1 .GT. N) return
        J2 = J2 + NKPB
        K2 = K1 - 1
        K1 = K2 + 1
        K2 = K1 + (NKPB - 1)
        K1S = K2 + 1
        KK = KD + NKPB
        do J = J1,N
          write (IUOUT,fmt3) J, (R(K),K=K1,K2)
          KK = KK + 1
          K1 = K1 + KK
          K2 = K2 + KK
        end do
        KD = KD + NKPB
      end do
    end if
    if (IR .ne. 0) then
      K1 = K1S
      J2 = J1 + IR - 1
      KK = 0
      K2 = K1
      write (IUOUT,fmt2) (J,J=J1,J2)
      write (IUOUT,*)
      do J = J1,J2
        write (IUOUT,fmt3) J, (R(K),K=K1,K2)
        KK = KK + 1
        K1 = K1 + KD + KK
        K2 = K1 + KK
      end do
    end if
    return
  end if

  IBL = M / NKPB
  IR = M - IBL * NKPB
  I2 = 0
  K2 = 0
  if (IBL .ne. 0) then
    do I = 1,IBL
      I1 = (I - 1) * N * NKPB + 1
      I2 = I1 + (NKPB - 1) * N
      K1 = K2 + 1
      K2 = K1 + (NKPB - 1)
      write (IUOUT,fmt2) (K,K=K1,K2)
      do J = 1,N
        write (IUOUT,fmt3) J, (R(IJ),IJ=I1,I2,N)
        I1 = I1 + 1
        I2 = I1 + (NKPB - 1) * N
      end do
    end do
  end if
  if (IR .ne. 0) then
    I1 = IBL * N * NKPB + 1
    I2 = I1 + (IR - 1) * N
    K1 = K2 + 1
    K2 = M
    write (IUOUT,fmt2) (K,K=K1,K2)
    write (IUOUT,*)
    do J = 1,N
      write (IUOUT,fmt3) J, (R(IJ),IJ=I1,I2,N)
      I1 = I1 + 1
      I2 = I1 + (IR - 1) * N
    end do
  end if
  write (IUOUT,*)
  return
end subroutine PRMATI

!========================================================================================!
subroutine preigf(io,e,norbs)
  implicit real * 8(a - h,o - z)
  dimension e(*)
  n = 6
  ntimes = norbs / n
  nrest = mod(norbs,n)
  if (ntimes .eq. 0) nrest = norbs
  j = 1
  n2 = n
  do k = 1,ntimes
    write (io,'("eig.   : ",2x,10f9.2)') (e(i),i=j,n2)
    j = j + n
    n2 = n2 + n
  end do
  if (nrest .gt. 0 .or. ntimes .eq. 0) then
    write (io,'("eig.   : ",2x,10f9.2)') (e(i),i=j,j + nrest - 1)
  end if
end subroutine preigf
!========================================================================================!
subroutine preigf0(io,e,norbs)
  implicit real * 8(a - h,o - z)
  dimension e(*)
  n = 6
  ntimes = norbs / n
  nrest = mod(norbs,n)
  if (ntimes .eq. 0) nrest = norbs
  j = 1
  n2 = n
  do k = 1,ntimes
    write (io,'("eig.   : ",2x,10f9.2)') (e(i),i=j,n2)
    j = j + n
    n2 = n2 + n
  end do
  if (nrest .gt. 0 .or. ntimes .eq. 0) then
    write (io,'("eig.   : ",2x,10f9.2)') (e(i),i=j,j + nrest - 1)
  end if
end subroutine preigf0
!=========================================================================================!
