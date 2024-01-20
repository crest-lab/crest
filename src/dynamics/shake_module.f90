!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2021 - 2022 Philipp Pracht
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
!
! Routines were adapted from the xtb code (github.com/grimme-lab/xtb)
! under the Open-source software LGPL-3.0 Licencse.
!================================================================================!

module shake_module
  use crest_parameters
  implicit none

  integer,private,parameter :: metal(1:86) = [&
     & 0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1, &
     & 1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1, &
     & 1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1, &
     & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
     & 1,0,0,0,0,0]

  !> Atomic radii
  !> M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
  !> in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
  !> edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
  !> corrected Nov. 17, 2010 for the 92nd edition.
  real(wp),private,parameter :: atomicRad(1:118) = aatoau*[ &
     & 0.32_wp,0.37_wp,1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp, &
     & 0.60_wp,0.62_wp,1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp, &
     & 1.00_wp,1.01_wp,2.00_wp,1.74_wp,1.59_wp,1.48_wp,1.44_wp,1.30_wp, &
     & 1.29_wp,1.24_wp,1.18_wp,1.17_wp,1.22_wp,1.20_wp,1.23_wp,1.20_wp, &
     & 1.20_wp,1.18_wp,1.17_wp,1.16_wp,2.15_wp,1.90_wp,1.76_wp,1.64_wp, &
     & 1.56_wp,1.46_wp,1.38_wp,1.36_wp,1.34_wp,1.30_wp,1.36_wp,1.40_wp, &
     & 1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp,2.38_wp,2.06_wp, &
     & 1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp,1.82_wp, &
     & 1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp,1.74_wp,1.64_wp, &
     & 1.58_wp,1.50_wp,1.41_wp,1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp, &
     & 1.44_wp,1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp,2.42_wp,2.11_wp, &
     & 2.01_wp,1.90_wp,1.84_wp,1.83_wp,1.80_wp,1.80_wp,1.73_wp,1.68_wp, &
     & 1.68_wp,1.68_wp,1.65_wp,1.67_wp,1.73_wp,1.76_wp,1.61_wp,1.57_wp, &
     & 1.49_wp,1.43_wp,1.41_wp,1.34_wp,1.29_wp,1.28_wp,1.21_wp,1.22_wp, &
     & 1.36_wp,1.43_wp,1.62_wp,1.75_wp,1.65_wp,1.57_wp]

  public :: shakedata
  !=========================================================================================!
  !data object that contains settings for a molecular dynamics simulation.
  type :: shakedata

    logical :: initialized = .false.

    integer :: shake_mode = 0 ! 0 = off/user, 1 = H only, 2 = all

    !>--- user defined SHAKE bonds
    integer :: nusr = 0 !number of user-defined shake bonds
    integer,allocatable :: conslistu(:,:)

    !>--- topological connectivity/WBO matrix
    real(wp),allocatable :: wbo(:,:)

    !>--- final SHAKE constrainment data
    integer :: ncons = 0
    integer,allocatable :: conslist(:,:)
    real(wp),allocatable :: distcons(:)
    real(wp),allocatable :: dro(:,:)
    real(wp),allocatable :: dr(:,:)

    !>--- SHAKE convergence parameter
    integer :: maxcyc = 250
    real(wp) :: tolshake = 1.d-7

    !>--- temporary coordinate space (to avoid reallocation every shake call)
    real(wp),allocatable :: xyzt(:,:)

    !>--- pointer to the freeze list
    logical,pointer :: freezeptr(:)

  end type shakedata

  public :: init_shake
  public :: do_shake


!========================================================================================!
!========================================================================================!
contains !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine init_shake(nat,at,xyz,shk,pr)
!********************************************
!* subroutine init_shake
!* initialize SHAKE algorithm by documenting 
!* which pairs of atoms to constrain
!********************************************
    implicit none
    integer,intent(in) :: nat
    integer,intent(in) :: at(nat)
    real(wp),intent(in) :: xyz(3,nat)
    type(shakedata) :: shk
    logical,intent(in) :: pr
    integer :: iat,jat,i,j,jmin
    real(wp) :: minrij,rij,wthr
    real(wp) :: drij(3)

    integer :: nconsu
    integer :: ij
    integer,allocatable :: list(:)
    real(wp) :: rco
    logical :: metalbond,rcut,checkfreeze
    integer :: nfrz

    integer :: n2
    integer,allocatable :: cons2(:,:)

!>--- if shake was already set up, return
    if(shk%initialized) return

!>--- reset counter
    nconsu = 0
    n2 = 0
    checkfreeze = associated(shk%freezeptr)
    !> gfortran workaround
    if(checkfreeze)then !> gfortran check
      nfrz = size(shk%freezeptr,1)
    else
      nfrz = 0  
    endif

!>--- count user-defined bonds
    if ((shk%nusr > 0) .and. allocated(shk%conslistu)) then
      do i = 1,shk%nusr
        nconsu = nconsu + 1
      end do
    end if

!>--- constrain all X-H only
    if (shk%shake_mode == 1) then
      ij = nat * (nat + 1) / 2
      allocate (cons2(2,ij),source=0)
      do i = 1,nat
        if (at(i) .eq. 1) then
          minrij = 1000.0_wp
          do j = 1,nat
            if (j .ne. i) then

              !> check if BOTH are frozen
              if(checkfreeze.and.j<=nfrz.and.i<=nfrz)then
                if( shk%freezeptr(i) .and. shk%freezeptr(j) ) cycle
              endif 
              
              rij = (xyz(1,i) - xyz(1,j))**2 + (xyz(2,i) - xyz(2,j))**2 + (xyz(3,i) - xyz(3,j))**2
              if (rij .lt. minrij) then
                minrij = rij
                jmin = j
              end if
            end if
          end do
          if (at(jmin) .eq. 1) then
            if (jmin .gt. i) then
              nconsu = nconsu + 1
              n2 = n2 + 1
              cons2(1,n2) = i
              cons2(2,n2) = jmin
            end if
          else
            nconsu = nconsu + 1
            n2 = n2 + 1
            cons2(1,n2) = i
            cons2(2,n2) = jmin
          end if
        end if
      end do
    end if

!>--- SHAKE all bonds
    if (shk%shake_mode == 2) then
      if (allocated(shk%wbo)) then
        ij = nat * (nat + 1) / 2
        allocate (cons2(2,ij),source=0)
        allocate (list(ij),source=0)
        do i = 1,nat
          do j = 1,nat
            if (i .eq. j) cycle

            !>--- if both are frozen, no SHAKE!
            if(checkfreeze.and.j<=nfrz.and.i<=nfrz)then
              if( shk%freezeptr(i) .and. shk%freezeptr(j) ) cycle
            endif

            rij = (xyz(1,i) - xyz(1,j))**2 + (xyz(2,i) - xyz(2,j))**2 + (xyz(3,i) - xyz(3,j))**2
            rco = (atomicRad(at(j)) + atomicRad(at(i))) * autoaa
            ij = shake_lin(i,j)
            !>--- to consider?
            rcut = (0.52917726_wp * sqrt(rij) .lt. 1.2_wp * rco) .and. (list(ij) .eq. 0)
            !>--- metal bond?
            metalbond = metal(at(i)) .eq. 1 .or. metal(at(j)) .eq. 1
            !>--- WBO threshold. if WBO > thr, the bond is constrained
            wthr = 0.5_wp
            !>--- modify wthr, e.g. K...O in phosphates have WBO around 0.15
            if (metalbond) wthr = 0.1_wp
            !>--- check relevant
            if (rcut .and. shk%wbo(i,j) .gt. wthr) then
              !>--- do not constrain M bonds except Li/Be
              metalbond = metal(at(i)) .eq. 1 .or. metal(at(j)) .eq. 1
              if (at(i) .eq. 3 .or. at(i) .eq. 4 .or. at(j) .eq. 3 .or. at(j) .eq. 4) &
                 &         metalbond = .false.
              if (metalbond) then
                if ((i .lt. j) .and. pr) &
                  & write (*,*) 'init_shake: metal bond ',i,j,'not constrained'
              else
                list(ij) = 1
                nconsu = nconsu + 1
                n2 = n2 + 1
                cons2(1,n2) = i
                cons2(2,n2) = j
              end if
            end if
          end do
        end do
        deallocate (list)
      else
        write (*,*) 'No bonding information provided!'
        write (*,*) 'Automatic SHAKE setup failed.'
        error stop
      end if
    end if

    shk%initialized = .true.

    shk%ncons = nconsu
    if (nconsu .lt. 1) then
      if (allocated(cons2)) deallocate (cons2)
      return
    end if

    allocate (shk%conslist(2,shk%ncons),shk%distcons(shk%ncons))
    allocate (shk%dro(3,shk%ncons),shk%dr(4,shk%ncons))

    if (shk%nusr > 0) then
      shk%conslist(1:2,1:shk%nusr) = shk%conslistu(1:2,1:shk%nusr)
      j = shk%nusr + 1
    else
      j = 1
    end if
    if (n2 > 0) then
      shk%conslist(1:2,j:shk%ncons) = cons2(1:2,1:n2)
    end if
    do i = 1,shk%ncons
      iat = shk%conslist(1,i)
      jat = shk%conslist(2,i)
      drij(1:3) = xyz(1:3,iat) - xyz(1:3,jat)
      shk%distcons(i) = drij(1)**2 + drij(2)**2 + drij(3)**2
    end do

    if (allocated(cons2)) deallocate (cons2)
    return
  end subroutine init_shake

!========================================================================================!

  subroutine do_shake(nat,xyzo,xyz,velo,acc,mass,tstep,shk,pr,iostat)
!**************************************************
!* subroutine do_shake
!* Calculation of the actual SHAKE constraint
!**************************************************
    implicit none
    integer,intent(in) :: nat
    type(shakedata) :: shk
    logical,intent(in) :: pr
    integer,intent(out),optional :: iostat
    real(wp) :: xyzo(3,nat),xyz(3,nat),velo(3,nat),acc(3,nat),mass(nat)
    real(wp) :: virsh(3),tstep
    real(wp),allocatable :: xyzt(:,:)
    real(wp) :: vel(3)
    integer :: i,icyc,io
    integer :: iat,jat
    logical :: conv
    real(wp) :: maxdev
    real(wp) :: r,dev,dist,denom
    real(wp) :: gcons,rmi,rmj
    real(wp) :: tau1,tau2
    integer :: jmaxdev

    conv = .false. !> assume not converged
    io = 1         !> and therefore unsuccessfull
    icyc = 0

    if(.not.allocated(shk%xyzt)) allocate(shk%xyzt(3,nat))
    shk%xyzt = xyz

    tau1 = 1.0_wp / tstep
    tau2 = tau1 * tau1

    do i = 1,shk%ncons
      iat = shk%conslist(1,i)
      jat = shk%conslist(2,i)
      shk%dro(1:3,i) = xyzo(1:3,iat) - xyzo(1:3,jat)
    end do

!>--- iterative SHAKE loop
    do
      maxdev = 0.d0

      do i = 1,shk%ncons
        iat = shk%conslist(1,i)
        jat = shk%conslist(2,i)
        shk%dr(1:3,i) = shk%xyzt(1:3,iat) - shk%xyzt(1:3,jat)
        shk%dr(4,i) = shk%dr(1,i)**2 + shk%dr(2,i)**2 + shk%dr(3,i)**2
        dist = shk%distcons(i)
        dev = abs(shk%dr(4,i) - dist) / dist
        if (dev .gt. maxdev) then
          maxdev = dev
          jmaxdev = i
        end if
      end do

      if (maxdev .lt. shk%tolshake) conv = .true.

      if (.not. conv) then
        do i = 1,shk%ncons
          iat = shk%conslist(1,i)
          jat = shk%conslist(2,i)
          dist = shk%distcons(i)
          rmi = 1.0_wp / mass(iat)
          rmj = 1.0_wp / mass(jat)
          denom = 2.0_wp * (rmi + rmj) * (shk%dr(1,i) * shk%dro(1,i) +  &
             &     shk%dr(2,i) * shk%dro(2,i) + &
             &                            shk%dr(3,i) * shk%dro(3,i))
          gcons = (dist - shk%dr(4,i)) / denom
          shk%xyzt(1:3,iat) = shk%xyzt(1:3,iat) + rmi * gcons * shk%dro(1:3,i)
          shk%xyzt(1:3,jat) = shk%xyzt(1:3,jat) - rmj * gcons * shk%dro(1:3,i)
        end do
      end if
      icyc = icyc + 1
      if (.not. conv .and. icyc .le. shk%maxcyc) cycle
      exit
    end do

    if (conv) then
      velo = velo + (shk%xyzt - xyz) * tau1
      acc = acc + (shk%xyzt - xyz) * tau2
      xyz = shk%xyzt
      io = 0 !> successful termination
    else if (pr) then
      write (*,*) 'SHAKE did not converge! maxdev=',maxdev
    end if

    if(present(iostat))then
      iostat = io
    endif 

    return
  end subroutine do_shake

!========================================================================================!
  integer function shake_lin(i1,i2)
!**********************
!* helper function lin
!**********************
    integer :: i1,i2,idum1,idum2
    idum1 = max(i1,i2)
    idum2 = min(i1,i2)
    shake_lin = idum2 + idum1 * (idum1 - 1) / 2
    return
  end function shake_lin

!========================================================================================!
!========================================================================================!
end module shake_module
