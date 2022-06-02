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

module atmasses
      use iso_fortran_env, only : wp => real64
      implicit none
      private :: wp
      public

      real(wp) ::  ams(107)
      data ams       /1.00790_wp,  4.00260_wp,  6.94000_wp,  9.01218_wp, &
      &  10.81000_wp, 12.01100_wp, 14.00670_wp, 15.99940_wp, 18.99840_wp, &
      &  20.17900_wp, 22.98977_wp, 24.30500_wp, 26.98154_wp, 28.08550_wp, &
      &  30.97376_wp, 32.06000_wp, 35.45300_wp, 39.94800_wp, 39.09830_wp, &
      &  40.08000_wp, 44.95590_wp, 47.90000_wp, 50.94150_wp, 51.99600_wp, &
      &  54.93800_wp, 55.84700_wp, 58.93320_wp, 58.71000_wp, 63.54600_wp, &
      &  65.38000_wp, 69.73500_wp, 72.59000_wp, 74.92160_wp, 78.96000_wp, &
      &  79.90400_wp, 83.80000_wp, 85.46780_wp, 87.62000_wp, 88.90590_wp, &
      &  91.22000_wp, 92.90640_wp, 95.94000_wp, 98.90620_wp, 101.0700_wp, &
      &  102.9055_wp, 106.4000_wp, 107.8680_wp, 112.4100_wp, 114.8200_wp, &
      &  118.6900_wp, 121.7500_wp, 127.6000_wp, 126.9045_wp, 131.3000_wp, &
      &  132.9054_wp, 137.3300_wp,  &
      &  138.91,140.12,140.91,144.24,147.00,150.36,151.97,157.25, &
      &  158.93,162.50,164.93,167.26,168.93,173.04,174.97, &
      &  178.4900_wp, 180.9479_wp, &
      &  183.8500_wp, 186.2070_wp, 190.2000_wp, 192.2200_wp, 195.0900_wp, &
      &  196.9665_wp, 200.5900_wp, 204.3700_wp, 207.2000_wp, 208.9804_wp, &
      &  209.,210.,222.,21*0.000_wp/

contains

!------------------------------------------------------------
! This is a subroutine that replaces masses with fitted or scaled
!  masses to increase the quality of frequencies in IR spectra calculations
!
subroutine freqmass(amss,param)
      implicit none
      real(wp),intent(inout) :: amss(107) 
      character(len=*) :: param
 
      select case( param )
       case( 'gfn2' )
          amss=amss*1.016_wp

!          open(newunit=ich,file='~/.mscal')
!          read(ich,*) scal
!          close(ich)

          amss( 1)=  1.032_wp
          amss( 6)= 11.79_wp
          amss( 7)= 14.82_wp
          amss( 8)= 17.11_wp
          amss( 9)= 16.17_wp
          !amss(14)= ams(14)*scal
          !amss(14)= ams(14)*1.163_wp
          amss(14)= 32.6634_wp
          !amss(15)= ams(15)*scal
          !amss(15)= ams(15)*0.847_wp
          amss(15)= 26.2348_wp
          amss(16)= 22.39_wp
          amss(17)= 30.23_wp
          amss(35)= 55.34_wp
       case( 'b3lyp-3c' )
          amss=amss*1.07_wp

          amss( 1)=  1.072_wp
          amss( 6)= 12.60_wp
          amss( 7)= 15.48_wp
          amss( 8)= 17.75_wp
          amss( 9)= 16.92_wp
          amss(16)= 24.99_wp
          amss(17)= 30.55_wp
          amss(35)= 55.83_wp
       case default
        continue
      end select

      return
end subroutine freqmass


!=========================================
! calculate total weight of molecule
!=========================================
function molweight(nat,at)
      implicit none
      real(wp) :: molweight
      integer :: nat
      integer :: at(nat)
      integer :: i
      molweight = 0.0_wp
      do i=1,nat
       molweight = molweight + ams(at(i))
      enddo
      return
end function molweight


end module atmasses
