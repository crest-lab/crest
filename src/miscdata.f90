!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2018-2020 Stefan Grimme
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

module miscdata
  use crest_parameters
  implicit none
  public

!&<
  !> Element symbols
  character(len=2),parameter :: PSE(118) = [ &
   & 'H ',                                                                                'He', &
   & 'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &
   & 'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &
   & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
   & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
   & 'Cs','Ba','La',                                                                            &
   &                'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',      &
   &                'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
   & 'Fr','Ra','Ac',                                                                            &
   &                'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',      &
   &                'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' ]
!&>

!&<
  !> Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
  !> 188-197), values for metals decreased by 10%.
  !> As used in D3
  real(wp), parameter :: rcov(1:118) = [ &
   & 0.32_wp,0.46_wp, & ! H,He
   & 1.20_wp,0.94_wp,0.77_wp,0.75_wp,0.71_wp,0.63_wp,0.64_wp,0.67_wp, & ! Li-Ne
   & 1.40_wp,1.25_wp,1.13_wp,1.04_wp,1.10_wp,1.02_wp,0.99_wp,0.96_wp, & ! Na-Ar
   & 1.76_wp,1.54_wp, & ! K,Ca
   &                 1.33_wp,1.22_wp,1.21_wp,1.10_wp,1.07_wp, & ! Sc-
   &                 1.04_wp,1.00_wp,0.99_wp,1.01_wp,1.09_wp, & ! -Zn
   &                 1.12_wp,1.09_wp,1.15_wp,1.10_wp,1.14_wp,1.17_wp, & ! Ga-Kr
   & 1.89_wp,1.67_wp, & ! Rb,Sr
   &                 1.47_wp,1.39_wp,1.32_wp,1.24_wp,1.15_wp, & ! Y-
   &                 1.13_wp,1.13_wp,1.08_wp,1.15_wp,1.23_wp, & ! -Cd
   &                 1.28_wp,1.26_wp,1.26_wp,1.23_wp,1.32_wp,1.31_wp, & ! In-Xe
   & 2.09_wp,1.76_wp, & ! Cs,Ba
   &         1.62_wp,1.47_wp,1.58_wp,1.57_wp,1.56_wp,1.55_wp,1.51_wp, & ! La-Eu
   &         1.52_wp,1.51_wp,1.50_wp,1.49_wp,1.49_wp,1.48_wp,1.53_wp, & ! Gd-Yb
   &                 1.46_wp,1.37_wp,1.31_wp,1.23_wp,1.18_wp, & ! Lu-
   &                 1.16_wp,1.11_wp,1.12_wp,1.13_wp,1.32_wp, & ! -Hg
   &                 1.30_wp,1.30_wp,1.36_wp,1.31_wp,1.38_wp,1.42_wp, & ! Tl-Rn
   & 2.01_wp,1.81_wp, & ! Fr,Ra
   &         1.67_wp,1.58_wp,1.52_wp,1.53_wp,1.54_wp,1.55_wp,1.49_wp, & ! Ac-Am
   &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
   &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
   &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
   &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp] & ! Nh-Og
   & * aatoau * 4.0_wp / 3.0_wp
!&>

!&<
   ! Radius used in QCxMS (in au)
   real(wp), parameter :: QCxMS_Rad(118) = aatoau *  [ &
   & 0.32_wp,0.37_wp, & ! H,He
   & 1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp,0.60_wp,0.62_wp, & ! Li-Ne
   & 1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp,1.00_wp,1.01_wp, & ! Na-Ar
   & 2.00_wp,1.74_wp, & ! K,Ca
   &                 1.59_wp,1.48_wp,1.44_wp,1.30_wp,1.29_wp, & ! Sc-
   &                 1.24_wp,1.18_wp,1.17_wp,1.22_wp,1.20_wp, & ! -Zn
   &                 1.23_wp,1.20_wp,1.20_wp,1.18_wp,1.17_wp,1.16_wp, & ! Ga-Kr
   & 2.15_wp,1.90_wp, & ! Rb,Sr
   &                 1.76_wp,1.64_wp,1.56_wp,1.46_wp,1.38_wp, & ! Y-
   &                 1.36_wp,1.34_wp,1.30_wp,1.36_wp,1.40_wp, & ! -Cd
   &                 1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp, & ! In-Xe
   & 2.38_wp,2.06_wp, & ! Cs,Ba
   &         1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp, & ! La-Eu
   &         1.82_wp,1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp, & ! Gd-Yb
   &                 1.74_wp,1.64_wp,1.58_wp,1.50_wp,1.41_wp, & ! Lu-
   &                 1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp, & ! -Hg
   &                 1.44_wp,1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp, & ! Tl-Rn
   & 2.42_wp,2.11_wp, & ! Fr,Ra
   &         2.01_wp,1.90_wp,1.84_wp,1.83_wp,1.80_wp,1.80_wp,& ! Ac-Pu
   ! from covalent 2009 covalent radii, such that it is complete up to 118
   &                                                         1.49_wp, & ! Am
   &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
   &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
   &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
   &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og   
!&>


!&<
  !> D3 pairwise van-der-Waals radii (only homoatomic pairs present here)
  real(wp), parameter :: vdw_D3(1:94) = aatoau * [&
   & 1.09155_wp, 0.86735_wp, 1.74780_wp, 1.54910_wp, &
   & 1.60800_wp, 1.45515_wp, 1.31125_wp, 1.24085_wp, &
   & 1.14980_wp, 1.06870_wp, 1.85410_wp, 1.74195_wp, &
   & 2.00530_wp, 1.89585_wp, 1.75085_wp, 1.65535_wp, &
   & 1.55230_wp, 1.45740_wp, 2.12055_wp, 2.05175_wp, &
   & 1.94515_wp, 1.88210_wp, 1.86055_wp, 1.72070_wp, &
   & 1.77310_wp, 1.72105_wp, 1.71635_wp, 1.67310_wp, &
   & 1.65040_wp, 1.61545_wp, 1.97895_wp, 1.93095_wp, &
   & 1.83125_wp, 1.76340_wp, 1.68310_wp, 1.60480_wp, &
   & 2.30880_wp, 2.23820_wp, 2.10980_wp, 2.02985_wp, &
   & 1.92980_wp, 1.87715_wp, 1.78450_wp, 1.73115_wp, &
   & 1.69875_wp, 1.67625_wp, 1.66540_wp, 1.73100_wp, &
   & 2.13115_wp, 2.09370_wp, 2.00750_wp, 1.94505_wp, &
   & 1.86900_wp, 1.79445_wp, 2.52835_wp, 2.59070_wp, &
   & 2.31305_wp, 2.31005_wp, 2.28510_wp, 2.26355_wp, &
   & 2.24480_wp, 2.22575_wp, 2.21170_wp, 2.06215_wp, &
   & 2.12135_wp, 2.07705_wp, 2.13970_wp, 2.12250_wp, &
   & 2.11040_wp, 2.09930_wp, 2.00650_wp, 2.12250_wp, &
   & 2.04900_wp, 1.99275_wp, 1.94775_wp, 1.87450_wp, &
   & 1.72280_wp, 1.67625_wp, 1.62820_wp, 1.67995_wp, &
   & 2.15635_wp, 2.13820_wp, 2.05875_wp, 2.00270_wp, &
   & 1.93220_wp, 1.86080_wp, 2.53980_wp, 2.46470_wp, &
   & 2.35215_wp, 2.21260_wp, 2.22970_wp, 2.19785_wp, &
   & 2.17695_wp, 2.21705_wp]
!&>

end module miscdata
