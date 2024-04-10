module crest_parameters
  use iso_fortran_env, only: wp => real64, sp => real32
  use iso_fortran_env, only: ap => real64 
  use iso_fortran_env, only: dp => int64
  use iso_fortran_env, only: int8,int16,int32,int64,real64,real32
  use iso_fortran_env, only: stdout => output_unit
  use iso_fortran_env, only: stderr => error_unit

  public :: wp,sp,ap,dp,stdout,stderr
  public :: int8,int16,int32,int64,real64,real32

  real(wp),parameter,public :: bohr = 0.52917726_wp
  real(wp),parameter,public :: angstrom = 1.0_wp / bohr
  real(wp),parameter,public :: autoaa = bohr
  real(wp),parameter,public :: aatoau = angstrom

  real(wp),parameter,public :: pi  = acos(0.0_wp)*2.0_wp
  real(wp),parameter,public :: radtodeg = 180.0_wp / pi
  real(wp),parameter,public :: degtorad = 1.0_wp / radtodeg

  real(wp),parameter,public :: amutokg = 1.660539040e-27_wp
  real(wp),parameter,public :: autokj = 2625.49964038_wp
  real(wp),parameter,public :: autokcal = 627.50947428_wp
  real(wp),parameter,public :: autoeV = 27.211324570273_wp
  real(wp),parameter,public :: kcaltoau = 1.0_wp/autokcal
  real(wp),parameter,public :: kcaltokj = autokj/autokcal
  real(wp),parameter,public :: autorcm = 219474.63067_wp
  real(wp),parameter,public :: rcmtoau = 1.0_wp/autorcm
  real(wp),parameter,public :: metokg  = 9.10938356e-31_wp
  real(wp),parameter,public :: kgtome  = 1.0_wp/metokg

  real(wp),parameter,public :: Rcal = 8.31446261815324_wp/kcaltokj
  real(wp),parameter,public :: kB = 3.166808578545117e-06_wp
  real(wp),parameter,public :: avogadro = 6.0221413e23_wp ! 1/mol
  real(wp),parameter,public :: planck = 6.62606957e-34_wp ! J*s
  real(wp),parameter,public :: hbar = planck/(2.0_wp*pi)

  real(wp),public,parameter :: lightspeed = 137.0359990740_wp 
  !> femtosectons to atomic time units
  real(wp), public, parameter :: fstoau = 41.3413733365614_wp
  !> Coulomb to atomic charge units (electrons)
  real(wp), public, parameter :: autoc = 1.6021766208e-19_wp
  !> Debye to atomic units
  real(wp), public, parameter :: autod = autoc * lightspeed * autoaa**2 * fstoau * 1.0e+16_wp
  real(wp), public, parameter :: dtoau = 1.0_wp / autod

  character(len=1),public,parameter :: sep = '/'
  character(len=12),public,parameter :: dev0 = ' 2>/dev/null'

end module crest_parameters
