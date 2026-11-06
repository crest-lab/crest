!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2022-2026 Christoph Plett, Philipp Pracht
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
! Ported from the xtb docking submodule (github.com/grimme-lab/xtb,
! src/docking/param.f90, src/sphereparam.f90), LGPL-3.0-or-later.
!================================================================================!

!> Parameters, run state and geometry helpers of the aISS docking (xTB-IFF)
module aiss_param
  use crest_parameters,only:wp,sp,stdout,kB
#ifdef WITH_TBLITE
  use dftd4_data,only:get_covalent_rad,get_r4r2_val
#endif
  implicit none
  private

  !> xtb-internal conversion factors, kept for bitwise-comparable screening
  real(wp),parameter,public :: au = 627.509541_wp
  real(wp),parameter,public :: autoang = 0.52917726_wp
  real(wp),parameter,public :: pi = 3.14159265358979_wp
  real(wp),parameter,public :: pi2 = 2.0_wp*pi

  !> probe atom (Kr) of the rare-gas grid
  integer,parameter,public :: probe_atom_type = 36

  !> directed docking potentials
  integer,parameter,public :: p_wall_pot = 1
  integer,parameter,public :: p_atom_pot = 2
  integer,parameter,public :: p_atom_att = 3
  integer,parameter,public :: p_atom_qcg = 4
  real(wp),parameter,public :: pot_scal = 1.0_wp
  real(wp),parameter,public :: pot_expo = 2.0_wp
  real(wp),parameter,public :: attractive_pot = -0.032_wp

  real(wp),public :: par_rep_scal,par_d3_a1,par_d3_a2,par_d3_s8,par_chrg_lp
  real(wp),public :: par_chrg_sig,par_chrg_pi,par_pos1_lp,par_pos2_lp,par_pos_sig
  real(wp),public :: par_es_damp,par_drude_fc,par_drude_damp,par_xh1,par_xh2
  real(wp),public :: r0ab6(94,94),r0ab8(94,94),rrab(94,94),r0scal(94)
  real(wp),public :: val_e(86),r0_atom(86),rcov(94),r2r4(94)
  integer,public  :: lpatom(86),sigatom(86)

  !> Drude force constant scaling (chemical hardness)
  real(wp),public :: gam(94)
  data gam/&
  &0.47259288,0.92203391,0.17452888,0.25700733,0.33949086,0.42195412,&
  &0.50438193,0.58691863,0.66931351,0.75191607,0.17964105,0.22157276,&
  &0.26348578,0.30539645,0.34734014,0.38924725,0.43115670,0.47308269,&
  &0.17105469,0.20276244,0.21007322,0.21739647,0.22471039,0.23201501,&
  &0.23933969,0.24665638,0.25398255,0.26128863,0.26859476,0.27592565,&
  &0.30762999,0.33931580,0.37235985,0.40273549,0.43445776,0.46611708,&
  &0.15585079,0.18649324,0.19356210,0.20063311,0.20770522,0.21477254,&
  &0.22184614,0.22891872,0.23598621,0.24305612,0.25013018,0.25719937,&
  &0.28784780,0.31848673,0.34912431,0.37976593,0.41040808,0.44105777,&
  &0.05019332,0.06762570,0.08504445,0.10247736,0.11991105,0.13732772,&
  &0.15476297,0.17218265,0.18961288,0.20704760,0.22446752,0.24189645,&
  &0.25932503,0.27676094,0.29418231,0.31159587,0.32902274,0.34592298,&
  &0.36388048,0.38130586,0.39877476,0.41614298,0.43364510,0.45104014,&
  &0.46848986,0.48584550,0.12526730,0.14268677,0.16011615,0.17755889,&
  &0.19497557,0.21240778,0.07263525,0.09422158,0.09920295,0.10418621,&
  &0.14235633,0.16394294,0.18551941,0.22370139/

  !> atomic masses used for all CMA and inertia evaluations of the docking
  real(wp),parameter :: ams(107) = [ &
  & 1.00790_wp,4.00260_wp,6.94000_wp,9.01218_wp, &
  & 10.81000_wp,12.01100_wp,14.00670_wp,15.99940_wp,18.99840_wp, &
  & 20.17900_wp,22.98977_wp,24.30500_wp,26.98154_wp,28.08550_wp, &
  & 30.97376_wp,32.06000_wp,35.45300_wp,39.94800_wp,39.09830_wp, &
  & 40.08000_wp,44.95590_wp,47.90000_wp,50.94150_wp,51.99600_wp, &
  & 54.93800_wp,55.84700_wp,58.93320_wp,58.71000_wp,63.54600_wp, &
  & 65.38000_wp,69.73500_wp,72.59000_wp,74.92160_wp,78.96000_wp, &
  & 79.90400_wp,83.80000_wp,85.46780_wp,87.62000_wp,88.90590_wp, &
  & 91.22000_wp,92.90640_wp,95.94000_wp,98.90620_wp,101.0700_wp, &
  & 102.9055_wp,106.4000_wp,107.8680_wp,112.4100_wp,114.8200_wp, &
  & 118.6900_wp,121.7500_wp,127.6000_wp,126.9045_wp,131.3000_wp, &
  & 132.9054_wp,137.3300_wp, &
  & 138.91_wp,140.12_wp,140.91_wp,144.24_wp,147.00_wp,150.36_wp,151.97_wp,157.25_wp, &
  & 158.93_wp,162.50_wp,164.93_wp,167.26_wp,168.93_wp,173.04_wp,174.97_wp, &
  & 178.4900_wp,180.9479_wp, &
  & 183.8500_wp,186.2070_wp,190.2000_wp,192.2200_wp,195.0900_wp, &
  & 196.9665_wp,200.5900_wp,204.3700_wp,207.2000_wp,208.9804_wp, &
  & 209.0_wp,210.0_wp,222.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
  & 0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp,0.0_wp, &
  & 0.0_wp,0.0_wp,0.0_wp]

!========================================================================================!

  !> ellipsoidal polynomial wall, centered at the origin
  type,public :: aiss_wall
    real(wp) :: radius(3) = 0.0_wp
    real(wp) :: center(3) = 0.0_wp
    integer  :: alpha = 30  !> xtb sphere_alpha, V = (|r/radius|)^alpha
    !> logfermi wall, V = kT*log(1+exp(beta*(|r|-R0))) (xtb sphere_temp/beta)
    logical  :: logfermi = .false.
    real(wp) :: temp = 300.0_wp
    real(wp) :: beta = 6.0_wp
    integer,allocatable :: list(:)
  end type aiss_wall

  !> run settings (xtb dock command line and $dock block)
  type,public :: aiss_config
    integer  :: maxgen = 10
    integer  :: maxparent = 100
    real(wp) :: stepr = 2.5_wp
    real(wp) :: stepa = 45.0_wp
    integer  :: mxcma = 1000
    integer  :: mxcent_clust = 500
    integer  :: n_opt = 15
    integer  :: nfrag1 = 0
    logical  :: cssym = .false.
    logical  :: fulle = .false.
    logical  :: debug = .false.
    logical  :: pocket_grid = .false.
    logical  :: angular_grid = .true.
    logical  :: stack_grid = .true.
    logical  :: docking_ens = .false.
    logical  :: qcg = .false.
    logical  :: samerand = .false.
    logical  :: wrfiles = .true.   !> write best.xyz/final_structures.xyz
    !> 'gfnff','gfn2','gfn1' final optimization level
    character(len=:),allocatable :: optlvl
    !> 'gfn2','gfn1','gfn0' wavefunction for the LMO/IFF precomputation
    character(len=:),allocatable :: lmolvl
    real(wp) :: etemp = 300.0_wp
    integer  :: chrg(2) = 0
    integer  :: uhf(2) = 0
    !> final optimizations: ANCopt level and optional implicit solvation
    integer  :: optlev = 0
    character(len=:),allocatable :: solvmodel,solvent
  contains
    procedure :: set_fast => aiss_config_set_fast
    procedure :: set_qcg => aiss_config_set_qcg
  end type aiss_config

  !> standalone docking input (CLI --dock / TOML [dock])
  type,public :: aiss_input
    type(aiss_config) :: cfg
    character(len=:),allocatable :: fileA,fileB,xcontrol
    logical :: chrg_set(2) = .false.
    logical :: uhf_set(2) = .false.
    logical :: nfinal_set = .false.
    character(len=:),allocatable :: fix       !> atom list of frozen atoms
    character(len=:),allocatable :: directed  !> atom list for directed docking
    integer :: directed_type = p_atom_att
    real(wp) :: directed_fc = 1.0_wp
  end type aiss_input

  !> IFF state of one docking run (xtb: TIFFData plus docking globals)
  type,public :: iff_data
    integer :: n = 0,n1 = 0,n2 = 0
    integer :: nlmo1 = 0,nlmo2 = 0
    integer,allocatable  :: at1(:),at2(:),at(:)
    integer,allocatable  :: lmo1(:),lmo2(:)
    integer,allocatable  :: neigh(:,:)
    real(wp),allocatable :: xyz1(:,:),xyz2(:,:),xyz(:,:)
    real(wp),allocatable :: rlmo1(:,:),rlmo2(:,:)
    real(wp),allocatable :: q1(:),q2(:),q(:)
    real(wp),allocatable :: qdr1(:),qdr2(:)
    real(wp),allocatable :: cn1(:),cn2(:)
    real(wp),allocatable :: z1(:),z2(:)
    real(wp),allocatable :: alp1(:),alp2(:),alpab(:,:)
    real(wp),allocatable :: c6ab(:,:),cprob(:)
    real(wp),allocatable :: gab1(:,:),gab2(:,:)
    real(wp),allocatable :: den1(:,:,:),den2(:,:,:)
    real(wp),allocatable :: qct1(:,:),qct2(:,:)
    !> per-fragment electronic properties from the precomputation
    real(wp) :: dipol(2) = 0.0_wp
    real(wp) :: ehomo(2) = 0.0_wp
    real(wp) :: elumo(2) = 0.0_wp
    real(wp) :: chrg(2) = 0.0_wp
    real(wp) :: etot(2) = 0.0_wp
    real(wp) :: shift_geo = 0.0_wp
    !> copies of the settings that iff_e depends on
    logical :: fulle = .false.
    logical :: debug = .false.
    logical :: qcg = .false.
    !> directed docking
    integer :: directed_type = p_atom_att
    integer :: ndirected = 0
    integer,allocatable  :: directed_atoms(:)
    real(wp) :: directed_fc = 1.0_wp
    real(wp),allocatable :: directed_val(:)
    real(wp),allocatable :: directed_expo(:)
    !> wall potentials entering the screening energy
    integer :: nwalls = 0
    type(aiss_wall),allocatable :: walls(:)
  contains
    procedure :: allocate => iff_data_allocate
  end type iff_data

  public :: set_iff_param
  public :: diptot,axis_docking,cmadock,rcma
  public :: rotmat,move,move2,rnorm
  public :: rand6,irand,crossover,sort6
  public :: mrec
  public :: wrc2,wrc3,wr_grid
  public :: rsp_packed
  public :: cavity_egrad,cavitye

  interface
    recursive subroutine qsort(a,first,last,ind)
      import :: wp
      real(wp) :: a(*)
      integer :: ind(*)
      integer :: first,last
    end subroutine qsort
  end interface
  public :: qsort

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine aiss_config_set_fast(self)
    !***********************************
    !* xtb dock --fast
    !***********************************
    class(aiss_config),intent(inout) :: self
    self%maxparent = 30
    self%maxgen = 7
    self%mxcma = 250
    self%stepr = 4.0_wp
    self%stepa = 60.0_wp
    self%n_opt = 4
  end subroutine aiss_config_set_fast

  subroutine aiss_config_set_qcg(self)
    !***********************************
    !* xtb dock --qcg
    !***********************************
    class(aiss_config),intent(inout) :: self
    self%maxparent = 50
    self%maxgen = 7
    self%n_opt = 5
    self%qcg = .true.
  end subroutine aiss_config_set_qcg

!========================================================================================!

  subroutine iff_data_allocate(self,n1,n2)
    class(iff_data),intent(inout) :: self
    integer,intent(in) :: n1,n2
    integer :: n
    self%n1 = n1
    self%n2 = n2
    self%n = n1+n2
    n = self%n
    allocate (self%at1(n1),self%at2(n2),source=0)
    allocate (self%xyz1(3,n1),self%q1(n1),self%cn1(n1),self%alp1(n1), &
    &         self%qdr1(n1),self%z1(n1),self%den1(2,4,n1),self%gab1(n1,n1),source=0.0_wp)
    allocate (self%xyz2(3,n2),self%q2(n2),self%cn2(n2),self%alp2(n2), &
    &         self%qdr2(n2),self%z2(n2),self%den2(2,4,n2),self%gab2(n2,n2),source=0.0_wp)
    allocate (self%rlmo1(4,10*n1),self%qct1(n1,2),self%rlmo2(4,10*n2),self%qct2(n2,2), &
    &         source=0.0_wp)
    allocate (self%lmo1(10*n1),self%lmo2(10*n2),source=0)
    allocate (self%at(n),source=0)
    allocate (self%neigh(0:n,n),source=0)
    allocate (self%xyz(3,n),self%q(n),self%c6ab(n,n),self%alpab(n2,n1),self%cprob(n1), &
    &         source=0.0_wp)
  end subroutine iff_data_allocate

!========================================================================================!

  subroutine set_iff_param(iunit)
    !*********************************************************
    !* Global IFF parameters. An xtb-style '.param.xtbiff'
    !* in the working directory overrides the defaults.
    !*********************************************************
    integer,intent(in) :: iunit
    character(len=*),parameter :: fname = '.param.xtbiff'
    logical :: ex
    integer :: i,j,ich

#ifdef WITH_TBLITE
    do i = 1,94
      rcov(i) = get_covalent_rad(i)
      r2r4(i) = get_r4r2_val(i)
    end do
#endif
    call valel(val_e)

    inquire (file=fname,exist=ex)
    if (ex) then
      if (iunit /= 0) write (iunit,*) 'reading parameter file ',fname
      open (newunit=ich,file=fname)
      read (ich,*) par_rep_scal
      read (ich,*) par_d3_a1
      read (ich,*) par_d3_a2
      read (ich,*) par_d3_s8
      read (ich,*) par_chrg_lp
      read (ich,*) par_chrg_sig
      read (ich,*) par_chrg_pi
      read (ich,*) par_pos1_lp
      read (ich,*) par_es_damp
      read (ich,*) par_xh1
      read (ich,*) par_xh2
      read (ich,*) par_drude_fc
      do i = 1,86
        read (ich,*) j,r0scal(i)
      end do
      close (ich)
    else
      if (iunit /= 0) write (iunit,*) 'taking internal default parameters'
      par_rep_scal = 0.1232398_wp
      par_d3_a1 = 0.4459020_wp
      par_d3_a2 = 4.4216615_wp
      par_d3_s8 = 2.0000000_wp
      par_chrg_lp = 3.1064458_wp
      par_chrg_sig = 2.7000000_wp
      par_chrg_pi = 0.8892689_wp
      par_pos1_lp = 3.7082486_wp
      par_es_damp = 0.2298506_wp
      par_xh1 = 0.3614582_wp
      par_xh2 = 6.7354042_wp
      par_drude_fc = 0.9811899_wp
      r0scal(1:86) = [ &
      & 0.65_wp,0.65_wp,0.85_wp,1.10_wp,1.20_wp,1.2606324_wp,1.37447_wp,1.1377692_wp, &
      & 0.93_wp,0.85_wp,0.85_wp,1.30_wp,1.35_wp,1.35_wp,1.60_wp,1.60_wp,1.60_wp,0.90_wp, &
      & 0.85_wp,1.40_wp,1.35_wp,1.35_wp,1.35_wp,1.35_wp,1.35_wp,1.35_wp,1.60_wp,1.60_wp, &
      & 1.60_wp,1.60_wp,1.55_wp,1.55_wp,1.55_wp,1.60_wp,1.60_wp,0.90_wp,0.95_wp,1.50_wp, &
      & 1.35_wp,1.35_wp,1.60_wp,1.35_wp,1.65_wp,1.35_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp, &
      & 1.35_wp,1.35_wp,1.55_wp,1.55_wp,1.60_wp,0.95_wp,0.95_wp,1.50_wp,1.35_wp,1.60_wp, &
      & 1.60_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp, &
      & 1.60_wp,1.60_wp,1.60_wp,1.35_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp,1.60_wp, &
      & 1.60_wp,1.60_wp,1.45_wp,1.45_wp,1.45_wp,1.60_wp,1.60_wp,0.95_wp]
    end if

    !> hand-made, not fitted
    par_drude_damp = 0.020_wp
    par_pos_sig = par_pos1_lp
    par_pos2_lp = 0.0_wp

    !> atoms with LP electrostatics
    lpatom = 0
    lpatom(6:9) = 1
    lpatom(14:17) = 1
    lpatom(32:35) = 1
    lpatom(50:53) = 1
    lpatom(82:85) = 1

    !> atoms with sigma-hole quadrupole electrostatics
    sigatom = 0
    sigatom(9) = 1
    sigatom(17) = 1
    sigatom(35) = 1
    sigatom(53) = 1
    sigatom(85) = 1
  end subroutine set_iff_param

!========================================================================================!

  subroutine valel(z)
    !***********************************
    !* valence electrons per element
    !***********************************
    real(wp),intent(out) :: z(86)
    z(1:2) = [1,2]
    z(3:10) = [1,2,3,4,5,6,7,8]
    z(11:18) = [1,2,3,4,5,6,7,8]
    z(19:36) = [1,2,3,4,5,6,7,8,9,10,11,2,3,4,5,6,7,8]
    z(37:54) = [1,2,3,4,5,6,7,8,9,10,11,2,3,4,5,6,7,8]
    z(55:56) = [1,2]
    z(57:81) = 3
    z(82:86) = [4,5,6,7,8]
  end subroutine valel

!========================================================================================!

  subroutine diptot(n1,n2,nl1,nl2,c1,c2,cl1,cl2,q1,q2,dipol)
    !***********************************************
    !* dipole of the complex from atomic and
    !* off-center LMO charges
    !***********************************************
    integer,intent(in)   :: n1,n2,nl1,nl2
    real(wp),intent(in)  :: c1(3,n1),c2(3,n2)
    real(wp),intent(in)  :: cl1(4,n1*10),cl2(4,n2*10)
    real(wp),intent(in)  :: q1(n1),q2(n2)
    real(wp),intent(out) :: dipol(3)
    integer :: i
    dipol = 0
    do i = 1,nl1
      dipol(1:3) = dipol(1:3)+cl1(1:3,i)*cl1(4,i)
    end do
    do i = 1,n1
      dipol(1:3) = dipol(1:3)+c1(1:3,i)*q1(i)
    end do
    do i = 1,nl2
      dipol(1:3) = dipol(1:3)+cl2(1:3,i)*cl2(4,i)
    end do
    do i = 1,n2
      dipol(1:3) = dipol(1:3)+c2(1:3,i)*q2(i)
    end do
  end subroutine diptot

!========================================================================================!

  subroutine axis_docking(numat,nat,coord,sum3,sumw,eig,evec)
    !***********************************************
    !* CMA (sum3), total mass (sumw), principal
    !* moments per mass (eig) and axes (evec)
    !***********************************************
    integer,intent(in)   :: numat,nat(numat)
    real(wp),intent(in)  :: coord(3,*)
    real(wp),intent(out) :: sum3(3),sumw,eig(3),evec(3,3)
    real(wp) :: t(6),x(numat),y(numat),z(numat)
    real(wp) :: sumwx,sumwy,sumwz,m
    integer :: i

    sumw = 0
    sumwx = 0.0_wp
    sumwy = 0.0_wp
    sumwz = 0.0_wp
    do i = 1,numat
      m = ams(nat(i))
      sumw = sumw+m
      sumwx = sumwx+m*coord(1,i)
      sumwy = sumwy+m*coord(2,i)
      sumwz = sumwz+m*coord(3,i)
    end do
    sum3(1) = sumwx/sumw
    sum3(2) = sumwy/sumw
    sum3(3) = sumwz/sumw

    do i = 1,numat
      x(i) = coord(1,i)-sum3(1)
      y(i) = coord(2,i)-sum3(2)
      z(i) = coord(3,i)-sum3(3)
    end do
    do i = 1,6
      t(i) = real(i,wp)*1.0e-10_wp
    end do
    do i = 1,numat
      m = ams(nat(i))
      t(1) = t(1)+m*(y(i)**2+z(i)**2)
      t(2) = t(2)-m*x(i)*y(i)
      t(3) = t(3)+m*(z(i)**2+x(i)**2)
      t(4) = t(4)-m*z(i)*x(i)
      t(5) = t(5)-m*y(i)*z(i)
      t(6) = t(6)+m*(x(i)**2+y(i)**2)
    end do
    call rsp_packed(t,3,eig,evec)
    eig = eig/sumw
  end subroutine axis_docking

!========================================================================================!

  subroutine cmadock(n,numat,nat,coord,sum3)
    !***********************************************
    !* center of mass of the first n atoms
    !***********************************************
    integer,intent(in)   :: n,numat,nat(numat)
    real(wp),intent(in)  :: coord(3,numat)
    real(wp),intent(out) :: sum3(3)
    real(wp) :: sumw,sumwx,sumwy,sumwz
    integer :: i
    sumw = 1.e-20_wp
    sumwx = 0.0_wp
    sumwy = 0.0_wp
    sumwz = 0.0_wp
    do i = 1,n
      sumw = sumw+ams(nat(i))
      sumwx = sumwx+ams(nat(i))*coord(1,i)
      sumwy = sumwy+ams(nat(i))*coord(2,i)
      sumwz = sumwz+ams(nat(i))*coord(3,i)
    end do
    sum3(1) = sumwx/sumw
    sum3(2) = sumwy/sumw
    sum3(3) = sumwz/sumw
  end subroutine cmadock

!========================================================================================!

  subroutine rcma(n1,xyz1,iz1,n2,xyz2,iz2,r,rmin)
    !***********************************************
    !* CMA distance and closest contact of A and B
    !***********************************************
    integer,intent(in)   :: n1,n2
    real(wp),intent(in)  :: xyz1(3,n1),xyz2(3,n2)
    integer,intent(in)   :: iz1(n1),iz2(n2)
    real(wp),intent(out) :: r,rmin
    real(wp) :: x1(3),x2(3),rr
    integer :: i,j
    call cmadock(n1,n1,iz1,xyz1,x1)
    call cmadock(n2,n2,iz2,xyz2,x2)
    r = sqrt((x1(1)-x2(1))**2+(x1(2)-x2(2))**2+(x1(3)-x2(3))**2)
    rmin = 1.d+42
    do i = 1,n1
      do j = 1,n2
        rr = sqrt((xyz1(1,i)-xyz2(1,j))**2+(xyz1(2,i)-xyz2(2,j))**2+(xyz1(3,i)-xyz2(3,j))**2)
        if (rr .lt. rmin) rmin = rr
      end do
    end do
  end subroutine rcma

!========================================================================================!

  subroutine rotmat(rxyz,rot)
    !***********************************************
    !* rotation matrix R = Rx(4)*Ry(5)*Rz(6)
    !***********************************************
    real(wp),intent(in)  :: rxyz(6)
    real(wp),intent(out) :: rot(3,3)
    real(wp) :: r1(3,3),r2(3,3),r3(3,3)
    r1 = 0.0_wp
    r1(1,1) = 1.0_wp
    r1(2,2) = cos(rxyz(4))
    r1(2,3) = -sin(rxyz(4))
    r1(3,2) = sin(rxyz(4))
    r1(3,3) = cos(rxyz(4))
    r2 = 0.0_wp
    r2(2,2) = 1.0_wp
    r2(1,1) = cos(rxyz(5))
    r2(1,3) = sin(rxyz(5))
    r2(3,1) = -sin(rxyz(5))
    r2(3,3) = cos(rxyz(5))
    r3 = 0.0_wp
    r3(3,3) = 1.0_wp
    r3(1,1) = cos(rxyz(6))
    r3(1,2) = -sin(rxyz(6))
    r3(2,1) = sin(rxyz(6))
    r3(2,2) = cos(rxyz(6))
    rot = matmul(r1,matmul(r2,r3))
  end subroutine rotmat

!========================================================================================!

  subroutine move(n2,nl2,xyz,cl,xyznew,clnew,rotm,coord)
    !***********************************************
    !* rotate (about the CMA) and shift molecule B
    !* together with its LMO charge centers
    !***********************************************
    integer,intent(in)   :: n2,nl2
    real(wp),intent(in)  :: xyz(3,n2),cl(4,n2*10)
    real(wp),intent(out) :: xyznew(3,n2),clnew(4,n2*10)
    real(wp),intent(in)  :: coord(6),rotm(3,3)
    integer :: k
    xyznew = xyz
    clnew = cl
    if (sum(abs(coord)) .lt. 1.d-12) return
    xyznew = matmul(rotm,xyznew)
    clnew(1:3,1:nl2) = matmul(rotm,clnew(1:3,1:nl2))
    do k = 1,3
      xyznew(k,1:n2) = xyznew(k,1:n2)+coord(k)
      clnew(k,1:nl2) = clnew(k,1:nl2)+coord(k)
    end do
  end subroutine move

  subroutine move2(n2,xyz,xyznew,coord)
    integer,intent(in)   :: n2
    real(wp),intent(in)  :: xyz(3,n2),coord(6)
    real(wp),intent(out) :: xyznew(3,n2)
    real(wp) :: rotm(3,3)
    integer :: k
    call rotmat(coord,rotm)
    xyznew = matmul(rotm,xyz)
    do k = 1,3
      xyznew(k,1:n2) = xyznew(k,1:n2)+coord(k)
    end do
  end subroutine move2

!========================================================================================!

  subroutine rnorm(i,j,found,r)
    !***********************************************
    !* distance of two genes (translation+angles)
    !***********************************************
    integer,intent(in)   :: i,j
    real(wp),intent(in)  :: found(7,*)
    real(wp),intent(out) :: r
    integer :: k
    r = 0
    do k = 1,3
      r = r+(found(k,i)-found(k,j))**2
    end do
    do k = 4,6
      r = r+sin(0.5*(found(k,i)-found(k,j)))**2
    end do
    r = sqrt(r)
  end subroutine rnorm

!========================================================================================!

  subroutine rand6(mut,r1,r2,d)
    !***********************************************
    !* random mutation of a gene d(1:6); with
    !* probability 1-mut, shifts by up to r1 and
    !* rotates by up to r2 (5% flips by 180 deg)
    !***********************************************
    real(wp),intent(in)    :: mut,r1,r2
    real(wp),intent(inout) :: d(6)
    real(wp) :: x,y,yy,f
    integer :: j,k

    call random_number(yy)
    if (yy .lt. mut) return

    do j = 1,3
      call random_number(y)
      f = 1.0
      if (y .lt. 0.5) f = -1.0
      call random_number(x)
      d(j) = d(j)+f*r1*x
    end do

    if (yy .gt. 0.95) then
      call random_number(y)
      call irand(3,k)
      f = 0.5
      if (y .lt. 0.5) f = -0.5
      d(k+3) = d(k+3)+f*pi2
      return
    end if

    do j = 4,6
      call random_number(y)
      f = 1.0
      if (y .lt. 0.5) f = -1.0
      call random_number(x)
      d(j) = d(j)+f*r2*x
      if (d(j) .gt. pi2) d(j) = d(j)-pi2
      if (d(j) .lt. 0) d(j) = d(j)+pi2
    end do
  end subroutine rand6

  subroutine irand(n,k)
    !> random integer 1 <= k <= n
    integer,intent(in)  :: n
    integer,intent(out) :: k
    real :: x,nx
    call random_number(x)
    nx = n*x+1
    k = int(nx)
    if (k .gt. n) k = n
  end subroutine irand

  subroutine crossover(r,d)
    !> random crossover weights (single precision, as in xtb)
    real(wp),intent(in)    :: r
    real(wp),intent(inout) :: d(6)
    integer :: j
    real :: x
    do j = 1,6
      call random_number(x)
      if (x .gt. r) then
        d(j) = x
      else
        d(j) = 0.0_wp
      end if
    end do
  end subroutine crossover

!========================================================================================!

  subroutine sort6(n,e)
    !> sort genes e(1:7,:) by energy e(7,:)
    integer,intent(in)     :: n
    real(wp),intent(inout) :: e(7,n)
    integer :: ind(n),i
    real(wp) :: tmp(n),ftmp(7,n)
    ftmp = e
    do i = 1,n
      ind(i) = i
      tmp(i) = e(7,i)
    end do
    call qsort(tmp,1,n,ind)
    do i = 1,n
      e(1:7,i) = ftmp(1:7,ind(i))
    end do
  end subroutine sort6

!========================================================================================!

  subroutine mrec(molcount,xyz,nat,at,molvec)
    !***********************************************
    !* split xyz into non-bonded fragments
    !* molcount : number of fragments
    !* molvec   : fragment index per atom
    !***********************************************
    integer,intent(in)   :: nat,at(nat)
    real(wp),intent(in)  :: xyz(3,nat)
    integer,intent(out)  :: molvec(nat),molcount
    real(wp) :: cn(nat),bond(nat,nat)
    integer :: i
    logical :: taken(nat)

    molvec = 0
    molcount = 1
    taken = .false.
    call xcoord(nat,at,xyz,cn,bond)
    do i = 1,nat
      if (.not.taken(i)) then
        molvec(i) = molcount
        taken(i) = .true.
        call neighbours(i,taken,nat,cn,bond,molvec,molcount)
        molcount = molcount+1
      end if
    end do
    molcount = molcount-1
  end subroutine mrec

  recursive subroutine neighbours(i,taken,nat,cn,bond,molvec,molcnt)
    integer,intent(in)     :: i,nat,molcnt
    real(wp),intent(inout) :: cn(nat),bond(nat,nat)
    integer,intent(inout)  :: molvec(nat)
    logical,intent(inout)  :: taken(nat)
    integer :: j,icn,k
    icn = nint(cn(i))
    do k = 1,icn
      j = maxloc(bond(:,i),1)
      bond(j,i) = 0.0_wp
      if (i .eq. j) cycle
      if (.not.taken(j)) then
        molvec(j) = molcnt
        taken(j) = .true.
        call neighbours(j,taken,nat,cn,bond,molvec,molcnt)
      end if
    end do
  end subroutine neighbours

  subroutine xcoord(natoms,iz,xyz,cn,bond)
    !***********************************************
    !* exponential counting CN with scaled D3 radii
    !***********************************************
    integer,intent(in)   :: natoms,iz(natoms)
    real(wp),intent(in)  :: xyz(3,natoms)
    real(wp),intent(out) :: bond(natoms,natoms),cn(natoms)
    integer :: i,iat
    real(wp) :: r,rco,damp
    bond = 0.0_wp
    cn = 0.0_wp
    do i = 1,natoms
      do iat = 1,natoms
        if (iat .eq. i) cycle
        r = sqrt(sum((xyz(:,iat)-xyz(:,i))**2))
        !> the 0.75 scaling keeps the grid-point clusters small
        rco = (rcov(iz(i))+rcov(iz(iat)))*0.75_wp
        damp = 1.0_wp/(1.0_wp+exp(-16.0_wp*(rco/r-1.0_wp)))
        bond(iat,i) = damp
        cn(i) = cn(i)+damp
      end do
    end do
  end subroutine xcoord

!========================================================================================!

  subroutine wrc2(fname,mode,n1,n2,at1,at2,xyz1,xyz2,nf,found)
    !***********************************************
    !* write the nf first genes as XYZ trajectory;
    !* mode /= 0 draws B as boron for visualization
    !***********************************************
    use strucrd,only:i2e
    character(len=*),intent(in) :: fname
    integer,intent(in)  :: nf,mode,n1,n2,at1(n1),at2(n2)
    real(wp),intent(in) :: found(7,*),xyz1(3,n1),xyz2(3,n2)
    real(wp) :: c2(3,n2)
    integer :: i,j,ich
    open (newunit=ich,file=fname)
    do i = 1,nf
      write (ich,*) n1+n2
      write (ich,'('' SCF done '',F16.8)') found(7,i)
      call move2(n2,xyz2,c2,found(1:6,i))
      do j = 1,n1
        write (ich,'(a2,3F24.10)') i2e(at1(j),'nc'),xyz1(1:3,j)*autoang
      end do
      do j = 1,n2
        if (mode .eq. 0) then
          write (ich,'(a2,3F24.10)') i2e(at2(j),'nc'),c2(1:3,j)*autoang
        else
          write (ich,'(a2,3F24.10)') i2e(5,'nc'),c2(1:3,j)*autoang
        end if
      end do
    end do
    close (ich)
  end subroutine wrc2

  subroutine wrc3(fname,n1,n2,at1,at2,xyz1,xyz2,i,found)
    !> best gene as XYZ; i == 0 converts its energy to kcal/mol
    use strucrd,only:i2e
    character(len=*),intent(in) :: fname
    integer,intent(in)  :: i,n1,n2,at1(n1),at2(n2)
    real(wp),intent(in) :: found(7,*),xyz1(3,n1),xyz2(3,n2)
    real(wp) :: c2(3,n2),f
    integer :: j,ich
    f = 1.0_wp
    if (i .eq. 0) f = au
    open (newunit=ich,file=fname)
    write (ich,*) n1+n2
    write (ich,'('' SCF done '',F16.8)') found(7,1)*f
    call move2(n2,xyz2,c2,found(1:6,1))
    do j = 1,n1
      write (ich,'(a2,3F24.10)') i2e(at1(j),'nc'),xyz1(1:3,j)*autoang
    end do
    do j = 1,n2
      write (ich,'(a2,3F24.10)') i2e(at2(j),'nc'),c2(1:3,j)*autoang
    end do
    close (ich)
  end subroutine wrc3

  subroutine wr_grid(fname,n1,n2,at1,element,xyz1,xyz2)
    !> debug output of a probe grid around molecule A
    use strucrd,only:i2e
    character(len=*),intent(in) :: fname
    integer,intent(in)  :: n1,n2,at1(n1),element
    real(wp),intent(in) :: xyz1(3,n1),xyz2(3,n2)
    integer :: j,ich
    open (newunit=ich,file=fname)
    write (ich,*) n1+n2
    write (ich,*)
    do j = 1,n1
      write (ich,'(a2,3F24.10)') i2e(at1(j),'nc'),xyz1(1:3,j)*autoang
    end do
    do j = 1,n2
      write (ich,'(a2,3F24.10)') i2e(element,'nc'),xyz2(1:3,j)*autoang
    end do
    close (ich)
  end subroutine wr_grid

!========================================================================================!

  subroutine rsp_packed(a,n,w,z)
    !***********************************************
    !* eigenvalues w (ascending) and vectors z of a
    !* packed lower-triangular symmetric matrix a
    !***********************************************
    integer,intent(in)   :: n
    real(wp),intent(in)  :: a(n*(n+1)/2)
    real(wp),intent(out) :: w(n),z(n,n)
    real(wp) :: ap(n*(n+1)/2),work(3*n)
    integer :: info
    external :: dspev
    ap = a
    call dspev('V','U',n,ap,w,z,n,work,info)
  end subroutine rsp_packed

!========================================================================================!

  subroutine cavity_egrad(nat,xyz,nwalls,walls,efix)
    !*****************************************************
    !* xtb wall energy with w = R0/radius, R0 = max(radius):
    !*  polynomial (|w*(r-c)|/R0)^alpha
    !*  logfermi   kT*log(1+exp(beta*(|w*(r-c)|-R0)))
    !*****************************************************
    integer,intent(in)  :: nat,nwalls
    real(wp),intent(in) :: xyz(3,nat)
    type(aiss_wall),intent(in) :: walls(:)
    real(wp),intent(inout) :: efix
    integer :: i,j,iat,nlist
    real(wp) :: r(3),w(3),dist,r0

    do i = 1,nwalls
      r0 = maxval(walls(i)%radius)
      w = r0/walls(i)%radius
      if (allocated(walls(i)%list)) then
        nlist = size(walls(i)%list,1)
      else
        nlist = nat
      end if
      do j = 1,nlist
        if (allocated(walls(i)%list)) then
          iat = walls(i)%list(j)
          if (iat > nat) cycle
        else
          iat = j
        end if
        r = w*(xyz(:,iat)-walls(i)%center)
        dist = sqrt(sum(r**2))
        if (walls(i)%logfermi) then
          efix = efix+kB*walls(i)%temp*log(1.0_wp+exp(walls(i)%beta*(dist-r0)))
        else
          efix = efix+(dist/r0)**walls(i)%alpha
        end if
      end do
    end do
  end subroutine cavity_egrad

  subroutine cavitye(n,xyz,rabc,e)
    !*****************************************************
    !* old QCG-mode ellipsoid, ((x/a)^2+(y/b)^2+(z/c)^2)^15
    !*****************************************************
    integer,intent(in)  :: n
    real(wp),intent(in) :: xyz(3,n),rabc(3)
    real(wp),intent(inout) :: e
    integer,parameter :: a = 15
    integer :: i
    real(wp) :: rx,ry,rz
    do i = 1,n
      rx = (xyz(1,i)/rabc(1))**2
      ry = (xyz(2,i)/rabc(2))**2
      rz = (xyz(3,i)/rabc(3))**2
      e = e+(rx+ry+rz)**a
    end do
  end subroutine cavitye

!========================================================================================!
end module aiss_param
