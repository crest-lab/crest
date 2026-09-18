module test_irmsd
  use testdrive,only:new_unittest,unittest_type,error_type,test_failed
  use crest_parameters,only:wp,aatoau
  use crest_testmol,only:get_testmol
  use strucrd,only:coord
  use irmsd_module,only:rmsd,irmsd
  implicit none
  private

  public :: collect_irmsd

  real(wp),parameter :: thr       = 1.0e-10_wp
  real(wp),parameter :: thr_loose = 1.0e-3_wp

!========================================================================================!
!> H2 diatomic geometries in Bohr (centered at origin)
!========================================================================================!

  integer,parameter :: nat_h2 = 2
  integer,parameter :: at_h2(nat_h2) = [1,1]
  !> bond = 1.0 Bohr
  real(wp),parameter :: xyz_h2_ref(3,nat_h2) = reshape([ &
    & 0.5_wp, 0.0_wp, 0.0_wp,  &
    &-0.5_wp, 0.0_wp, 0.0_wp], &
    & [3,nat_h2])
  !> bond = 2.0 Bohr → analytical RMSD vs ref = 0.5 Bohr
  real(wp),parameter :: xyz_h2_mol(3,nat_h2) = reshape([ &
    & 1.0_wp, 0.0_wp, 0.0_wp,  &
    &-1.0_wp, 0.0_wp, 0.0_wp], &
    & [3,nat_h2])

!========================================================================================!
!> Fluoxetine (40 atoms, Angstrom) — struc1 and struc2 from the iRMSD paper example.
!> Both represent the same conformer; struc2 has completely scrambled atom order
!> and a random rotation applied. iRMSD must return ≈ 0.
!========================================================================================!

  integer,parameter :: nat_fluo = 40

  !> struc1 atom types (Z)
  integer,parameter :: at_fluo1(nat_fluo) = [ &
    & 6,9,1,1,6,1,6,9,9,7, &
    & 6,1,8,6,6,6,6,6,1,6, &
    & 1,6,1,1,1,1,1,1,1,6, &
    & 1,1,1,1,6,6,1,6,6,6]

  !> struc1 coordinates, Angstrom
  real(wp),parameter :: xyz_fluo1_ang(3,nat_fluo) = reshape([ &
    &-0.0198_wp, 0.2158_wp, 0.5308_wp, &
    &-5.0246_wp,-0.2464_wp, 0.5593_wp, &
    &-2.0207_wp, 0.0761_wp, 3.2888_wp, &
    & 2.0221_wp, 5.0303_wp, 0.4892_wp, &
    &-0.4032_wp, 0.2003_wp, 1.8748_wp, &
    & 1.4542_wp, 5.1890_wp,-1.1937_wp, &
    & 2.1773_wp, 2.6686_wp,-0.7626_wp, &
    &-4.6775_wp, 0.8294_wp, 2.4172_wp, &
    &-4.3667_wp,-1.3218_wp, 2.3289_wp, &
    & 3.0555_wp, 3.8290_wp,-0.9233_wp, &
    & 2.3371_wp, 5.0467_wp,-0.5601_wp, &
    & 3.7051_wp, 1.1947_wp,-0.3434_wp, &
    & 1.3305_wp, 0.2076_wp, 0.3364_wp, &
    &-2.3576_wp, 0.1386_wp,-0.0950_wp, &
    & 3.3367_wp,-3.7475_wp,-2.0944_wp, &
    &-2.7426_wp, 0.0481_wp, 1.2484_wp, &
    &-4.1855_wp,-0.1659_wp, 1.6292_wp, &
    & 2.4651_wp,-1.1750_wp,-1.3319_wp, &
    & 1.7906_wp, 2.6512_wp, 0.2660_wp, &
    & 2.8963_wp, 1.3473_wp,-1.0679_wp, &
    & 1.3203_wp, 2.7752_wp,-1.4405_wp, &
    & 1.8806_wp, 0.1908_wp,-0.9905_wp, &
    &-0.7865_wp, 0.2662_wp,-1.5101_wp, &
    & 3.3527_wp, 1.3925_wp,-2.0651_wp, &
    & 1.1282_wp, 0.3962_wp,-1.7541_wp, &
    & 2.8529_wp,-0.6554_wp,-3.3999_wp, &
    & 2.9547_wp,-4.3108_wp,-0.0530_wp, &
    & 3.8626_wp, 3.7250_wp,-0.3037_wp, &
    & 0.3581_wp, 0.2294_wp, 2.6522_wp, &
    &-1.7511_wp, 0.1165_wp, 2.2345_wp, &
    & 2.9886_wp, 5.9146_wp,-0.7002_wp, &
    & 2.1980_wp,-2.0723_wp, 0.6278_wp, &
    & 3.6193_wp,-2.9009_wp,-4.0530_wp, &
    &-3.1084_wp, 0.1174_wp,-0.8845_wp, &
    & 2.5162_wp,-2.2260_wp,-0.4015_wp, &
    &-1.0075_wp, 0.2216_wp,-0.4508_wp, &
    & 3.6502_wp,-4.7430_wp,-2.3936_wp, &
    & 2.9433_wp,-3.5025_wp,-0.7811_wp, &
    & 2.8862_wp,-1.4373_wp,-2.6475_wp, &
    & 3.3158_wp,-2.7117_wp,-3.0256_wp], &
    & [3,nat_fluo])

  !> struc2 atom types (Z) — scrambled order
  integer,parameter :: at_fluo2(nat_fluo) = [ &
    & 6,1,6,6,6,6,1,1,1,1, &
    & 1,1,6,1,6,6,6,1,9,6, &
    & 1,6,6,9,1,1,1,1,1,1, &
    & 6,9,1,1,6,6,7,8,6,6]

  !> struc2 coordinates, Angstrom — scrambled + rotated
  real(wp),parameter :: xyz_fluo2_ang(3,nat_fluo) = reshape([ &
    &-0.2470_wp,-1.0143_wp,-0.4213_wp, &
    & 0.0828_wp, 4.2564_wp,-4.4517_wp, &
    &-2.1193_wp, 2.6178_wp, 0.0144_wp, &
    & 1.2089_wp,-1.8204_wp, 1.8160_wp, &
    & 0.0619_wp,-0.0804_wp, 0.5643_wp, &
    & 0.9632_wp,-2.7436_wp, 0.7925_wp, &
    &-5.9431_wp, 1.8117_wp, 2.4085_wp, &
    &-3.2933_wp, 0.8099_wp, 0.1366_wp, &
    & 2.4555_wp, 4.5828_wp,-3.8102_wp, &
    &-0.0328_wp,-3.0384_wp,-1.1072_wp, &
    &-3.9856_wp, 3.0777_wp, 1.8789_wp, &
    &-2.4640_wp, 1.2332_wp, 1.6473_wp, &
    &-0.4255_wp, 3.2037_wp,-2.6383_wp, &
    & 1.7581_wp, 2.5331_wp,-0.1112_wp, &
    & 1.3640_wp, 2.9121_wp,-1.0521_wp, &
    & 1.5502_wp,-4.1306_wp, 0.8575_wp, &
    &-5.0208_wp, 1.3357_wp, 2.0622_wp, &
    &-4.4695_wp, 0.9920_wp, 2.9447_wp, &
    & 1.1253_wp,-4.8340_wp, 1.9430_wp, &
    & 2.2338_wp, 3.5914_wp,-1.9111_wp, &
    &-1.7744_wp, 3.4140_wp, 0.6847_wp, &
    & 0.4427_wp, 3.8850_wp,-3.4947_wp, &
    & 0.7599_wp,-0.5019_wp, 1.6995_wp, &
    & 1.2525_wp,-4.8990_wp,-0.2270_wp, &
    &-1.3360_wp, 1.0931_wp,-1.2358_wp, &
    & 1.7763_wp,-2.1151_wp, 2.6974_wp, &
    & 0.9785_wp, 0.2106_wp, 2.4926_wp, &
    &-5.3029_wp, 0.4654_wp, 1.4587_wp, &
    &-0.8129_wp,-0.7679_wp,-1.3112_wp, &
    &-2.6932_wp, 3.0927_wp,-0.7916_wp, &
    &-0.9181_wp, 1.8428_wp,-0.5615_wp, &
    & 2.9112_wp,-4.1100_wp, 0.9336_wp, &
    &-1.4519_wp, 3.0470_wp,-2.9556_wp, &
    & 3.2752_wp, 3.7325_wp,-1.6299_wp, &
    & 1.7755_wp, 4.0721_wp,-3.1352_wp, &
    & 0.0230_wp, 2.6924_wp,-1.4078_wp, &
    &-4.2396_wp, 2.2887_wp, 1.2797_wp, &
    &-0.2303_wp, 1.2522_wp, 0.5528_wp, &
    &-3.0176_wp, 1.6462_wp, 0.7923_wp, &
    & 0.1996_wp,-2.3349_wp,-0.3081_wp], &
    & [3,nat_fluo])

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for rmsd and irmsd
!========================================================================================!
!========================================================================================!

  subroutine collect_irmsd(testsuite)
    !***********************************
    !* Register all irmsd test cases.  *
    !***********************************
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
!&<
    testsuite = [ &
      new_unittest("RMSD self comparison            ",test_rmsd_self),           &
      new_unittest("RMSD H2 bond-stretch (known val)",test_rmsd_h2_bondstretch), &
      new_unittest("iRMSD self comparison           ",test_irmsd_self),          &
      new_unittest("iRMSD scrambled atom order      ",test_irmsd_scrambled)      &
    ]
!&>
  end subroutine collect_irmsd

!========================================================================================!

  subroutine test_rmsd_self(error)
    !*****************************************************
    !* RMSD of a structure compared to itself must be 0. *
    !*****************************************************
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp) :: rmsdval
    call get_testmol('caffeine',mol)
    rmsdval = rmsd(mol,mol)
    if (abs(rmsdval) > 1.0e-6_wp) &
      call test_failed(error,'RMSD(mol,mol) should be 0, got: '//to_str(rmsdval))
  end subroutine test_rmsd_self

!========================================================================================!

  subroutine test_rmsd_h2_bondstretch(error)
    !*************************************************************
    !* Two H2 molecules with different bond lengths (1 and 2 Bohr).
    !* Both centered at origin; no rotation needed for alignment.
    !* Analytical RMSD = sqrt((0.5^2 + 0.5^2)/2) = 0.5 Bohr.
    !*************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: ref,mol
    real(wp) :: rmsdval
    real(wp),parameter :: expected = 0.5_wp

    ref%nat = nat_h2
    allocate(ref%at(nat_h2),ref%xyz(3,nat_h2))
    ref%at  = at_h2
    ref%xyz = xyz_h2_ref

    mol%nat = nat_h2
    allocate(mol%at(nat_h2),mol%xyz(3,nat_h2))
    mol%at  = at_h2
    mol%xyz = xyz_h2_mol

    rmsdval = rmsd(ref,mol)
    if (abs(rmsdval - expected) > thr) &
      call test_failed(error,'Expected RMSD = 0.5 Bohr, got: '//to_str(rmsdval))
  end subroutine test_rmsd_h2_bondstretch

!========================================================================================!

  subroutine test_irmsd_self(error)
    !*******************************************************
    !* iRMSD of a structure compared to itself must be 0.  *
    !*******************************************************
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp) :: rmsdval
    call get_testmol('caffeine',mol)
    rmsdval = irmsd(mol,mol,topocheck=.false.)
    if (rmsdval > thr_loose) &
      call test_failed(error,'iRMSD(mol,mol) should be ~0, got: '//to_str(rmsdval))
  end subroutine test_irmsd_self

!========================================================================================!

  subroutine test_irmsd_scrambled(error)
    !***********************************************************************
    !* Two copies of the same fluoxetine conformer (40 atoms).             *
    !* struc2 has fully randomised atom order and a random rotation.       *
    !* iRMSD must find the correct permutation and return ≈ 0 Å.          *
    !***********************************************************************
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol1,mol2
    real(wp) :: rmsdval

    mol1%nat = nat_fluo
    allocate(mol1%at(nat_fluo),mol1%xyz(3,nat_fluo))
    mol1%at  = at_fluo1
    mol1%xyz = xyz_fluo1_ang * aatoau

    mol2%nat = nat_fluo
    allocate(mol2%at(nat_fluo),mol2%xyz(3,nat_fluo))
    mol2%at  = at_fluo2
    mol2%xyz = xyz_fluo2_ang * aatoau

    rmsdval = irmsd(mol1,mol2,topocheck=.false.)
    if (rmsdval > thr_loose) &
      call test_failed(error,'iRMSD of identical fluoxetine (scrambled) should be ~0, got: '//to_str(rmsdval))
  end subroutine test_irmsd_scrambled

!========================================================================================!

  pure function to_str(x) result(s)
    !> Minimal real→string helper for error messages.
    real(wp),intent(in) :: x
    character(len=32) :: s
    write(s,'(es16.8)') x
    s = adjustl(s)
  end function to_str

end module test_irmsd
