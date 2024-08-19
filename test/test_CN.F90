module test_CN
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_testmol
  use strucrd
  use crest_cn_module
  use miscdata,only:RCOV  !> D3 covalent radii
  implicit none
  private

  public :: collect_CN

  real(wp),parameter :: thr = 1.e-5

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for calculating coordination numbers
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_CN(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
    new_unittest("CN(exp) calculation           ",test_cn_exp), &
    new_unittest("CN(gfn) calculation           ",test_cn_gfn), &
    new_unittest("CN(erf) calculation           ",test_cn_erf), &
!    new_unittest("CN(erf/EN) calculation        ",test_cn_erf_en), &
    new_unittest("CN(erf/D4) calculation        ",test_cn_d4),  &
    new_unittest("CN calculation via coord type ",test_cn_coord)  &
    ]
!&>
  end subroutine collect_CN

!========================================================================================!

  subroutine test_cn_exp(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: cn(:)
    integer :: io
    logical :: wr,pr
    real(wp),allocatable :: bond(:,:)
!&<
    real(wp),parameter :: cn_ref(24) = reshape([&
    & 4.055891127145659_wp, &
    & 3.288889882015843_wp, &
    & 3.513090304082208_wp, &
    & 2.220376901140114_wp, &
    & 3.637790405819880_wp, &
    & 3.482378921657057_wp, &
    & 3.173753862846316_wp, &
    & 1.067653963029709_wp, &
    & 3.162194988565337_wp, &
    & 3.219594664750161_wp, &
    & 1.066642405807029_wp, &
    & 3.149745868662636_wp, &
    & 4.079048157993642_wp, &
    & 4.093458099252320_wp, &
    & 0.998258561306359_wp, &
    & 0.998521112348117_wp, &
    & 0.998518471567465_wp, &
    & 1.000789208143010_wp, &
    & 0.998145405602563_wp, &
    & 0.998002933550897_wp, &
    & 0.998000732079435_wp, &
    & 0.998213943231262_wp, &
    & 0.998095785577331_wp, &
    & 0.998094156726710_wp  ],shape(cn_ref))
!&>

    !> setup
    call get_testmol('caffeine',mol)
!    call mol%write('coord')
    allocate (cn(mol%nat))
    call calculate_CN(mol%nat,mol%at,mol%xyz,cn)
!   write(*,'(F25.15)') cn
!   write(*,'((F20.15,"_wp,")," &")') cn
    if (any(abs(cn-cn_ref) > thr)) then
      call test_failed(error,"Coordination number does not match")
      print'(es21.14)',cn
      print'("---")'
      print'(es21.14)',cn_ref
      print'("---")'
      print'(es21.14)',cn-cn_ref
    end if


    deallocate (cn)
  end subroutine test_cn_exp

!========================================================================================!

  subroutine test_cn_erf(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: cn(:)
    integer :: io
    logical :: wr,pr
    real(wp),allocatable :: bond(:,:)
!&<
    real(wp),parameter :: cn_ref(24) = reshape([&
    & 3.622643648420508_wp, &
    & 3.405718465347867_wp, &
    & 3.507245465631688_wp, &
    & 2.486310140446890_wp, &
    & 4.082546940409358_wp, &
    & 4.081594016649126_wp, &
    & 3.380301062250944_wp, &
    & 1.234253454564075_wp, &
    & 3.343427434675905_wp, &
    & 3.487619595812383_wp, &
    & 1.249505697640035_wp, &
    & 3.316227492282645_wp, &
    & 3.698556251266844_wp, &
    & 3.736560853778984_wp, &
    & 0.835607601692369_wp, &
    & 0.837012229249737_wp, &
    & 0.837010957128861_wp, &
    & 0.850441048113596_wp, &
    & 0.830336741440498_wp, &
    & 0.832966244796499_wp, &
    & 0.832940283784449_wp, &
    & 0.830257575136063_wp, &
    & 0.833483436263467_wp, &
    & 0.833442705428422_wp  &
    &  ],shape(cn_ref))
!&>
    !> setup
    call get_testmol('caffeine',mol)
    allocate (cn(mol%nat))
    call calculate_CN(mol%nat,mol%at,mol%xyz,cn,cntype='erf')
    !write(*,'((F20.15,"_wp,")," &")') cn
    if (any(abs(cn-cn_ref) > thr)) then
      call test_failed(error,"Coordination number does not match")
      print'(es21.14)',cn
      print'("---")'
      print'(es21.14)',cn_ref
      print'("---")'
      print'(es21.14)',cn-cn_ref
    end if


    deallocate (cn)
  end subroutine test_cn_erf

!========================================================================================!

  subroutine test_cn_gfn(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: cn(:),dcn(:,:,:)
    integer :: io
    logical :: wr,pr
    real(wp),allocatable :: bond(:,:)
!&<
    real(wp),parameter :: cn_ref(24) = reshape([&
    & 4.080787717742611_wp, &
    & 3.605947418297714_wp, &
    & 3.707376933085117_wp, &
    & 2.507369342111704_wp, &
    & 4.084010812746403_wp, &
    & 4.014235318357683_wp, &
    & 3.479735006842208_wp, &
    & 1.228166602213276_wp, &
    & 3.513131818377105_wp, &
    & 3.565453679633609_wp, &
    & 1.243466609796757_wp, &
    & 3.486018352311430_wp, &
    & 4.140527219855912_wp, &
    & 4.172156631475451_wp, &
    & 0.997681762757485_wp, &
    & 0.996836574428603_wp, &
    & 0.996829830734909_wp, &
    & 1.009883794482724_wp, &
    & 0.998871095489539_wp, &
    & 0.995131864627251_wp, &
    & 0.995121246092223_wp, &
    & 0.999619334837355_wp, &
    & 0.996233442007880_wp, &
    & 0.996201406389150_wp &
    &  ],shape(cn_ref))
!&>
    !> setup
    call get_testmol('caffeine',mol)
    allocate (cn(mol%nat),dcn(3,mol%nat,mol%nat))
    call calculate_CN(mol%nat,mol%at,mol%xyz,cn,cntype='gfn',dcndr=dcn)
    !write(*,'((F20.15,"_wp,")," &")') cn
    !write(*,'(3(F20.15,"_wp,")," &")') dcn
    if (any(abs(cn-cn_ref) > thr)) then
      call test_failed(error,"Coordination number does not match")
      print'(es21.14)',cn
      print'("---")'
      print'(es21.14)',cn_ref
      print'("---")'
      print'(es21.14)',cn-cn_ref
    end if

    deallocate (dcn,cn)
  end subroutine test_cn_gfn

!========================================================================================!

  subroutine test_cn_erf_en(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: cn(:)
    integer :: io
    logical :: wr,pr
    real(wp),allocatable :: bond(:,:)
!&<
    real(wp),parameter :: cn_ref(24) = reshape([&
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp, &
    & 0.000000000_wp  &
    &  ],shape(cn_ref))
!&>
    !> setup
    call get_testmol('caffeine',mol)
    allocate (cn(mol%nat))
    call calculate_CN(mol%nat,mol%at,mol%xyz,cn,cntype='erf_en')
    write(*,'((F20.15,"_wp,")," &")') cn
!    if (any(abs(cn-cn_ref) > thr)) then
!      call test_failed(error,"Coordination number does not match")
!      print'(es21.14)',cn
!      print'("---")'
!      print'(es21.14)',cn_ref
!      print'("---")'
!      print'(es21.14)',cn-cn_ref
!    end if
!

    deallocate (cn)
  end subroutine test_cn_erf_en

!========================================================================================!

  subroutine test_cn_d4(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: cn(:)
    integer :: io
    logical :: wr,pr
    real(wp),allocatable :: bond(:,:)
!&<
    real(wp),parameter :: cn_ref(24) = reshape([&
    & 3.687259184149240_wp, &
    & 2.867262205332734_wp, &
    & 3.150724166180586_wp, &
    & 1.893569777782766_wp, &
    & 3.203211463839986_wp, &
    & 3.080099878320893_wp, &
    & 2.766851050499062_wp, &
    & 0.858277785022996_wp, &
    & 2.742941624230694_wp, &
    & 2.714022969151142_wp, &
    & 0.857488513179795_wp, &
    & 2.737922550140081_wp, &
    & 3.693785224461224_wp, &
    & 3.700081658703103_wp, &
    & 0.924156979713941_wp, &
    & 0.924169842395564_wp, &
    & 0.924168470344702_wp, &
    & 0.925321663943406_wp, &
    & 0.924123074464423_wp, &
    & 0.923951547270359_wp, &
    & 0.923952329637323_wp, &
    & 0.924191441370459_wp, &
    & 0.923896277881526_wp, &
    & 0.923905844557051_wp &
    &  ],shape(cn_ref))
!&>
    !> setup
    call get_testmol('caffeine',mol)
    allocate (cn(mol%nat))
    call calculate_CN(mol%nat,mol%at,mol%xyz,cn,cntype='d4',cnthr=sqrt(1600.0_wp))
    !write(*,'((F20.15,"_wp,")," &")') cn
    if (any(abs(cn-cn_ref) > thr)) then
      call test_failed(error,"Coordination number does not match")
      print'(es21.14)',cn
      print'("---")'
      print'(es21.14)',cn_ref
      print'("---")'
      print'(es21.14)',cn-cn_ref
    end if

    deallocate (cn)
  end subroutine test_cn_d4

!========================================================================================!

  subroutine test_cn_coord(error)
    type(error_type),allocatable,intent(out) :: error
    type(coord) :: mol
    real(wp),allocatable :: cn(:)
    integer :: io
    logical :: wr,pr
    real(wp),allocatable :: dcndr(:,:,:)
!&<
    real(wp),parameter :: cn_ref(24) = reshape([&
    & 3.687259184149240_wp, &
    & 2.867262205332734_wp, &
    & 3.150724166180586_wp, &
    & 1.893569777782766_wp, &
    & 3.203211463839986_wp, &
    & 3.080099878320893_wp, &
    & 2.766851050499062_wp, &
    & 0.858277785022996_wp, &
    & 2.742941624230694_wp, &
    & 2.714022969151142_wp, &
    & 0.857488513179795_wp, &
    & 2.737922550140081_wp, &
    & 3.693785224461224_wp, &
    & 3.700081658703103_wp, &
    & 0.924156979713941_wp, &
    & 0.924169842395564_wp, &
    & 0.924168470344702_wp, &
    & 0.925321663943406_wp, &
    & 0.924123074464423_wp, &
    & 0.923951547270359_wp, &
    & 0.923952329637323_wp, &
    & 0.924191441370459_wp, &
    & 0.923896277881526_wp, &
    & 0.923905844557051_wp &
    &  ],shape(cn_ref))
!&>
    !> setup
    call get_testmol('caffeine',mol)
    call mol%get_cn(cn,'cov')
    !write(*,'((F20.15,"_wp,")," &")') cn
    if (any(abs(cn-cn_ref) > thr)) then
      call test_failed(error,"Coordination number does not match")
      print'(es21.14)',cn
      print'("---")'
      print'(es21.14)',cn_ref
      print'("---")'
      print'(es21.14)',cn-cn_ref
    end if

    deallocate (cn)
  end subroutine test_cn_coord

!========================================================================================!
!========================================================================================!
end module test_CN
