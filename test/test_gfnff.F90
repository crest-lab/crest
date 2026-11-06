module test_gfnff
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use crest_testmol
  implicit none
  private

  public :: collect_gfnff

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
  real(wp),parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for using gfnff in crest
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_gfnff(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
#ifdef WITH_GFNFF
    new_unittest("Compiled gfnff subproject     ",test_compiled_gfnff), &
    new_unittest("GFN-FF singlepoint            ",test_gfnff_sp), &
    new_unittest("GFN-FF singlepoint (cation)   ",test_gfnff_sp_cation), &
    new_unittest("GFN-FF singlepoint (anion)    ",test_gfnff_sp_anion), &
    new_unittest("GFN-FF singlepoint (ALPB)     ",test_gfnff_sp_alpb) &
#else
    new_unittest("Compiled gfnff subproject",test_compiled_gfnff,should_fail=.true.) &
#endif
    ]
!&>
  end subroutine collect_gfnff

  subroutine test_compiled_gfnff(error)
    type(error_type),allocatable,intent(out) :: error
#ifndef WITH_GFNFF
    write(*,'("       ...")') 'gfnff not compiled, expecting fail.'
    allocate (error)
#endif
  end subroutine test_compiled_gfnff

  subroutine test_gfnff_sp(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -4.672792615407980_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &   0.005301577528380_wp,    0.000273969963676_wp,    0.000002235966901_wp, & 
    &   0.008166049302148_wp,   -0.008220839652479_wp,   -0.000025577423357_wp, & 
    &  -0.003078315962938_wp,   -0.009433015660225_wp,    0.000033248957362_wp, & 
    &   0.009919813036231_wp,    0.008086621308606_wp,   -0.000022035638450_wp, & 
    &  -0.015632607591379_wp,   -0.026672383478048_wp,    0.000004837609574_wp, & 
    &   0.014525647437915_wp,   -0.001976836648137_wp,    0.000067168707615_wp, & 
    &   0.006146656755674_wp,    0.009561079703764_wp,   -0.000017505319062_wp, & 
    &  -0.008820848527167_wp,   -0.001068632036479_wp,   -0.000078000869434_wp, & 
    &  -0.000983352637051_wp,    0.014873587850193_wp,    0.000032976002056_wp, & 
    &  -0.006683050923820_wp,    0.007422835964584_wp,   -0.000019221292477_wp, & 
    &   0.012839290444322_wp,   -0.012743002931932_wp,   -0.000039643202321_wp, & 
    &  -0.023422684595081_wp,    0.021005864224867_wp,   -0.000002459565985_wp, & 
    &  -0.001884047168835_wp,   -0.003906629596872_wp,   -0.000013746283758_wp, & 
    &  -0.003754971577627_wp,    0.003730231633401_wp,   -0.000073759305057_wp, & 
    &   0.000742834486661_wp,    0.003621860437988_wp,    0.000003807409655_wp, & 
    &   0.001069305640750_wp,   -0.000350573335716_wp,    0.003705264819091_wp, & 
    &   0.001070928608593_wp,   -0.000349783768115_wp,   -0.003711454166132_wp, & 
    &  -0.002984448998131_wp,    0.000421235644680_wp,    0.000013800906585_wp, & 
    &   0.004499275275005_wp,    0.000660471751639_wp,    0.000002343070942_wp, & 
    &   0.000371387054650_wp,    0.001498977490928_wp,    0.003776574488090_wp, & 
    &   0.000381320992120_wp,    0.001507606004227_wp,   -0.003766956258310_wp, & 
    &  -0.001010470679296_wp,   -0.004606701841933_wp,    0.000057155046565_wp, & 
    &   0.001617337942360_wp,   -0.001636909416951_wp,    0.003219098019704_wp, & 
    &   0.001603374156518_wp,   -0.001699033611664_wp,   -0.003148151679798_wp & 
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfnff')
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-7_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient of energy does not match")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_gfnff_sp

!=======================================================================================!

  subroutine test_gfnff_sp_cation(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -3.329316204948430_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &    0.006025113616110_wp,  -0.000476041539293_wp,   0.000002594729221_wp, &
    &   0.002976419911185_wp,   0.000679240987425_wp,  -0.000027434135534_wp, &
    &  -0.007741249228758_wp,  -0.010698615223025_wp,   0.000033493626850_wp, &
    &   0.015675798394840_wp,   0.001864038271447_wp,  -0.000014511182788_wp, &
    &  -0.034151774354292_wp,  -0.007653966260220_wp,   0.000003680323479_wp, &
    &   0.022999436444320_wp,  -0.004730111785708_wp,   0.000060519347357_wp, &
    &   0.006126360638611_wp,   0.009306188454984_wp,  -0.000009596058284_wp, &
    &  -0.013550831209723_wp,  -0.002922580821436_wp,  -0.000095354457018_wp, &
    &  -0.002031960041829_wp,   0.013328804263952_wp,   0.000044228493753_wp, &
    &  -0.007664605104266_wp,   0.014627940401588_wp,  -0.000016489263802_wp, &
    &   0.018754969893106_wp,  -0.017773375742745_wp,  -0.000041517059516_wp, &
    &  -0.007820308319431_wp,   0.003713980391848_wp,  -0.000001085232997_wp, &
    &  -0.003458102422879_wp,  -0.003295384699656_wp,  -0.000018978176550_wp, &
    &  -0.003546064904121_wp,   0.005171739411741_wp,  -0.000084843666171_wp, &
    &   0.001188602766869_wp,   0.002934130726579_wp,   0.000003784320013_wp, &
    &   0.001506247126530_wp,   0.000290262093640_wp,   0.002739113228752_wp, &
    &   0.001506604239896_wp,   0.000291343554975_wp,  -0.002741664315015_wp, &
    &  -0.003214661896317_wp,  -0.001822644901899_wp,   0.000012841210822_wp, &
    &   0.003609901458510_wp,   0.002382972247914_wp,   0.000004307221005_wp, &
    &   0.000282820078859_wp,   0.001032394332879_wp,   0.002665731125889_wp, &
    &   0.000293959715940_wp,   0.001040386594915_wp,  -0.002655070180927_wp, &
    &   0.000911301550457_wp,  -0.004660465773164_wp,   0.000059733158755_wp, &
    &   0.000667168759801_wp,  -0.001290577497866_wp,   0.001915244081988_wp, &
    &   0.000654852886582_wp,  -0.001339657488875_wp,  -0.001838727139279_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfnff')
    sett%chrg = 1
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-7_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient of energy does not match")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_gfnff_sp_cation

!========================================================================================!

  subroutine test_gfnff_sp_anion(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    !> Reference values updated for gfnff v0.1.1: the neighbour-list and
    !> Wigner-Seitz cell (WSC) rework changed energetics for charged systems.
    !> Both energy and gradient references were regenerated from the v0.1.1 run.
    real(wp),parameter :: e_ref = -5.828532820056295_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &  3.94230543214767e-03_wp,  3.52030063768979e-04_wp,  1.81683057095039e-06_wp, &
    &  1.51705975544829e-02_wp, -8.59644901931333e-03_wp, -1.30289855263585e-05_wp, &
    & -3.32620665073756e-03_wp, -1.27838095920201e-02_wp,  3.07797052876764e-05_wp, &
    &  4.93740063275844e-03_wp,  4.06512386391377e-03_wp, -2.14659925847723e-05_wp, &
    & -2.31495929321227e-02_wp, -2.87916805922025e-02_wp, -3.22455640084147e-06_wp, &
    &  2.05728120146893e-02_wp,  1.07842246816042e-02_wp,  5.58592835861545e-05_wp, &
    & -3.57568867600167e-03_wp, -1.87642796677282e-03_wp, -2.53540979489806e-05_wp, &
    &  5.05671257009310e-03_wp,  4.08216306711967e-03_wp, -5.09893176526027e-05_wp, &
    & -2.48954902302761e-03_wp,  1.56793979443191e-02_wp,  3.18591292280463e-05_wp, &
    & -4.61712459034876e-03_wp,  7.38892241260220e-03_wp, -2.69984802863847e-05_wp, &
    &  6.67104429871043e-03_wp, -7.53862115125029e-03_wp, -3.64266870952139e-05_wp, &
    & -2.31759261766794e-02_wp,  1.91389704080971e-02_wp,  6.34113639207092e-06_wp, &
    &  6.05908705535911e-06_wp, -3.03220605329460e-03_wp, -1.40325040899702e-05_wp, &
    & -3.98348638058101e-03_wp,  1.54809264373304e-03_wp, -5.64650422844910e-05_wp, &
    &  1.71517910778865e-03_wp,  2.42726129127660e-03_wp,  2.86445331615998e-06_wp, &
    &  1.89583331406671e-03_wp, -5.18309297446946e-04_wp,  3.05090837898127e-03_wp, &
    &  1.89876480153798e-03_wp, -5.18627636530015e-04_wp, -3.05758305240559e-03_wp, &
    & -2.71548970167826e-03_wp,  1.07706281565603e-03_wp,  1.41909698337930e-05_wp, &
    &  3.43414680170101e-03_wp, -1.61582876165697e-04_wp,  3.10793153391831e-06_wp, &
    & -4.93795501029164e-04_wp,  9.08335390055891e-04_wp,  3.28781757277625e-03_wp, &
    & -4.84391692211688e-04_wp,  9.16713856580635e-04_wp, -3.27867969575159e-03_wp, &
    & -1.25224049764651e-03_wp, -3.28460279557860e-03_wp,  4.53019960043941e-05_wp, &
    &  1.98710257887371e-03_wp, -6.04718036417461e-04_wp,  2.85313550792252e-03_wp, &
    &  1.97553362815908e-03_wp, -6.61263421734836e-04_wp, -2.79973448340641e-03_wp  &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfnff')
    sett%chrg = -1
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-7_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient of energy does not match")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_gfnff_sp_anion

!========================================================================================!

  subroutine test_gfnff_sp_alpb(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -4.689906438468960_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &   0.005946318956743_wp,   0.000224845414918_wp,   0.000002359259418_wp, &
    &   0.008276028708202_wp,  -0.008052841871053_wp,  -0.000024497556066_wp, &
    &  -0.002933422656266_wp,  -0.008965806356044_wp,   0.000032565115728_wp, &
    &   0.009327430055380_wp,   0.007163077799008_wp,  -0.000020616233726_wp, &
    &  -0.015612398376933_wp,  -0.026601046126910_wp,   0.000003393882782_wp, &
    &   0.014224603776207_wp,  -0.002328340250813_wp,   0.000066592790248_wp, &
    &   0.006359125001507_wp,   0.009728710502539_wp,  -0.000012692526671_wp, &
    &  -0.008060303189818_wp,  -0.001017324006654_wp,  -0.000067764705906_wp, &
    &  -0.000928875315221_wp,   0.014721274585459_wp,   0.000064282322941_wp, &
    &  -0.007032637800329_wp,   0.007686466080511_wp,  -0.000018779454375_wp, &
    &   0.012172269349249_wp,  -0.012198147638366_wp,  -0.000031535160939_wp, &
    &  -0.023075283175922_wp,   0.020590486278622_wp,   0.000000950293160_wp, &
    &  -0.002527752580930_wp,  -0.004378687256113_wp,  -0.000014408178734_wp, &
    &  -0.003644475629980_wp,   0.004533754258488_wp,  -0.000098967058627_wp, &
    &   0.000763589312794_wp,   0.003493537197421_wp,   0.000003659958242_wp, &
    &   0.001177972834069_wp,  -0.000489791575692_wp,   0.003518465980525_wp, &
    &   0.001179451150639_wp,  -0.000489262746804_wp,  -0.003525390673803_wp, &
    &  -0.002858204604008_wp,  -0.000053706012265_wp,   0.000013389518933_wp, &
    &   0.004179436123396_wp,   0.000474457663764_wp,   0.000002349459544_wp, &
    &   0.000262934613792_wp,   0.001522201563008_wp,   0.003711712203599_wp, &
    &   0.000272333461465_wp,   0.001531295083315_wp,  -0.003702175465263_wp, &
    &  -0.001005498943326_wp,  -0.004218548273957_wp,   0.000052754696087_wp, &
    &   0.001779350047276_wp,  -0.001420257920538_wp,   0.003106515658900_wp, &
    &   0.001758008882015_wp,  -0.001456346391844_wp,  -0.003062164125997_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfnff')
    sett%solvmodel = 'alpb'
    sett%solvent = 'water'
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3(F20.15,"_wp,")," &")') grad
    call check(error,io,0)
    if (allocated(error)) return

    call check(error,energy,e_ref,thr=1e-7_wp)
    if (allocated(error)) return

    if (any(abs(grad-g_ref) > thr2)) then
      call test_failed(error,"Gradient of energy does not match")
      print'(3es21.14)',grad
      print'("---")'
      print'(3es21.14)',g_ref
      print'("---")'
      print'(3es21.14)',grad-g_ref
    end if

    deallocate (grad)
  end subroutine test_gfnff_sp_alpb

!========================================================================================!
!========================================================================================!
end module test_gfnff
