module test_gfnff
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use testmol
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
    real(wp),parameter :: e_ref = -5.722199796041350_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &   0.004468921480142_wp,  -0.000083123134016_wp,   0.000001151477114_wp, &
    &   0.026784240021401_wp,  -0.003317301163864_wp,   0.000007330730819_wp, &
    &  -0.004004820272010_wp,  -0.030509909771367_wp,   0.000021392716860_wp, &
    &  -0.000751769525558_wp,   0.005889324801531_wp,  -0.000022156837684_wp, &
    &  -0.027788611671027_wp,  -0.033864525810969_wp,  -0.000006636279691_wp, &
    &   0.026412943164910_wp,   0.030830470039288_wp,   0.000033495772036_wp, &
    &  -0.017689873797001_wp,  -0.023514680185829_wp,  -0.000029972119517_wp, &
    &   0.025299322245125_wp,   0.011384898053826_wp,   0.000001507574957_wp, &
    &  -0.009184660893020_wp,   0.023372121393537_wp,  -0.000006131704134_wp, &
    &  -0.001241048489363_wp,   0.005561107788044_wp,  -0.000016875281359_wp, &
    &   0.002544027763531_wp,  -0.004003045629465_wp,  -0.000031756330103_wp, &
    &  -0.029084734620249_wp,   0.020333106964323_wp,   0.000017800756382_wp, &
    &  -0.000069352904917_wp,  -0.003305508039344_wp,  -0.000014616855913_wp, &
    &  -0.003685271682472_wp,   0.001732297413329_wp,  -0.000073775627784_wp, &
    &   0.001715180720526_wp,   0.002427255238024_wp,   0.000002864449221_wp, &
    &   0.001895835509322_wp,  -0.000518306207627_wp,   0.003050902872775_wp, &
    &   0.001898766998976_wp,  -0.000518624553391_wp,  -0.003057577541061_wp, &
    &  -0.002685446613991_wp,   0.000991558499999_wp,   0.000012103424716_wp, &
    &   0.003434143710645_wp,  -0.000161577974365_wp,   0.000003107935122_wp, &
    &  -0.000493795408284_wp,   0.000908332349694_wp,   0.003287812553032_wp, &
    &  -0.000484391594052_wp,   0.000916710811928_wp,  -0.003278674680856_wp, &
    &  -0.001252234734786_wp,  -0.003284602467482_wp,   0.000045301985594_wp, &
    &   0.001987099764679_wp,  -0.000604716538950_wp,   0.002853130443235_wp, &
    &   0.001975530827471_wp,  -0.000661261876856_wp,  -0.002799729433761_wp &
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
