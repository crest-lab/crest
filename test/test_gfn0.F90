module test_gfn0
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use testmol
  implicit none
  private

  public :: collect_gfn0

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
  real(wp),parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for using gfn0 in crest
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_gfn0(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
#ifdef WITH_GFN0
    new_unittest("Compiled gfn0 subproject      ",test_compiled_gfn0), &
    new_unittest("GFN0-xTB singlepoint          ",test_gfn0_sp), &
    new_unittest("GFN0-xTB singlepoint (cation) ",test_gfn0_sp_cation), &
    new_unittest("GFN0-xTB singlepoint (anion)  ",test_gfn0_sp_anion), &
    new_unittest("GFN0-xTB singlepoint (S1)     ",test_gfn0_sp_uhf), &
    new_unittest("GFN0-xTB singlepoint (ALPB)   ",test_gfn0_sp_alpb) &
#else
    new_unittest("Compiled gfn0 subproject",test_compiled_gfn0,should_fail=.true.) &
#endif
    ]
!&>
  end subroutine collect_gfn0

  subroutine test_compiled_gfn0(error)
    type(error_type),allocatable,intent(out) :: error
#ifndef WITH_GFN0
    write (*,'("       ...")') 'gfn0 not compiled, expecting fail.'
    allocate (error)
#endif
  end subroutine test_compiled_gfn0

  subroutine test_gfn0_sp(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -40.908850360158375_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &  -0.008276000880612_wp,   0.000674981651598_wp,   0.000004847527110_wp, &
    &   0.000668264513185_wp,  -0.042047718526109_wp,  -0.000053634043553_wp, &
    &   0.029870640452823_wp,   0.018041726002428_wp,   0.000047108740385_wp, &
    &  -0.014397691809728_wp,  -0.000783036074432_wp,  -0.000025116661451_wp, &
    &  -0.025652094129182_wp,   0.001170992902634_wp,  -0.000039879097860_wp, &
    &  -0.005414581736854_wp,   0.019478157220852_wp,   0.000025849983768_wp, &
    &   0.002171994638264_wp,   0.002763826717346_wp,  -0.000027761994727_wp, &
    &  -0.006339371025545_wp,   0.009614108548993_wp,  -0.000095144233133_wp, &
    &  -0.007195087486222_wp,  -0.023395259838885_wp,   0.000142400425061_wp, &
    &  -0.003510117314291_wp,   0.001640619882830_wp,  -0.000054834076982_wp, &
    &   0.022508811615001_wp,  -0.018116741176823_wp,  -0.000038916654210_wp, &
    &   0.010845831057935_wp,   0.023727337503019_wp,  -0.000033624275410_wp, &
    &   0.009140006411910_wp,   0.001070150449043_wp,  -0.000005140577297_wp, &
    &  -0.008520164815828_wp,  -0.003504121840219_wp,   0.000001694971373_wp, &
    &  -0.001910960997009_wp,   0.002489353607793_wp,   0.000004185369847_wp, &
    &   0.001243151419495_wp,   0.000934470115371_wp,   0.001848577600329_wp, &
    &   0.001245637087981_wp,   0.000931685537315_wp,  -0.001851704582091_wp, &
    &  -0.002704195204404_wp,   0.009583395574317_wp,   0.000027163737872_wp, &
    &   0.005501417797794_wp,   0.002334158463011_wp,   0.000002931622293_wp, &
    &   0.000100709327854_wp,  -0.000135677302777_wp,   0.001510921694029_wp, &
    &   0.000119765751454_wp,  -0.000116613513992_wp,  -0.001513953698281_wp, &
    &  -0.000586838167009_wp,  -0.006799931661327_wp,   0.000077003603818_wp, &
    &   0.000543218118930_wp,   0.000226162233238_wp,   0.000696572111601_wp, &
    &   0.000547655374056_wp,   0.000217973524777_wp,  -0.000649547492491_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn0')
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
!    write (*,'(F25.15)') energy
!    write (*,'(3(F20.15,"_wp,")," &")') grad
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
  end subroutine test_gfn0_sp

!=======================================================================================!

  subroutine test_gfn0_sp_cation(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -39.184258118641779_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &  -0.007842230338035_wp,  -0.000716107768510_wp,   0.000006027708001_wp, &
    &  -0.002549502219587_wp,  -0.006712715165025_wp,  -0.000063584928083_wp, &
    &   0.007492657772969_wp,   0.010735993410871_wp,   0.000043822832430_wp, &
    &  -0.004701920582963_wp,  -0.008992393110043_wp,  -0.000016900999295_wp, &
    &  -0.065308643658280_wp,   0.014998387461506_wp,  -0.000025263603064_wp, &
    &   0.026348051443162_wp,   0.000695960584767_wp,   0.000036422163237_wp, &
    &   0.011552078545256_wp,   0.013928105520511_wp,  -0.000015015405495_wp, &
    &  -0.016301123411167_wp,   0.004407346436043_wp,  -0.000110455543716_wp, &
    &  -0.007000882328314_wp,  -0.023249527069618_wp,   0.000116829671479_wp, &
    &  -0.008894681217503_wp,   0.018490096056011_wp,  -0.000046417895620_wp, &
    &   0.031760376572335_wp,  -0.027586600027256_wp,  -0.000037169135464_wp, &
    &   0.027470215463378_wp,  -0.007606885366080_wp,  -0.000031542739820_wp, &
    &   0.010038999410521_wp,   0.003472278445810_wp,  -0.000001309295826_wp, &
    &  -0.007885663618776_wp,  -0.000845723166647_wp,  -0.000014835557759_wp, &
    &  -0.001742861631250_wp,   0.000967263497864_wp,   0.000003516663274_wp, &
    &   0.000834212167689_wp,   0.002093154209008_wp,   0.000668398641072_wp, &
    &   0.000835357962140_wp,   0.002088629391442_wp,  -0.000667079195056_wp, &
    &  -0.002684468315013_wp,   0.007905570755284_wp,   0.000020202284591_wp, &
    &   0.004694468163203_wp,   0.004647604594929_wp,   0.000001966682944_wp, &
    &   0.001739391942329_wp,  -0.000902814088036_wp,   0.000582413882469_wp, &
    &   0.001758162748149_wp,  -0.000892066810248_wp,  -0.000586310681465_wp, &
    &   0.001826353796634_wp,  -0.007288556920403_wp,   0.000080112644597_wp, &
    &  -0.000722574033665_wp,   0.000176063546670_wp,  -0.001079601340131_wp, &
    &  -0.000715774633213_wp,   0.000186935581150_wp,   0.001135773146697_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn0')
    sett%chrg = 1
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
!    write (*,'(F25.15)') energy
!    write (*,'(3(F20.15,"_wp,")," &")') grad
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
  end subroutine test_gfn0_sp_cation

!========================================================================================!

  subroutine test_gfn0_sp_anion(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -42.344166332135032_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &  -0.010426382020451_wp,   0.001628829595313_wp,   0.000007210312846_wp, &
    &   0.033181572393851_wp,  -0.019170218624431_wp,  -0.000019048355272_wp, &
    &   0.023154208212837_wp,  -0.027018852959561_wp,   0.000030204155402_wp, &
    &  -0.033913203181934_wp,   0.011899345317063_wp,  -0.000034930531440_wp, &
    &  -0.026343002344950_wp,  -0.024284354151725_wp,  -0.000034684043749_wp, &
    &   0.001672756331309_wp,   0.052826607173373_wp,  -0.000034110755900_wp, &
    &  -0.036523532541808_wp,  -0.051427567776587_wp,  -0.000040815831358_wp, &
    &   0.056448934120884_wp,   0.034347894252949_wp,   0.000037020902616_wp, &
    &  -0.034867985005472_wp,  -0.008497410919343_wp,   0.000082164597833_wp, &
    &   0.004376880255650_wp,   0.008002069192913_wp,  -0.000057248113415_wp, &
    &   0.018129870164541_wp,  -0.014929003125888_wp,  -0.000042026145794_wp, &
    &  -0.004008388700057_wp,   0.029821884580246_wp,  -0.000000204657431_wp, &
    &   0.012333625323093_wp,   0.002013273343179_wp,  -0.000015895752098_wp, &
    &  -0.009301577149329_wp,  -0.006375597665755_wp,   0.000002449135366_wp, &
    &  -0.000738347066971_wp,   0.002524191851042_wp,   0.000004389653171_wp, &
    &   0.003476207329027_wp,   0.000273988079977_wp,   0.001398873626651_wp, &
    &   0.003481545977019_wp,   0.000270469106296_wp,  -0.001403512318146_wp, &
    &  -0.002469087328290_wp,   0.009366488704670_wp,   0.000032438663601_wp, &
    &   0.004640209194904_wp,   0.000420047295936_wp,   0.000003914960794_wp, &
    &  -0.001789731838366_wp,  -0.000233566591710_wp,   0.001663848929359_wp, &
    &  -0.001768091322568_wp,  -0.000209591676999_wp,  -0.001667189552238_wp, &
    &  -0.001887028684596_wp,  -0.005039986652338_wp,   0.000064725306285_wp, &
    &   0.001567359989035_wp,   0.001904232833389_wp,   0.000932148271364_wp, &
    &   0.001573187892642_wp,   0.001886828817991_wp,  -0.000909722458446_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn0')
    sett%chrg = -1
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
!    write (*,'(F25.15)') energy
!    write (*,'(3(F20.15,"_wp,")," &")') grad
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
  end subroutine test_gfn0_sp_anion

!========================================================================================!

  subroutine test_gfn0_sp_uhf(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -40.794220208842034_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &  -0.009683192083122_wp,   0.000062096462778_wp,   0.000008503079162_wp, &
    &   0.029950134754847_wp,   0.016154999834677_wp,  -0.000029559358726_wp, &
    &   0.001026094788509_wp,  -0.034289045648831_wp,   0.000026746432153_wp, &
    &  -0.023692035047363_wp,   0.004808003969286_wp,  -0.000027757841729_wp, &
    &  -0.066059837811864_wp,  -0.010598998081046_wp,  -0.000020112909173_wp, &
    &   0.033469390804893_wp,   0.034069110547216_wp,  -0.000024058369133_wp, &
    &  -0.027158226668960_wp,  -0.040298394690399_wp,  -0.000028307867095_wp, &
    &   0.045267701500952_wp,   0.028495521163272_wp,   0.000015296720032_wp, &
    &  -0.034731876264284_wp,  -0.008324613000098_wp,   0.000055864206958_wp, &
    &  -0.000950977413604_wp,   0.024709825590934_wp,  -0.000048981080232_wp, &
    &   0.028639632796385_wp,  -0.025444838911568_wp,  -0.000043224546602_wp, &
    &   0.012588974124660_wp,  -0.001417380069055_wp,   0.000000603696102_wp, &
    &   0.012881941138038_wp,   0.004366453541656_wp,  -0.000011028577992_wp, &
    &  -0.008496043097398_wp,  -0.003373490964412_wp,  -0.000015764466586_wp, &
    &  -0.001893487311906_wp,   0.002489088628558_wp,   0.000004520971875_wp, &
    &   0.001977504303537_wp,   0.001276443128655_wp,   0.001493591257119_wp, &
    &   0.001981343921494_wp,   0.001271441117791_wp,  -0.001496582541107_wp, &
    &  -0.002698130179752_wp,   0.009632226603761_wp,   0.000025851813993_wp, &
    &   0.005515249124662_wp,   0.002334826380529_wp,   0.000000147972166_wp, &
    &   0.000780537496596_wp,  -0.000202638555834_wp,   0.002005623023644_wp, &
    &   0.000801617743338_wp,  -0.000186719577198_wp,  -0.002010886741046_wp, &
    &  -0.000586042185047_wp,  -0.006771411480583_wp,   0.000076716303119_wp, &
    &   0.000532419821636_wp,   0.000625485053760_wp,   0.000482264816859_wp, &
    &   0.000537305743752_wp,   0.000612008956153_wp,  -0.000439465993761_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn0')
    sett%uhf = 2
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
!    write (*,'(F25.15)') energy
!    write (*,'(3(F20.15,"_wp,")," &")') grad
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
  end subroutine test_gfn0_sp_uhf

!========================================================================================!

  subroutine test_gfn0_sp_alpb(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -40.914343794671524_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &  -0.008898560311789_wp,   0.000892429056758_wp,   0.000004556059972_wp, &
    &  -0.000468776166310_wp,  -0.043191973189122_wp,  -0.000054862901428_wp, &
    &   0.024663000326671_wp,   0.019360907908987_wp,   0.000047913413987_wp, &
    &  -0.010385602759747_wp,   0.000142625634889_wp,  -0.000026453191079_wp, &
    &  -0.025619391962604_wp,  -0.002551430389827_wp,  -0.000033834302670_wp, &
    &  -0.005774090861161_wp,   0.019454186312275_wp,   0.000023090583090_wp, &
    &   0.007253376564080_wp,   0.005876146709966_wp,  -0.000007350539307_wp, &
    &  -0.010672008807614_wp,   0.007463393391104_wp,  -0.000107883944900_wp, &
    &  -0.007498735532863_wp,  -0.026202194147258_wp,   0.000176804338630_wp, &
    &  -0.009300643359264_wp,   0.006178344563789_wp,  -0.000051191689584_wp, &
    &   0.026956312467136_wp,  -0.021882994536491_wp,  -0.000046452481884_wp, &
    &   0.013312947174854_wp,   0.024105982533520_wp,  -0.000039072228588_wp, &
    &   0.010195763973726_wp,   0.002279939496888_wp,  -0.000004593954548_wp, &
    &  -0.008623042956319_wp,  -0.005290515866303_wp,  -0.000002934904397_wp, &
    &  -0.001195830395235_wp,   0.000342647459422_wp,   0.000002968426401_wp, &
    &   0.001562132758955_wp,   0.002094592384456_wp,   0.000349004075247_wp, &
    &   0.001563163374063_wp,   0.002090581886295_wp,  -0.000349436798086_wp, &
    &  -0.001603073760636_wp,   0.011233732948092_wp,   0.000027699940260_wp, &
    &   0.004037613498852_wp,   0.004181512677995_wp,   0.000002534288046_wp, &
    &  -0.000126245919963_wp,  -0.001310141049523_wp,  -0.000127946288841_wp, &
    &  -0.000105213385149_wp,  -0.001291320273641_wp,   0.000123231121617_wp, &
    &   0.001658673162282_wp,  -0.006431790301682_wp,   0.000072326100654_wp, &
    &  -0.000461422151234_wp,   0.001213582872177_wp,  -0.000877981885269_wp, &
    &  -0.000470344970731_wp,   0.001241753917236_wp,   0.000899866762676_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn0')
    sett%solvmodel = 'alpb'
    sett%solvent = 'water'
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
!    write (*,'(F25.15)') energy
!    write (*,'(3(F20.15,"_wp,")," &")') grad
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
  end subroutine test_gfn0_sp_alpb

!========================================================================================!
!========================================================================================!
end module test_gfn0
