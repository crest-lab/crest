module test_gfn0occ
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use testmol
  implicit none
  private

  public :: collect_gfn0occ

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
  real(wp),parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for using gfn0 in crest
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_gfn0occ(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
#ifdef WITH_GFN0
    new_unittest("Compiled gfn0 subproject      ",test_compiled_gfn0), &
    new_unittest("GFN0*-xTB singlepoint (S0)    ",test_gfn0_sp_s0), &
    new_unittest("GFN0*-xTB singlepoint (S1)    ",test_gfn0_sp_s1), &
    new_unittest("GFN0*-xTB singlepoint (S2)    ",test_gfn0_sp_s2) &
#else
    new_unittest("Compiled gfn0 subproject",test_compiled_gfn0,should_fail=.true.) &
#endif
    ]
!&>
  end subroutine collect_gfn0occ

  subroutine test_compiled_gfn0(error)
    type(error_type),allocatable,intent(out) :: error
#ifndef WITH_GFN0
    write (*,'("       ...")') 'gfn0 not compiled, expecting fail.'
    allocate (error)
#endif
  end subroutine test_compiled_gfn0

  subroutine test_gfn0_sp_s0(error)
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
    sett%id = jobtype%gfn0occ
    allocate(sett%config(2))
    sett%config = (/ 2, 0 /)
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
  end subroutine test_gfn0_sp_s0

!=======================================================================================!

  subroutine test_gfn0_sp_s1(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -40.794220064299274_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &  -0.009683192083122_wp,   0.000062096462778_wp,   0.000008503079162_wp, &
    &   0.029950134754846_wp,   0.016154999834677_wp,  -0.000029559358726_wp, &
    &   0.001026094788509_wp,  -0.034289045648831_wp,   0.000026746432153_wp, &
    &  -0.023692035047362_wp,   0.004808003969285_wp,  -0.000027757841729_wp, &
    &  -0.066059837811864_wp,  -0.010598998081046_wp,  -0.000020112909173_wp, &
    &   0.033469390804894_wp,   0.034069110547217_wp,  -0.000024058369133_wp, &
    &  -0.027158226668960_wp,  -0.040298394690399_wp,  -0.000028307867095_wp, &
    &   0.045267701500953_wp,   0.028495521163272_wp,   0.000015296720032_wp, &
    &  -0.034731876264285_wp,  -0.008324613000097_wp,   0.000055864206958_wp, &
    &  -0.000950977413603_wp,   0.024709825590934_wp,  -0.000048981080232_wp, &
    &   0.028639632796385_wp,  -0.025444838911568_wp,  -0.000043224546602_wp, &
    &   0.012588974124660_wp,  -0.001417380069055_wp,   0.000000603696102_wp, &
    &   0.012881941138038_wp,   0.004366453541656_wp,  -0.000011028577992_wp, &
    &  -0.008496043097398_wp,  -0.003373490964412_wp,  -0.000015764466586_wp, &
    &  -0.001893487311906_wp,   0.002489088628558_wp,   0.000004520971875_wp, &
    &   0.001977504303537_wp,   0.001276443128655_wp,   0.001493591257119_wp, &
    &   0.001981343921494_wp,   0.001271441117791_wp,  -0.001496582541107_wp, &
    &  -0.002698130179752_wp,   0.009632226603760_wp,   0.000025851813993_wp, &
    &   0.005515249124662_wp,   0.002334826380529_wp,   0.000000147972166_wp, &
    &   0.000780537496596_wp,  -0.000202638555834_wp,   0.002005623023643_wp, &
    &   0.000801617743338_wp,  -0.000186719577198_wp,  -0.002010886741046_wp, &
    &  -0.000586042185047_wp,  -0.006771411480583_wp,   0.000076716303119_wp, &
    &   0.000532419821636_wp,   0.000625485053760_wp,   0.000482264816859_wp, &
    &   0.000537305743752_wp,   0.000612008956153_wp,  -0.000439465993761_wp &
    & ], shape(g_ref))
!&>

    !> setup
    sett%id = jobtype%gfn0occ
    allocate(sett%config(2))
    sett%config = (/ 1, 1 /)
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
  end subroutine test_gfn0_sp_s1

!========================================================================================!

  subroutine test_gfn0_sp_s2(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -40.679590034280530_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &  -0.011090400703265_wp,  -0.000550797248329_wp,   0.000012158644771_wp, &
    &   0.059231969368526_wp,   0.074358105923631_wp,  -0.000005484741059_wp, &
    &  -0.027818712379222_wp,  -0.086619905186659_wp,   0.000006384034335_wp, &
    &  -0.032986301433234_wp,   0.010398875688673_wp,  -0.000030398873745_wp, &
    &  -0.106468035839301_wp,  -0.022368853570731_wp,  -0.000000346685034_wp, &
    &   0.072353687766465_wp,   0.048659928733743_wp,  -0.000073966588130_wp, &
    &  -0.056488580529763_wp,  -0.083360755297732_wp,  -0.000028853791246_wp, &
    &   0.096875046218719_wp,   0.047377108381327_wp,   0.000125738042617_wp, &
    &  -0.062268657067894_wp,   0.006746080150362_wp,  -0.000030672823990_wp, &
    &   0.001608123526999_wp,   0.047779250424054_wp,  -0.000043128177393_wp, &
    &   0.034770487416331_wp,  -0.032773029038149_wp,  -0.000047532282769_wp, &
    &   0.014332308826253_wp,  -0.026562470871541_wp,   0.000034831823648_wp, &
    &   0.016623917102329_wp,   0.007662792461814_wp,  -0.000016916559234_wp, &
    &  -0.008471949368056_wp,  -0.003242862418837_wp,  -0.000033223734021_wp, &
    &  -0.001876014175256_wp,   0.002488823075292_wp,   0.000004856576058_wp, &
    &   0.002711853392059_wp,   0.001618417904629_wp,   0.001138609703103_wp, &
    &   0.002717046952465_wp,   0.001611198450635_wp,  -0.001141465289808_wp, &
    &  -0.002692066285949_wp,   0.009681054218680_wp,   0.000024539819171_wp, &
    &   0.005529080674439_wp,   0.002335495633181_wp,  -0.000002635712312_wp, &
    &   0.001460375790183_wp,  -0.000269596453729_wp,   0.002500333777585_wp, &
    &   0.001483479833311_wp,  -0.000256822352583_wp,  -0.002507829221423_wp, &
    &  -0.000585237202597_wp,  -0.006742898208140_wp,   0.000076429068566_wp, &
    &   0.000521621763150_wp,   0.001024811502812_wp,   0.000267958165803_wp, &
    &   0.000526956353308_wp,   0.001006048097595_wp,  -0.000229385175492_wp &
    & ], shape(g_ref))
!&>

    !> setup
    sett%id = jobtype%gfn0occ
    allocate(sett%config(2))
    sett%config = (/ 0, 2 /)
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
  end subroutine test_gfn0_sp_s2

!========================================================================================!
!========================================================================================!
end module test_gfn0occ
