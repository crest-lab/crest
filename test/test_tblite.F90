module test_tblite
  use testdrive,only:new_unittest,unittest_type,error_type,check,test_failed
  use crest_parameters
  use crest_calculator
  use strucrd
  use testmol
  implicit none
  private

  public :: collect_tblite

  real(wp),parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
  real(wp),parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for using tblite in crest
!========================================================================================!
!========================================================================================!

!> Collect all exported unit tests
  subroutine collect_tblite(testsuite)
    !> Collection of tests
    type(unittest_type),allocatable,intent(out) :: testsuite(:)

!&<
    testsuite = [ &
#ifdef WITH_TBLITE
    new_unittest("Compiled tblite subproject    ",test_compiled_tblite), &
    new_unittest("GFN2-xTB singlepoint          ",test_gfn2_sp), &
    new_unittest("GFN2-xTB singlepoint (cation) ",test_gfn2_sp_cation), &
    new_unittest("GFN2-xTB singlepoint (anion)  ",test_gfn2_sp_anion), &
    new_unittest("GFN2-xTB singlepoint (S1)     ",test_gfn2_sp_uhf), &
    new_unittest("GFN2-xTB singlepoint (ALPB)   ",test_gfn2_sp_alpb), &
    new_unittest("GFN1-xTB singlepoint          ",test_gfn1_sp) &
#else
    new_unittest("Compiled tblite subproject",test_compiled_tblite,should_fail=.true.) &
#endif
    ]
!&>
  end subroutine collect_tblite

  subroutine test_compiled_tblite(error)
    type(error_type),allocatable,intent(out) :: error
#ifndef WITH_TBLITE
    write(*,'("       ...")') 'tblite not compiled, expecting fail.'
    allocate (error)
#endif
  end subroutine test_compiled_tblite

  subroutine test_gfn2_sp(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -42.147463158380759_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &   -0.004487752242140_wp, -0.000711346240072_wp,  0.000004422164841_wp, &
    &    0.005634362546189_wp, -0.026328061332228_wp, -0.000054754510167_wp, &
    &    0.004946048801527_wp,  0.019170217873591_wp,  0.000052031097841_wp, &
    &    0.004141916834604_wp,  0.003789795457963_wp, -0.000032613472094_wp, &
    &   -0.034492420000772_wp, -0.008306534022261_wp, -0.000038549844651_wp, &
    &    0.006099034314247_wp, -0.004026436591110_wp,  0.000034614598904_wp, &
    &    0.017464391491536_wp,  0.007913720841403_wp,  0.000037512385299_wp, &
    &   -0.014424186416761_wp,  0.007079404110621_wp, -0.000112173010722_wp, &
    &   -0.000507487343595_wp, -0.011316409953097_wp,  0.000072913115165_wp, &
    &   -0.015572958723902_wp,  0.012695690090090_wp, -0.000028265341194_wp, &
    &    0.028409899406550_wp, -0.023838232143263_wp, -0.000039688080918_wp, &
    &    0.002527690133545_wp,  0.013657975721102_wp, -0.000010799050903_wp, &
    &   -0.001619011135134_wp, -0.002968794206681_wp,  0.000002893676649_wp, &
    &   -0.006515329457992_wp,  0.007906403753381_wp, -0.000058364122249_wp, &
    &   -0.001453539200851_wp,  0.002783638728611_wp,  0.000004398794037_wp, &
    &    0.002596908343625_wp,  0.000907279408505_wp,  0.003981683485842_wp, &
    &    0.002598744285592_wp,  0.000908311001962_wp, -0.003984457648426_wp, &
    &   -0.003774165613087_wp,  0.008367958950496_wp,  0.000028978030648_wp, &
    &    0.007868100320925_wp,  0.002139539119133_wp,  0.000000731417166_wp, &
    &    0.000932101707044_wp, -0.000165802818748_wp,  0.003249098830930_wp, &
    &    0.000957167451635_wp, -0.000141981017010_wp, -0.003253611423897_wp, &
    &   -0.002069570321846_wp, -0.009289127414772_wp,  0.000103587670973_wp, &
    &    0.000358684569710_wp, -0.000098291491328_wp,  0.002383904434443_wp, &
    &    0.000381370249352_wp, -0.000128917826287_wp, -0.002343493197516_wp  &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn2')
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
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
  end subroutine test_gfn2_sp

!=======================================================================================!

  subroutine test_gfn2_sp_cation(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -41.662902644562010_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &   -0.008680441178755_wp, -0.003805585932997_wp,  0.000003833547822_wp, & 
    &    0.008556211724477_wp,  0.007298329032302_wp, -0.000060715808107_wp, &
    &   -0.017210129653521_wp,  0.008000848866626_wp,  0.000045675030514_wp, &
    &    0.014180937675588_wp, -0.003569605927023_wp, -0.000026595098533_wp, &
    &   -0.065969546794406_wp, -0.001431750313422_wp, -0.000024664910677_wp, &
    &    0.035369189836543_wp, -0.022587142672832_wp,  0.000057299230498_wp, &
    &    0.027217881045369_wp,  0.025975982754290_wp,  0.000020410444564_wp, &
    &   -0.024543524792713_wp, -0.000639673261627_wp, -0.000111734218168_wp, &
    &   -0.002428707956668_wp, -0.010444903854503_wp,  0.000058564290422_wp, &
    &   -0.016364194518134_wp,  0.044172268722857_wp, -0.000063262605340_wp, &
    &    0.037221292516429_wp, -0.035587227765337_wp, -0.000030296425373_wp, &
    &    0.004976226192841_wp, -0.023995007612999_wp,  0.000028348988032_wp, &
    &    0.005579294937249_wp,  0.005081360488652_wp,  0.000017005851806_wp, &
    &   -0.004901080903271_wp,  0.007389912260113_wp, -0.000058731488109_wp, &
    &   -0.001138968489320_wp,  0.003005596945492_wp,  0.000006328264437_wp, &
    &    0.003099117352777_wp,  0.001897847123878_wp,  0.003761878500882_wp, &
    &    0.003099436103951_wp,  0.001894452000258_wp, -0.003760475465119_wp, &
    &   -0.004253493889487_wp,  0.006389910199458_wp,  0.000023261411435_wp, &
    &    0.006392948917928_wp,  0.002861653171794_wp, -0.000007274304368_wp, &
    &    0.001014741765200_wp, -0.002193018249538_wp,  0.003407479151927_wp, &
    &    0.001044587604884_wp, -0.002181552984101_wp, -0.003420727645756_wp, &
    &   -0.002157293947451_wp, -0.008075804379714_wp,  0.000093116929524_wp, &
    &   -0.000060141992576_wp,  0.000289154141166_wp,  0.002614714491853_wp, &
    &   -0.000044341556934_wp,  0.000253957247207_wp, -0.002573438164165_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn2')
    sett%chrg = 1
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
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
  end subroutine test_gfn2_sp_cation

!========================================================================================!

  subroutine test_gfn2_sp_anion(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -42.311955686578870_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &   -0.011470105970055_wp,  0.000454323102443_wp,  0.000004025088780_wp, & 
    &    0.035872697562889_wp, -0.008460078242858_wp, -0.000021111515370_wp, &
    &    0.003692743414081_wp, -0.020258524834661_wp,  0.000037721453388_wp, &
    &   -0.019087091757337_wp,  0.017138450608344_wp, -0.000045771350144_wp, &
    &   -0.020213616463304_wp, -0.039924769953061_wp, -0.000028414555700_wp, &
    &    0.009408652907122_wp,  0.033864938287135_wp, -0.000033430463173_wp, &
    &   -0.016485709515356_wp, -0.051197854678100_wp,  0.000036590001373_wp, &
    &    0.045058407382885_wp,  0.033979204637975_wp, -0.000005048366790_wp, &
    &   -0.036758692818602_wp,  0.001990415645982_wp,  0.000018141382464_wp, &
    &    0.021925167453428_wp,  0.005846953244880_wp, -0.000041044406011_wp, &
    &    0.003436690227325_wp, -0.004233961986466_wp, -0.000060952000306_wp, &
    &   -0.023154168166028_wp,  0.024946658253272_wp,  0.000030358110993_wp, &
    &    0.004916230569813_wp,  0.001138591313257_wp, -0.000009402454000_wp, &
    &   -0.008313799516116_wp,  0.002428139980713_wp, -0.000040049501514_wp, &
    &   -0.001208556228635_wp,  0.000961400750060_wp,  0.000003367872204_wp, &
    &    0.004124763613055_wp,  0.001669700587915_wp,  0.003007813202271_wp, &
    &    0.004131893302950_wp,  0.001670334153738_wp, -0.003008973034920_wp, &
    &   -0.002632048440229_wp,  0.006629738397082_wp,  0.000032670619566_wp, &
    &    0.007424243930173_wp,  0.002858491363991_wp, -0.000000133946349_wp, &
    &    0.000371778457024_wp, -0.001850778398136_wp,  0.002001649485315_wp, &
    &    0.000408679516879_wp, -0.001820688458466_wp, -0.002012504665727_wp, &
    &   -0.001161553361651_wp, -0.010738152646510_wp,  0.000113042740553_wp, &
    &   -0.000156331938699_wp,  0.001465258504155_wp,  0.001088496060505_wp, &
    &   -0.000130274161610_wp,  0.001442210367316_wp, -0.001067039757407_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn2')
    sett%chrg = -1
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
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
  end subroutine test_gfn2_sp_anion

!========================================================================================!

  subroutine test_gfn2_sp_uhf(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -42.015481983629080_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &   -0.008640563198637_wp, -0.001276731011454_wp,  0.000007958615032_wp, &  
    &    0.037512683660300_wp,  0.023726123308498_wp, -0.000011120132727_wp, &
    &   -0.023304344662829_wp, -0.030936043041433_wp,  0.000027310525144_wp, &
    &   -0.003125727457313_wp,  0.011153715001891_wp, -0.000039680963846_wp, &
    &   -0.071897035549331_wp, -0.035684240793218_wp, -0.000021364112392_wp, &
    &    0.053937213287434_wp,  0.027607622002516_wp, -0.000019246037097_wp, &
    &   -0.015752681349571_wp, -0.038841382300575_wp,  0.000004645104546_wp, &
    &    0.035638413983781_wp,  0.024134341589115_wp,  0.000001459196707_wp, &
    &   -0.037275602976986_wp, -0.000596375689797_wp,  0.000022849188633_wp, &
    &    0.011882587422559_wp,  0.038230111262650_wp, -0.000064502857749_wp, &
    &    0.022686028312851_wp, -0.023788517748961_wp, -0.000048189366831_wp, &
    &   -0.006326790373970_wp, -0.007705697383338_wp,  0.000050970659367_wp, &
    &    0.002507168942608_wp,  0.000804113685594_wp, -0.000002670812276_wp, &
    &   -0.007100255305176_wp,  0.007182777740523_wp, -0.000064268848580_wp, &
    &   -0.002073508104906_wp,  0.002718196006540_wp,  0.000005207667148_wp, &
    &    0.003175843964289_wp,  0.001471376588308_wp,  0.003688753934759_wp, &
    &    0.003180536309837_wp,  0.001468559650412_wp, -0.003691578785923_wp, &
    &   -0.004351450534000_wp,  0.008596740689697_wp,  0.000025895635956_wp, &
    &    0.007700374408025_wp,  0.002614064547437_wp, -0.000002553113761_wp, &
    &    0.001391017646603_wp, -0.000983369331829_wp,  0.003331940615500_wp, &
    &    0.001421163490921_wp, -0.000962246268976_wp, -0.003341342331401_wp, &
    &   -0.002160782577365_wp, -0.009722835050767_wp,  0.000104945398812_wp, &
    &    0.000477407267670_wp,  0.000411637518681_wp,  0.002100726847295_wp, &
    &    0.000498303393207_wp,  0.000378059028485_wp, -0.002066146026319_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn2')
    sett%uhf = 2
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
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
  end subroutine test_gfn2_sp_uhf

!========================================================================================!

  subroutine test_gfn2_sp_alpb(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -42.150818113966622_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &   -0.004041852805369_wp, -0.001459116937018_wp,  0.000004212952013_wp, & 
    &    0.005431782185647_wp, -0.026769409144275_wp, -0.000053933849249_wp, &
    &    0.005232547908947_wp,  0.019316152153876_wp,  0.000051828956087_wp, &
    &    0.003984299183311_wp,  0.004654018323578_wp, -0.000031672641741_wp, &
    &   -0.034791239964625_wp, -0.009664440508457_wp, -0.000039268100546_wp, &
    &    0.006792385075909_wp, -0.003127178836761_wp,  0.000032791144117_wp, &
    &    0.015306198183764_wp,  0.005195424722641_wp,  0.000035383591075_wp, &
    &   -0.011644038936970_wp,  0.009526167376125_wp, -0.000103656945404_wp, &
    &   -0.001073718470205_wp, -0.010467260905049_wp,  0.000071170633309_wp, &
    &   -0.012208257598543_wp,  0.007792225857756_wp, -0.000026835529946_wp, &
    &    0.024297093142032_wp, -0.020057317289964_wp, -0.000036190038701_wp, &
    &    0.002288207714506_wp,  0.015784323612706_wp, -0.000015109378387_wp, &
    &   -0.001839503827186_wp, -0.003557956433963_wp,  0.000003971003294_wp, &
    &   -0.006475695156935_wp,  0.007869132396390_wp, -0.000061014145496_wp, &
    &   -0.001323378271276_wp,  0.002797437673167_wp,  0.000004384755727_wp, &
    &    0.002794743156382_wp,  0.000752555767144_wp,  0.003757718239733_wp, &
    &    0.002796550406646_wp,  0.000753562108628_wp, -0.003760695667456_wp, &
    &   -0.003950802127381_wp,  0.008247313874474_wp,  0.000028580331243_wp, &
    &    0.007699255283409_wp,  0.001924889526582_wp,  0.000000878165662_wp, &
    &    0.000911988480912_wp, -0.000170582007313_wp,  0.003208806647171_wp, &
    &    0.000936315963916_wp, -0.000147196304551_wp, -0.003212846696159_wp, &
    &   -0.002208060498652_wp, -0.009259192346010_wp,  0.000102993956974_wp, &
    &    0.000531577240056_wp,  0.000047927691044_wp,  0.002360617394228_wp, &
    &    0.000553603731705_wp,  0.000018519629250_wp, -0.002322114777547_wp &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn2')
    sett%solvmodel = 'alpb'
    sett%solvent = 'water'
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
    !write(*,'(F25.15)') energy
    !write(*,'(3F25.15)') grad
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
  end subroutine test_gfn2_sp_alpb

!========================================================================================!

  subroutine test_gfn1_sp(error)
    type(error_type),allocatable,intent(out) :: error
    type(calcdata) :: calc
    type(calculation_settings) :: sett
    type(coord) :: mol
    real(wp) :: energy
    real(wp),allocatable :: grad(:,:)
    integer :: io
!&<
    real(wp),parameter :: e_ref = -44.509702451339784_wp
    real(wp),parameter :: g_ref(3,24) = reshape([&
    &    0.003874633213987_wp,  0.001058841664050_wp,  0.000004010132495_wp, &  
    &   -0.006820697815903_wp, -0.025467667208147_wp, -0.000047472704275_wp, &
    &    0.010650204538643_wp,  0.008835966788566_wp,  0.000047997051770_wp, &
    &    0.000318472512888_wp,  0.007424576243087_wp, -0.000028497241728_wp, &
    &   -0.026207800400472_wp, -0.009013375235391_wp, -0.000038669187497_wp, &
    &   -0.001750473728420_wp,  0.003668870753592_wp,  0.000034119493952_wp, &
    &    0.019485861483131_wp,  0.004742782186204_wp,  0.000040473978448_wp, &
    &   -0.011724493022376_wp,  0.007481847734935_wp, -0.000116124549323_wp, &
    &   -0.001224423088802_wp, -0.019133403500940_wp,  0.000124617291725_wp, &
    &   -0.016766618899760_wp,  0.012104298176132_wp, -0.000034576323330_wp, &
    &    0.024844625147691_wp, -0.019287644128199_wp, -0.000049197092067_wp, &
    &    0.011083811277529_wp,  0.017004904025028_wp, -0.000003221055588_wp, &
    &   -0.006459232625454_wp, -0.004286162846774_wp, -0.000004086343732_wp, &
    &   -0.006590040129048_wp,  0.013220910391538_wp, -0.000090165069414_wp, &
    &   -0.001326880708616_wp,  0.002124670840594_wp,  0.000004093024610_wp, &
    &    0.001722713105321_wp,  0.001171656543542_wp,  0.002684536818471_wp, &
    &    0.001725380313342_wp,  0.001170987522855_wp, -0.002688028482816_wp, &
    &   -0.001550167832156_wp,  0.003931581737480_wp,  0.000025224933007_wp, &
    &    0.006357182837352_wp,  0.002289813675783_wp,  0.000003480620835_wp, &
    &    0.000517412225517_wp, -0.000506087340212_wp,  0.001942472921733_wp, &
    &    0.000540625438461_wp, -0.000485737254273_wp, -0.001945673784378_wp, &
    &   -0.001267730219018_wp, -0.007658771744659_wp,  0.000084917355339_wp, &
    &    0.000277986074076_wp, -0.000182942929362_wp,  0.001400429395618_wp, &
    &    0.000289650302091_wp, -0.000209916095429_wp, -0.001350661183856_wp  &
    & ], shape(g_ref))
!&>

    !> setup
    call sett%create('gfn1')
    call calc%add(sett)
    call get_testmol('caffeine',mol)
    allocate (grad(3,mol%nat))

    !> calculation
    call engrad(mol,calc,energy,grad,io)
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
  end subroutine test_gfn1_sp

!========================================================================================!
!========================================================================================!
end module test_tblite
