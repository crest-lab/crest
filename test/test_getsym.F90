module test_getsym
  use testdrive,only:new_unittest,unittest_type,error_type,test_failed
  use iso_fortran_env,only:wp => real64
  use symmetry_i,only:getsym
  implicit none
  private
  public :: collect_getsym

!========================================================================================!
!> Hardcoded geometries (converted from XYZ/Angstrom to Bohr)
!========================================================================================!

  integer,parameter :: nat_c1 = 5
  integer,parameter :: at_c1(nat_c1) = [6,1,9,17,35]
  real(wp),parameter :: xyz_c1(3,nat_c1) = reshape([ &
    1.7686891663_wp, -0.0460904202_wp, -0.1376098564_wp, &
    3.8326480397_wp, -0.0460904202_wp, -0.1376098564_wp, &
    1.0806965762_wp, 1.4516120199_wp, -1.3799725025_wp, &
    1.0806965762_wp, 0.2809833775_wp, 1.7806322355_wp, &
    1.0806965762_wp, -1.8708666579_wp, -0.8134704049_wp], &
    [3,nat_c1])

  integer,parameter :: nat_ci = 8
  integer,parameter :: at_ci(nat_ci) = [17,6,6,17,9,1,9,1]
  real(wp),parameter :: xyz_ci(3,nat_ci) = reshape([ &
    1.9984798631_wp, 0.0892139703_wp, 0.0537060165_wp, &
    5.3641387828_wp, 0.1444128704_wp, 0.0626444210_wp, &
    6.3141419002_wp, 2.8167312751_wp, 0.4807085316_wp, &
    9.6797819227_wp, 2.8713254628_wp, 0.4934830802_wp, &
    6.0032441582_wp, -0.6077548189_wp, -1.7546107067_wp, &
    6.0030173911_wp, -1.1269192772_wp, 1.5632759366_wp, &
    5.6729956206_wp, 3.5699572111_wp, 2.2967542346_wp, &
    5.6772286072_wp, 4.0874209158_wp, -1.0213780731_wp], &
    [3,nat_ci])

  integer,parameter :: nat_cs = 5
  integer,parameter :: at_cs(nat_cs) = [6,1,1,17,35]
  real(wp),parameter :: xyz_cs(3,nat_cs) = reshape([ &
    1.7686891663_wp, -0.0460904202_wp, -0.1376098564_wp, &
    3.8326480397_wp, -0.0460904202_wp, -0.1376098564_wp, &
    1.0806965762_wp, 1.4516120199_wp, -1.3799725025_wp, &
    1.0806965762_wp, 0.2809833775_wp, 1.7806322355_wp, &
    1.0806965762_wp, -1.8708666579_wp, -0.8134704049_wp], &
    [3,nat_cs])

  integer,parameter :: nat_c2 = 4
  integer,parameter :: at_c2(nat_c2) = [8,8,1,1]
  real(wp),parameter :: xyz_c2(3,nat_c2) = reshape([ &
    -1.3413435337_wp, 0.0553519912_wp, -0.0129809434_wp, &
    1.3413435337_wp, -0.0553519912_wp, 0.0129809434_wp, &
    -1.7986927696_wp, -1.6645598353_wp, 0.3904096780_wp, &
    1.7986927696_wp, -1.6645598353_wp, -0.3904096780_wp], &
    [3,nat_c2])

  integer,parameter :: nat_c2h = 6
  integer,parameter :: at_c2h(nat_c2h) = [17,6,6,17,1,1]
  real(wp),parameter :: xyz_c2h(3,nat_c2h) = reshape([ &
    1.8977196661_wp, -0.2014448049_wp, -0.0710348050_wp, &
    5.1323638737_wp, -0.1056356904_wp, -0.0379268033_wp, &
    6.4183602960_wp, 1.9359677229_wp, 0.6678481097_wp, &
    9.6529856062_wp, 2.0317768374_wp, 0.7009750087_wp, &
    5.9826083489_wp, -1.8636101096_wp, -0.6456438277_wp, &
    5.5680969235_wp, 3.6939421421_wp, 1.2755840314_wp], &
    [3,nat_c2h])

  integer,parameter :: nat_c2v = 3
  integer,parameter :: at_c2v(nat_c2v) = [8,1,1]
  real(wp),parameter :: xyz_c2v(3,nat_c2v) = reshape([ &
    1.7412503430_wp, -0.0503045094_wp, 0.1594739877_wp, &
    3.5702973618_wp, 0.0293096522_wp, 0.1225676364_wp, &
    1.2143191104_wp, 1.5408070930_wp, -0.5779349407_wp], &
    [3,nat_c2v])

  integer,parameter :: nat_c3 = 8
  integer,parameter :: at_c3(nat_c3) = [8,15,8,8,8,1,1,1]
  real(wp),parameter :: xyz_c3(3,nat_c3) = reshape([ &
    2.1838664179_wp, 0.0466527838_wp, 1.9125045432_wp, &
    -0.2066665717_wp, -0.0189806738_wp, 0.1004614502_wp, &
    -2.6626254324_wp, -0.0861740813_wp, 1.3747444618_wp, &
    0.2477461752_wp, 2.3413303983_wp, -1.6955100891_wp, &
    0.3753052019_wp, -2.3513858200_wp, -1.6952502706_wp, &
    1.8620047872_wp, -0.9714351061_wp, 3.3672205705_wp, &
    -0.5208961465_wp, 3.8036152890_wp, -0.9694469441_wp, &
    -1.1456946322_wp, -2.8663836211_wp, -2.5184443558_wp], &
    [3,nat_c3])

  integer,parameter :: nat_c3v = 4
  integer,parameter :: at_c3v(nat_c3v) = [7,1,1,1]
  real(wp),parameter :: xyz_c3v(3,nat_c3v) = reshape([ &
    2.0094213774_wp, 0.1277832805_wp, -0.0094864251_wp, &
    3.9318964557_wp, 0.0784425314_wp, -0.1081112316_wp, &
    1.4206016142_wp, -1.6941772652_wp, -0.2138414083_wp, &
    1.4205827169_wp, 1.0563191091_wp, -1.5903557120_wp], &
    [3,nat_c3v])

  integer,parameter :: nat_c4v = 6
  integer,parameter :: at_c4v(nat_c4v) = [54,8,9,9,9,9]
  real(wp),parameter :: xyz_c4v(3,nat_c4v) = reshape([ &
    0.0000000000_wp, 0.0000000000_wp, 0.0000000000_wp, &
    2.8345891869_wp, 0.0000000000_wp, 0.0000000000_wp, &
    0.0000000000_wp, 2.2676713496_wp, 2.2676713496_wp, &
    0.0000000000_wp, 2.2676713496_wp, -2.2676713496_wp, &
    0.0000000000_wp, -2.2676713496_wp, 2.2676713496_wp, &
    0.0000000000_wp, -2.2676713496_wp, -2.2676713496_wp], &
    [3,nat_c4v])

  integer,parameter :: nat_c5 = 30
  integer,parameter :: at_c5(nat_c5) = [6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,1,9,1,9,1,9,9,1,1,9]
  real(wp),parameter :: xyz_c5(3,nat_c5) = reshape([ &
    4.1506027534_wp, -2.0577782147_wp, -0.7793164178_wp, &
    6.1034755769_wp, -0.2721278465_wp, -0.2078341373_wp, &
    5.6677527186_wp, 2.2879313702_wp, -0.1333417100_wp, &
    3.2287955057_wp, 3.3565846641_wp, -0.6214740288_wp, &
    2.1555044741_wp, 5.7218229465_wp, 0.1429213557_wp, &
    -0.4097575886_wp, 6.0946420523_wp, 0.3095743804_wp, &
    -2.1972322554_wp, 4.1450305843_wp, -0.2688976403_wp, &
    -4.7566213713_wp, 3.8026161128_wp, 0.5439005762_wp, &
    -5.9057076660_wp, 1.4728953771_wp, 0.5727859198_wp, &
    -4.6279273429_wp, -0.7822543594_wp, -0.2081874783_wp, &
    -2.3260847824_wp, -0.3558899934_wp, -1.3645953292_wp, &
    -0.4212128650_wp, -2.2172109020_wp, -1.5190221584_wp, &
    1.9366611069_wp, -0.9752564488_wp, -1.6417754489_wp, &
    1.4889571281_wp, 1.6533686360_wp, -1.5646128982_wp, &
    -1.1457234397_wp, 2.0362307059_wp, -1.3934152450_wp, &
    -0.7053856104_wp, -4.6162955879_wp, -0.5285699878_wp, &
    -3.2269454839_wp, -5.1877810615_wp, 0.2742834659_wp, &
    -5.0813249127_wp, -3.3749399931_wp, 0.4262836582_wp, &
    1.6314970454_wp, -5.8938226986_wp, -0.0423215999_wp, &
    3.9271679611_wp, -4.6839996881_wp, -0.1604477348_wp, &
    7.9404810111_wp, -0.9711247579_wp, 0.3544797089_wp, &
    7.1758252770_wp, 3.5221796708_wp, 0.4851426936_wp, &
    3.4066392826_wp, 7.2100758904_wp, 0.7757496641_wp, &
    -1.0969646698_wp, 7.8643592820_wp, 1.0684833919_wp, &
    -5.7621122661_wp, 5.4024319347_wp, 1.3250328734_wp, &
    -7.7790778971_wp, 1.3128978782_wp, 1.3759015475_wp, &
    -3.6400790044_wp, -7.0730048337_wp, 0.9489170013_wp, &
    -6.8954493888_wp, -3.8908471871_wp, 1.2160970832_wp, &
    1.5754712934_wp, -7.8309868158_wp, 0.6087556096_wp, &
    5.6057219651_wp, -5.7074640447_wp, 0.4016119744_wp], &
    [3,nat_c5])

  integer,parameter :: nat_c5v = 30
  integer,parameter :: at_c5v(nat_c5v) = [6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,1,1,1,1,1,1,1,1,1,1]
  real(wp),parameter :: xyz_c5v(3,nat_c5v) = reshape([ &
    4.1506027534_wp, -2.0577782147_wp, -0.7793164178_wp, &
    6.1034755769_wp, -0.2721278465_wp, -0.2078341373_wp, &
    5.6677527186_wp, 2.2879313702_wp, -0.1333417100_wp, &
    3.2287955057_wp, 3.3565846641_wp, -0.6214740288_wp, &
    2.1555044741_wp, 5.7218229465_wp, 0.1429213557_wp, &
    -0.4097575886_wp, 6.0946420523_wp, 0.3095743804_wp, &
    -2.1972322554_wp, 4.1450305843_wp, -0.2688976403_wp, &
    -4.7566213713_wp, 3.8026161128_wp, 0.5439005762_wp, &
    -5.9057076660_wp, 1.4728953771_wp, 0.5727859198_wp, &
    -4.6279273429_wp, -0.7822543594_wp, -0.2081874783_wp, &
    -2.3260847824_wp, -0.3558899934_wp, -1.3645953292_wp, &
    -0.4212128650_wp, -2.2172109020_wp, -1.5190221584_wp, &
    1.9366611069_wp, -0.9752564488_wp, -1.6417754489_wp, &
    1.4889571281_wp, 1.6533686360_wp, -1.5646128982_wp, &
    -1.1457234397_wp, 2.0362307059_wp, -1.3934152450_wp, &
    -0.7053856104_wp, -4.6162955879_wp, -0.5285699878_wp, &
    -3.2269454839_wp, -5.1877810615_wp, 0.2742834659_wp, &
    -5.0813249127_wp, -3.3749399931_wp, 0.4262836582_wp, &
    1.6314970454_wp, -5.8938226986_wp, -0.0423215999_wp, &
    3.9271679611_wp, -4.6839996881_wp, -0.1604477348_wp, &
    7.9404810111_wp, -0.9711247579_wp, 0.3544797089_wp, &
    7.1758252770_wp, 3.5221796708_wp, 0.4851426936_wp, &
    3.4066392826_wp, 7.2100758904_wp, 0.7757496641_wp, &
    -1.0969646698_wp, 7.8643592820_wp, 1.0684833919_wp, &
    -5.7621122661_wp, 5.4024319347_wp, 1.3250328734_wp, &
    -7.7790778971_wp, 1.3128978782_wp, 1.3759015475_wp, &
    -3.6400790044_wp, -7.0730048337_wp, 0.9489170013_wp, &
    -6.8954493888_wp, -3.8908471871_wp, 1.2160970832_wp, &
    1.5754712934_wp, -7.8309868158_wp, 0.6087556096_wp, &
    5.6057219651_wp, -5.7074640447_wp, 0.4016119744_wp], &
    [3,nat_c5v])

  integer,parameter :: nat_cinfv = 2
  integer,parameter :: at_cinfv(nat_cinfv) = [9,1]
  real(wp),parameter :: xyz_cinfv(3,nat_cinfv) = reshape([ &
    1.6545119139_wp, 0.1226999173_wp, -0.0701655310_wp, &
    3.4293237929_wp, 0.1226999173_wp, -0.0701655310_wp], &
    [3,nat_cinfv])

  integer,parameter :: nat_d2 = 22
  integer,parameter :: at_d2(nat_d2) = [6,6,6,6,6,6,6,6,6,6,6,6,1,1,1,1,1,1,1,1,1,1]
  real(wp),parameter :: xyz_d2(3,nat_d2) = reshape([ &
    2.5621851661_wp, -0.4439911530_wp, -0.0436148790_wp, &
    0.9253421914_wp, -2.5127877252_wp, 0.1502332269_wp, &
    -1.6739005039_wp, -2.1369023017_wp, 0.1913347701_wp, &
    -2.6332955601_wp, 0.3046994403_wp, 0.0348465497_wp, &
    -0.9891960372_wp, 2.3654268820_wp, -0.1896340166_wp, &
    1.6529812357_wp, 2.0499560027_wp, -0.2395794781_wp, &
    3.3960268186_wp, 4.2370305331_wp, -0.4989632859_wp, &
    2.7654063135_wp, 6.6565791713_wp, 0.4215601039_wp, &
    4.3970903358_wp, 8.7186861130_wp, 0.1505166858_wp, &
    6.7269714694_wp, 8.4139488781_wp, -1.0255921624_wp, &
    7.4141325801_wp, 6.0515266607_wp, -1.9418447711_wp, &
    5.7626631253_wp, 3.9991329112_wp, -1.6961614777_wp, &
    4.5822457015_wp, -0.8214828436_wp, -0.0238294464_wp, &
    1.6914371624_wp, -4.4146269943_wp, 0.2723851236_wp, &
    -2.9477459873_wp, -3.7402404322_wp, 0.3445726616_wp, &
    -4.6633527468_wp, 0.6148223946_wp, 0.0725843804_wp, &
    -1.8237746829_wp, 4.2357833138_wp, -0.3557031484_wp, &
    0.9789537216_wp, 6.9700847353_wp, 1.3885707564_wp, &
    3.8407738620_wp, 10.5624352010_wp, 0.8651544144_wp, &
    8.0008169527_wp, 10.0118257001_wp, -1.2267913028_wp, &
    9.2341467080_wp, 5.7935979419_wp, -2.8609508664_wp, &
    6.3573788339_wp, 2.1980160446_wp, -2.4864260457_wp], &
    [3,nat_d2])

  integer,parameter :: nat_d2d = 7
  integer,parameter :: at_d2d(nat_d2d) = [6,6,6,1,1,1,1]
  real(wp),parameter :: xyz_d2d(3,nat_d2d) = reshape([ &
    -2.4443487644_wp, -0.0000022937_wp, -0.0000001336_wp, &
    -0.0000010639_wp, -0.0000043018_wp, -0.0000001767_wp, &
    2.4443458222_wp, -0.0000057331_wp, 0.0000012450_wp, &
    -3.5163354416_wp, -0.4079430602_wp, -1.6837648529_wp, &
    -3.5163315932_wp, 0.4079495522_wp, 1.6837640156_wp, &
    3.5163309025_wp, -1.6837439760_wp, 0.4080523160_wp, &
    3.5163297587_wp, 1.6837342434_wp, -0.4080524134_wp], &
    [3,nat_d2d])

  integer,parameter :: nat_d2h = 6
  integer,parameter :: at_d2h(nat_d2h) = [6,6,1,1,1,1]
  real(wp),parameter :: xyz_d2h(3,nat_d2h) = reshape([ &
    2.0246903644_wp, 0.0223932546_wp, 0.1508190420_wp, &
    4.5488164464_wp, 0.0223932546_wp, 0.1508190420_wp, &
    0.9669161662_wp, -0.8803667097_wp, 1.6587637977_wp, &
    0.9669161662_wp, 0.9251343216_wp, -1.3571446109_wp, &
    5.6066095419_wp, -0.8803667097_wp, 1.6587637977_wp, &
    5.6066095419_wp, 0.9251343216_wp, -1.3571446109_wp], &
    [3,nat_d2h])

  integer,parameter :: nat_d3 = 19
  integer,parameter :: at_d3(nat_d3) = [8,6,6,8,8,8,26,8,6,6,8,8,8,8,6,6,8,8,8]
  real(wp),parameter :: xyz_d3(3,nat_d3) = reshape([ &
    -5.5946338021_wp, -0.8229064614_wp, 5.2053622115_wp, &
    -3.8359811652_wp, -0.1763064438_wp, 3.8074275063_wp, &
    -2.7745085917_wp, 2.5933383340_wp, 3.8848133004_wp, &
    -3.5603543255_wp, 4.1325496693_wp, 5.4589836994_wp, &
    -1.1089414828_wp, 3.0199843392_wp, 2.2154904239_wp, &
    -2.7131559285_wp, -1.5806752624_wp, 2.2233057633_wp, &
    -0.0134444680_wp, 0.0555107607_wp, -0.0126742790_wp, &
    2.6937923405_wp, -1.3396015306_wp, 2.3709864194_wp, &
    3.8108202955_wp, -3.3564374539_wp, 1.7154489832_wp, &
    2.7271524135_wp, -4.5428355546_wp, -0.7798409646_wp, &
    3.4884935170_wp, -6.6150723077_wp, -1.5481502414_wp, &
    1.0689832135_wp, -3.1742433657_wp, -1.8393468833_wp, &
    5.5813984523_wp, -4.3738859186_wp, 2.8528582093_wp, &
    2.3849000934_wp, 2.0543474787_wp, -2.2983085914_wp, &
    1.4354620793_wp, 3.1595785624_wp, -4.2017472655_wp, &
    -1.4637774680_wp, 2.6289929240_wp, -4.5453833855_wp, &
    -2.6147523355_wp, 3.4372639067_wp, -6.4129074540_wp, &
    -2.4130858761_wp, 1.3426279459_wp, -2.7590959144_wp, &
    2.5858890563_wp, 4.5360541927_wp, -5.7014434355_wp], &
    [3,nat_d3])

  integer,parameter :: nat_d3d = 8
  integer,parameter :: at_d3d(nat_d3d) = [6,6,1,1,1,1,1,1]
  real(wp),parameter :: xyz_d3d(3,nat_d3d) = reshape([ &
    2.0104607267_wp, -0.0348843443_wp, 0.1472096651_wp, &
    4.8678589080_wp, -0.0348843443_wp, 0.1472096651_wp, &
    1.2842956688_wp, 0.8372809540_wp, 1.8753642061_wp, &
    1.2842956688_wp, -1.9675828410_wp, 0.0384370294_wp, &
    1.2842956688_wp, 1.0256866487_wp, -1.4721911374_wp, &
    5.5940239659_wp, -1.0954364399_wp, 1.7665915703_wp, &
    5.5940239659_wp, 1.8978330497_wp, 0.2559823008_wp, &
    5.5940428632_wp, -0.9070307453_wp, -1.5809637731_wp], &
    [3,nat_d3d])

  integer,parameter :: nat_d3h = 4
  integer,parameter :: at_d3h(nat_d3h) = [9,5,9,9]
  real(wp),parameter :: xyz_d3h(3,nat_d3h) = reshape([ &
    1.9631230873_wp, -0.1245707461_wp, 0.0334670497_wp, &
    4.6923411455_wp, -0.1245707461_wp, 0.0334670497_wp, &
    6.0569501747_wp, -1.5824944513_wp, -1.8268927310_wp, &
    6.0569501747_wp, 1.3333718563_wp, 1.8938268303_wp], &
    [3,nat_d3h])

  integer,parameter :: nat_d4 = 22
  integer,parameter :: at_d4(nat_d4) = [8,6,25,6,8,6,8,6,8,6,8,25,6,8,6,8,6,8,6,8,6,8]
  real(wp),parameter :: xyz_d4(3,nat_d4) = reshape([ &
    2.2841686586_wp, 0.2090037094_wp, 0.1362114591_wp, &
    4.5444511790_wp, 0.1935646469_wp, 0.0922564294_wp, &
    8.0752343675_wp, 0.1854010301_wp, 0.0151933980_wp, &
    8.0120608232_wp, 1.9881619585_wp, -3.0202358814_wp, &
    7.9811071093_wp, 3.1479246757_wp, -4.9606066662_wp, &
    8.0751020867_wp, 3.2437337902_wp, 1.7781188997_wp, &
    8.0761981279_wp, 5.2013388772_wp, 2.9070034920_wp, &
    11.6050160012_wp, 0.2565681159_wp, -0.0161004666_wp, &
    13.8652040353_wp, 0.3039435499_wp, -0.0245286451_wp, &
    8.1379165831_wp, -1.5381803737_wp, 3.0963540497_wp, &
    8.1691726532_wp, -2.6353931561_wp, 5.0726107336_wp, &
    8.0741761209_wp, -4.5149147597_wp, -2.6937101043_wp, &
    9.8948705502_wp, -3.0413819168_wp, -5.3358495828_wp, &
    11.0515152193_wp, -2.1019990602_wp, -7.0370377319_wp, &
    11.1000622834_wp, -5.4598912056_wp, -1.1375395380_wp, &
    13.0420960272_wp, -6.0674948464_wp, -0.1537481175_wp, &
    8.0747430387_wp, -7.5720380951_wp, -4.4554072840_wp, &
    8.0736658949_wp, -9.5304179698_wp, -5.5823076639_wp, &
    6.2539541231_wp, -6.0621280242_wp, -0.0951666076_wp, &
    5.0967614335_wp, -7.0624734456_wp, 1.5704946904_wp, &
    5.0491781296_wp, -3.6432029957_wp, -4.2940624675_wp, &
    3.1068042352_wp, -3.0971855292_wp, -5.3130027939_wp], &
    [3,nat_d4])

  integer,parameter :: nat_d4h = 5
  integer,parameter :: at_d4h(nat_d4h) = [54,9,9,9,9]
  real(wp),parameter :: xyz_d4h(3,nat_d4h) = reshape([ &
    0.0000000000_wp, 0.0000000000_wp, 0.0000000000_wp, &
    0.0000000000_wp, 2.2676713496_wp, 2.2676713496_wp, &
    0.0000000000_wp, 2.2676713496_wp, -2.2676713496_wp, &
    0.0000000000_wp, -2.2676713496_wp, 2.2676713496_wp, &
    0.0000000000_wp, -2.2676713496_wp, -2.2676713496_wp], &
    [3,nat_d4h])

  integer,parameter :: nat_d5d = 21
  integer,parameter :: at_d5d(nat_d5d) = [6,6,6,6,6,1,1,1,1,1,26,6,6,6,6,6,1,1,1,1,1]
  real(wp),parameter :: xyz_d5d(3,nat_d5d) = reshape([ &
    1.8752508225_wp, 1.3624547413_wp, -3.2840605457_wp, &
    1.8752508225_wp, -1.3624547413_wp, -3.2840605457_wp, &
    -0.7162817903_wp, -2.2044978052_wp, -3.2840605457_wp, &
    -2.3179380645_wp, 0.0000000000_wp, -3.2840605457_wp, &
    -0.7162817903_wp, 2.2044978052_wp, -3.2840605457_wp, &
    3.5276840377_wp, 2.5630166456_wp, -3.2695285518_wp, &
    3.5276840377_wp, -2.5630166456_wp, -3.2695285518_wp, &
    -1.3474503159_wp, -4.1470417750_wp, -3.2695285518_wp, &
    -4.3604485463_wp, 0.0000000000_wp, -3.2695285518_wp, &
    -1.3474503159_wp, 4.1470417750_wp, -3.2695285518_wp, &
    -0.0000000000_wp, -0.0000000000_wp, -0.0000000000_wp, &
    -1.8752508225_wp, -1.3624547413_wp, 3.2840605457_wp, &
    -1.8752508225_wp, 1.3624547413_wp, 3.2840605457_wp, &
    0.7162817903_wp, 2.2044978052_wp, 3.2840605457_wp, &
    2.3179380645_wp, -0.0000000000_wp, 3.2840605457_wp, &
    0.7162817903_wp, -2.2044978052_wp, 3.2840605457_wp, &
    -3.5276840377_wp, -2.5630166456_wp, 3.2695285518_wp, &
    -3.5276840377_wp, 2.5630166456_wp, 3.2695285518_wp, &
    1.3474503159_wp, 4.1470417750_wp, 3.2695285518_wp, &
    4.3604485463_wp, -0.0000000000_wp, 3.2695285518_wp, &
    1.3474503159_wp, -4.1470417750_wp, 3.2695285518_wp], &
    [3,nat_d5d])

  integer,parameter :: nat_d5h = 10
  integer,parameter :: at_d5h(nat_d5h) = [6,6,6,6,6,1,1,1,1,1]
  real(wp),parameter :: xyz_d5h(3,nat_d5h) = reshape([ &
    -1.8663524711_wp, 1.2366992873_wp, 0.0051139268_wp, &
    0.6026542002_wp, 2.1530019351_wp, 0.0043775201_wp, &
    2.2371465685_wp, 0.0880473317_wp, 0.0060344776_wp, &
    0.7782640422_wp, -2.1045045776_wp, 0.0074568022_wp, &
    -1.7577693293_wp, -1.3946159551_wp, 0.0071930800_wp, &
    -3.5767552803_wp, 2.3747834075_wp, 0.0043379967_wp, &
    1.1564942170_wp, 4.1315770941_wp, 0.0029612919_wp, &
    4.2898362320_wp, 0.1727660771_wp, 0.0060819138_wp, &
    1.4930016348_wp, -4.0307549451_wp, 0.0085934742_wp, &
    -3.3686738275_wp, -2.6698661099_wp, 0.0083161327_wp], &
    [3,nat_d5h])

  integer,parameter :: nat_d6h = 12
  integer,parameter :: at_d6h(nat_d6h) = [6,6,6,6,6,6,1,1,1,1,1,1]
  real(wp),parameter :: xyz_d6h(3,nat_d6h) = reshape([ &
    2.6136424084_wp, -0.4184042613_wp, 0.0100533430_wp, &
    0.9574486383_wp, -2.4690216681_wp, -0.0142107405_wp, &
    -1.6464616806_wp, -2.0607841334_wp, -0.0269852891_wp, &
    -2.5945939691_wp, 0.3983353698_wp, -0.0087116374_wp, &
    -0.9386647606_wp, 2.4494252082_wp, 0.0204657339_wp, &
    1.6651132775_wp, 2.0409987009_wp, 0.0256813780_wp, &
    4.6419610471_wp, -0.7366152434_wp, 0.0176122475_wp, &
    1.6964071421_wp, -4.3853740339_wp, -0.0241884944_wp, &
    -2.9363320415_wp, -3.6582263184_wp, -0.0521942356_wp, &
    -4.6232338612_wp, 0.7170565780_wp, -0.0156847268_wp, &
    -1.6777366480_wp, 4.3659287520_wp, 0.0383047485_wp, &
    2.9553426863_wp, 3.6388377283_wp, 0.0426322214_wp], &
    [3,nat_d6h])

  integer,parameter :: nat_d7h = 14
  integer,parameter :: at_d7h(nat_d7h) = [6,6,6,6,6,6,6,1,1,1,1,1,1,1]
  real(wp),parameter :: xyz_d7h(3,nat_d7h) = reshape([ &
    -3.0250909355_wp, -0.1819499680_wp, -0.0116177785_wp, &
    -2.0295202814_wp, 2.2508838774_wp, -0.0016445251_wp, &
    0.4931844451_wp, 2.9892921059_wp, -0.0030129336_wp, &
    2.6435520954_wp, 1.4774985135_wp, -0.0143554973_wp, &
    2.8021645025_wp, -1.1465948089_wp, -0.0184410880_wp, &
    0.8494762113_wp, -2.9064326781_wp, -0.0175336090_wp, &
    -1.7439494254_wp, -2.4770856286_wp, -0.0170237108_wp, &
    -5.0696637259_wp, -0.3055858821_wp, -0.0130586120_wp, &
    -3.4008719473_wp, 3.7725554171_wp, 0.0063520778_wp, &
    0.8275050278_wp, 5.0100568745_wp, 0.0037676503_wp, &
    4.4322977666_wp, 2.4760577858_wp, -0.0182010363_wp, &
    4.6978588149_wp, -1.9224255016_wp, -0.0222170612_wp, &
    1.4250021608_wp, -4.8724090140_wp, -0.0195888430_wp, &
    -2.9223717563_wp, -4.1524896636_wp, -0.0214410923_wp], &
    [3,nat_d7h])

  integer,parameter :: nat_d8h = 16
  integer,parameter :: at_d8h(nat_d8h) = [6,6,6,6,6,6,6,6,1,1,1,1,1,1,1,1]
  real(wp),parameter :: xyz_d8h(3,nat_d8h) = reshape([ &
    -3.1953385351_wp, -1.2992913342_wp, -0.0013443022_wp, &
    -3.1782563860_wp, 1.3405447404_wp, -0.0007915785_wp, &
    3.1781269385_wp, -1.3404545521_wp, -0.0014447561_wp, &
    3.1950752347_wp, 1.2993923915_wp, -0.0013802919_wp, &
    -1.3405368625_wp, -3.1782313687_wp, 0.0013485110_wp, &
    -1.2993730227_wp, 3.1954063105_wp, 0.0015601677_wp, &
    1.2993044510_wp, -3.1951175719_wp, 0.0009977162_wp, &
    1.3404668185_wp, 3.1782610779_wp, 0.0007763660_wp, &
    -5.1398876580_wp, -2.0896976163_wp, 0.0007417474_wp, &
    -5.1121989694_wp, 2.1560890106_wp, -0.0002704245_wp, &
    5.1124825940_wp, -2.1559398160_wp, 0.0009894776_wp, &
    5.1396923174_wp, 2.0898768932_wp, -0.0006786475_wp, &
    -2.1559780569_wp, -5.1130969396_wp, -0.0013939168_wp, &
    -2.0898182036_wp, 5.1398167921_wp, 0.0002326967_wp, &
    2.0899305811_wp, -5.1394090818_wp, 0.0007171611_wp, &
    2.1561490482_wp, 5.1124899088_wp, -0.0000599261_wp], &
    [3,nat_d8h])

  integer,parameter :: nat_dinfh = 3
  integer,parameter :: at_dinfh(nat_dinfh) = [8,6,8]
  real(wp),parameter :: xyz_dinfh(3,nat_dinfh) = reshape([ &
    1.9258198936_wp, -0.1348886508_wp, -0.0823920590_wp, &
    4.1878220648_wp, -0.1348886508_wp, -0.0823920590_wp, &
    6.4498242360_wp, -0.1348886508_wp, -0.0823920590_wp], &
    [3,nat_dinfh])

  integer,parameter :: nat_s4 = 13
  integer,parameter :: at_s4(nat_s4) = [6,6,6,6,6,9,1,9,1,9,1,1,9]
  real(wp),parameter :: xyz_s4(3,nat_s4) = reshape([ &
    2.3072044201_wp, -1.2009776440_wp, -0.8184403846_wp, &
    2.4025600003_wp, 1.1865779309_wp, 0.8250733233_wp, &
    0.0157792131_wp, 0.0560870714_wp, 0.0471297695_wp, &
    -2.3250056402_wp, -0.7273555854_wp, 1.2698014694_wp, &
    -2.3238340100_wp, 0.9671429333_wp, -1.0852886106_wp, &
    2.8351938993_wp, -2.9734273653_wp, 0.0561248659_wp, &
    2.8306207621_wp, -1.0219260937_wp, -2.7883475887_wp, &
    2.9927025718_wp, 2.9301526371_wp, -0.0707324488_wp, &
    2.9884506880_wp, 0.9771584818_wp, 2.7741368482_wp, &
    -2.8216634602_wp, 0.1301454382_wp, 3.0596366711_wp, &
    -2.9424736514_wp, -2.6640791987_wp, 1.0395383412_wp, &
    -2.8250271727_wp, 2.9347824661_wp, -0.8342762895_wp, &
    -2.9360107880_wp, 0.1393106099_wp, -2.8526360714_wp], &
    [3,nat_s4])

  integer,parameter :: nat_td = 5
  integer,parameter :: at_td(nat_td) = [6,1,1,1,1]
  real(wp),parameter :: xyz_td(3,nat_td) = reshape([ &
    1.9771448552_wp, -0.1316194246_wp, -0.1036703752_wp, &
    4.0411037285_wp, -0.1316194246_wp, -0.1036703752_wp, &
    1.2891522650_wp, 0.5743255638_wp, 1.7096730195_wp, &
    1.2891522650_wp, -2.0549826742_wp, -0.3989967740_wp, &
    1.2891522650_wp, 1.0858177339_wp, -1.6217062684_wp], &
    [3,nat_td])

  integer,parameter :: nat_oh = 7
  integer,parameter :: at_oh(nat_oh) = [9,16,9,9,9,9,9]
  real(wp),parameter :: xyz_oh(3,nat_oh) = reshape([ &
    1.5011606389_wp, 0.0480557353_wp, -0.1059758411_wp, &
    4.8021909252_wp, 0.0480557353_wp, -0.1059758411_wp, &
    4.8021909252_wp, -3.1466396619_wp, -0.9371340825_wp, &
    4.8021909252_wp, 0.8791950795_wp, -3.3006523411_wp, &
    8.1032212114_wp, 0.0480557353_wp, -0.1059758411_wp, &
    4.8021909252_wp, 3.2427322353_wp, 0.7251635031_wp, &
    4.8021909252_wp, -0.7831025060_wp, 3.0887006589_wp], &
    [3,nat_oh])

  integer,parameter :: nat_ih = 60
  integer,parameter :: at_ih(nat_ih) = [6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, &
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6]
  real(wp),parameter :: xyz_ih(3,nat_ih) = reshape([ &
    0.0000000000_wp, 2.3362684079_wp, 6.2880636797_wp, &
    -2.2219399773_wp, 0.7218753796_wp, 6.2880636797_wp, &
    -1.3732639748_wp, -1.8901040699_wp, 6.2880636797_wp, &
    1.3732639748_wp, -1.8901040699_wp, 6.2880636797_wp, &
    2.2219399773_wp, 0.7218753796_wp, 6.2880636797_wp, &
    4.3552517994_wp, 1.4152158947_wp, 4.9017605947_wp, &
    5.7285157742_wp, -0.4748881751_wp, 3.4578208628_wp, &
    4.9136658693_wp, -2.9829326877_wp, 3.4578208628_wp, &
    2.6917258919_wp, -3.7048080673_wp, 4.9017605947_wp, &
    1.3184619172_wp, -5.5949121372_wp, 3.4578208628_wp, &
    -1.3184619172_wp, -5.5949121372_wp, 3.4578208628_wp, &
    -2.6917258919_wp, -3.7048080673_wp, 4.9017605947_wp, &
    -4.9136658693_wp, -2.9829326877_wp, 3.4578208628_wp, &
    -5.7285157742_wp, -0.4748881751_wp, 3.4578208628_wp, &
    -4.3552517994_wp, 1.4152158947_wp, 4.9017605947_wp, &
    -4.3552517994_wp, 3.7514843026_wp, 3.4578208628_wp, &
    -5.7285157742_wp, 3.3051309920_wp, 1.1215524550_wp, &
    -6.5771917768_wp, 0.6931515425_wp, 1.1215524550_wp, &
    -6.5771917768_wp, -0.6931515425_wp, -1.1215524550_wp, &
    -5.7285157742_wp, -3.3051309920_wp, -1.1215524550_wp, &
    -4.9136658693_wp, -4.4268724195_wp, 1.1215524550_wp, &
    -2.6917258919_wp, -6.0410764752_wp, 1.1215524550_wp, &
    -1.3732639748_wp, -6.4694773877_wp, -1.1215524550_wp, &
    1.3732639748_wp, -6.4694773877_wp, -1.1215524550_wp, &
    2.6917258919_wp, -6.0410764752_wp, 1.1215524550_wp, &
    4.9136658693_wp, -4.4268724195_wp, 1.1215524550_wp, &
    5.7285157742_wp, -3.3051309920_wp, -1.1215524550_wp, &
    6.5771917768_wp, -0.6931515425_wp, -1.1215524550_wp, &
    6.5771917768_wp, 0.6931515425_wp, 1.1215524550_wp, &
    5.7285157742_wp, 3.3051309920_wp, 1.1215524550_wp, &
    4.9136658693_wp, 4.4268724195_wp, -1.1215524550_wp, &
    4.9136658693_wp, 2.9829326877_wp, -3.4578208628_wp, &
    5.7285157742_wp, 0.4748881751_wp, -3.4578208628_wp, &
    4.3552517994_wp, -1.4152158947_wp, -4.9017605947_wp, &
    4.3552517994_wp, -3.7514843026_wp, -3.4578208628_wp, &
    2.2219399773_wp, -5.3014376700_wp, -3.4578208628_wp, &
    0.0000000000_wp, -4.5795622904_wp, -4.9017605947_wp, &
    -2.2219399773_wp, -5.3014376700_wp, -3.4578208628_wp, &
    -4.3552517994_wp, -3.7514843026_wp, -3.4578208628_wp, &
    -4.3552517994_wp, -1.4152158947_wp, -4.9017605947_wp, &
    -5.7285157742_wp, 0.4748881751_wp, -3.4578208628_wp, &
    -4.9136658693_wp, 2.9829326877_wp, -3.4578208628_wp, &
    -4.9136658693_wp, 4.4268724195_wp, -1.1215524550_wp, &
    -2.6917258919_wp, 6.0410764752_wp, -1.1215524550_wp, &
    -1.3184619172_wp, 5.5949121372_wp, -3.4578208628_wp, &
    -2.6917258919_wp, 3.7048080673_wp, -4.9017605947_wp, &
    -1.3732639748_wp, 1.8901040699_wp, -6.2880636797_wp, &
    -2.2219399773_wp, -0.7218753796_wp, -6.2880636797_wp, &
    0.0000000000_wp, -2.3362684079_wp, -6.2880636797_wp, &
    2.2219399773_wp, -0.7218753796_wp, -6.2880636797_wp, &
    1.3732639748_wp, 1.8901040699_wp, -6.2880636797_wp, &
    2.6917258919_wp, 3.7048080673_wp, -4.9017605947_wp, &
    1.3184619172_wp, 5.5949121372_wp, -3.4578208628_wp, &
    2.6917258919_wp, 6.0410764752_wp, -1.1215524550_wp, &
    1.3732639748_wp, 6.4694773877_wp, 1.1215524550_wp, &
    -1.3732639748_wp, 6.4694773877_wp, 1.1215524550_wp, &
    -2.2219399773_wp, 5.3014376700_wp, 3.4578208628_wp, &
    0.0000000000_wp, 4.5795622904_wp, 4.9017605947_wp, &
    2.2219399773_wp, 5.3014376700_wp, 3.4578208628_wp, &
    4.3552517994_wp, 3.7514843026_wp, 3.4578208628_wp], &
    [3,nat_ih])

!========================================================================================!
!========================================================================================!
contains  !> Unit tests for getsym point-group detection
!========================================================================================!
!========================================================================================!

  subroutine collect_getsym(testsuite)
    type(unittest_type),allocatable,intent(out) :: testsuite(:)
!&<
    testsuite = [ &
      new_unittest("symmetry c1                   ",test_c1),    &
      new_unittest("symmetry ci                   ",test_ci),    &
      new_unittest("symmetry cs                   ",test_cs),    &
      new_unittest("symmetry c2                   ",test_c2),    &
      new_unittest("symmetry c2h                  ",test_c2h),   &
      new_unittest("symmetry c2v                  ",test_c2v),   &
      new_unittest("symmetry c3                   ",test_c3),    &
      new_unittest("symmetry c3v                  ",test_c3v),   &
      new_unittest("symmetry c4v                  ",test_c4v),   &
      new_unittest("symmetry c5                   ",test_c5),    &
      new_unittest("symmetry c5v                  ",test_c5v),   &
      new_unittest("symmetry cinfv                ",test_cinfv), &
      new_unittest("symmetry d2                   ",test_d2),    &
      new_unittest("symmetry d2d                  ",test_d2d),   &
      new_unittest("symmetry d2h                  ",test_d2h),   &
      new_unittest("symmetry d3                   ",test_d3),    &
      new_unittest("symmetry d3d                  ",test_d3d),   &
      new_unittest("symmetry d3h                  ",test_d3h),   &
      new_unittest("symmetry d4                   ",test_d4),    &
      new_unittest("symmetry d4h                  ",test_d4h),   &
      new_unittest("symmetry d5d                  ",test_d5d),   &
      new_unittest("symmetry d5h                  ",test_d5h),   &
      new_unittest("symmetry d6h                  ",test_d6h),   &
      new_unittest("symmetry d7h                  ",test_d7h),   &
      new_unittest("symmetry d8h                  ",test_d8h),   &
      new_unittest("symmetry dinfh                ",test_dinfh), &
      new_unittest("symmetry s4                   ",test_s4),    &
      new_unittest("symmetry td                   ",test_td),    &
      new_unittest("symmetry oh                   ",test_oh),    &
      new_unittest("symmetry ih                   ",test_ih)     &
    ]
!&>
  end subroutine collect_getsym

!========================================================================================!

  subroutine test_c1(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c1,at_c1,xyz_c1,sfsym)
    if (sfsym /= 'c1 ') call test_failed(error,'expected c1 , got: '//sfsym)
  end subroutine test_c1

  subroutine test_ci(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_ci,at_ci,xyz_ci,sfsym)
    if (sfsym /= 'ci ') call test_failed(error,'expected ci , got: '//sfsym)
  end subroutine test_ci

  subroutine test_cs(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_cs,at_cs,xyz_cs,sfsym)
    if (sfsym /= 'cs ') call test_failed(error,'expected cs , got: '//sfsym)
  end subroutine test_cs

  subroutine test_c2(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c2,at_c2,xyz_c2,sfsym)
    if (sfsym /= 'c2 ') call test_failed(error,'expected c2 , got: '//sfsym)
  end subroutine test_c2

  subroutine test_c2h(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c2h,at_c2h,xyz_c2h,sfsym)
    if (sfsym /= 'c2h') call test_failed(error,'expected c2h, got: '//sfsym)
  end subroutine test_c2h

  subroutine test_c2v(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c2v,at_c2v,xyz_c2v,sfsym)
    if (sfsym /= 'c2v') call test_failed(error,'expected c2v, got: '//sfsym)
  end subroutine test_c2v

  subroutine test_c3(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c3,at_c3,xyz_c3,sfsym)
    if (sfsym /= 'c3 ') call test_failed(error,'expected c3 , got: '//sfsym)
  end subroutine test_c3

  subroutine test_c3v(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c3v,at_c3v,xyz_c3v,sfsym)
    if (sfsym /= 'c3v') call test_failed(error,'expected c3v, got: '//sfsym)
  end subroutine test_c3v

  subroutine test_c4v(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c4v,at_c4v,xyz_c4v,sfsym)
    if (sfsym /= 'c4v') call test_failed(error,'expected c4v, got: '//sfsym)
  end subroutine test_c4v

  subroutine test_c5(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c5,at_c5,xyz_c5,sfsym)
    if (sfsym /= 'c5 ') call test_failed(error,'expected c5 , got: '//sfsym)
  end subroutine test_c5

  subroutine test_c5v(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_c5v,at_c5v,xyz_c5v,sfsym)
    if (sfsym /= 'c5v') call test_failed(error,'expected c5v, got: '//sfsym)
  end subroutine test_c5v

  subroutine test_cinfv(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_cinfv,at_cinfv,xyz_cinfv,sfsym)
    if (sfsym /= 'cin') call test_failed(error,'expected cin, got: '//sfsym)
  end subroutine test_cinfv

  subroutine test_d2(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d2,at_d2,xyz_d2,sfsym)
    if (sfsym /= 'd2 ') call test_failed(error,'expected d2 , got: '//sfsym)
  end subroutine test_d2

  subroutine test_d2d(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d2d,at_d2d,xyz_d2d,sfsym)
    if (sfsym /= 'd2d') call test_failed(error,'expected d2d, got: '//sfsym)
  end subroutine test_d2d

  subroutine test_d2h(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d2h,at_d2h,xyz_d2h,sfsym)
    if (sfsym /= 'd2h') call test_failed(error,'expected d2h, got: '//sfsym)
  end subroutine test_d2h

  subroutine test_d3(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d3,at_d3,xyz_d3,sfsym)
    if (sfsym /= 'd3 ') call test_failed(error,'expected d3 , got: '//sfsym)
  end subroutine test_d3

  subroutine test_d3d(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d3d,at_d3d,xyz_d3d,sfsym)
    if (sfsym /= 'd3d') call test_failed(error,'expected d3d, got: '//sfsym)
  end subroutine test_d3d

  subroutine test_d3h(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d3h,at_d3h,xyz_d3h,sfsym)
    if (sfsym /= 'd3h') call test_failed(error,'expected d3h, got: '//sfsym)
  end subroutine test_d3h

  subroutine test_d4(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d4,at_d4,xyz_d4,sfsym)
    if (sfsym /= 'd4 ') call test_failed(error,'expected d4 , got: '//sfsym)
  end subroutine test_d4

  subroutine test_d4h(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d4h,at_d4h,xyz_d4h,sfsym)
    if (sfsym /= 'd4h') call test_failed(error,'expected d4h, got: '//sfsym)
  end subroutine test_d4h

  subroutine test_d5d(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d5d,at_d5d,xyz_d5d,sfsym)
    if (sfsym /= 'd5d') call test_failed(error,'expected d5d, got: '//sfsym)
  end subroutine test_d5d

  subroutine test_d5h(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d5h,at_d5h,xyz_d5h,sfsym)
    if (sfsym /= 'd5h') call test_failed(error,'expected d5h, got: '//sfsym)
  end subroutine test_d5h

  subroutine test_d6h(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d6h,at_d6h,xyz_d6h,sfsym)
    if (sfsym /= 'd6h') call test_failed(error,'expected d6h, got: '//sfsym)
  end subroutine test_d6h

  subroutine test_d7h(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d7h,at_d7h,xyz_d7h,sfsym)
    if (sfsym /= 'd7h') call test_failed(error,'expected d7h, got: '//sfsym)
  end subroutine test_d7h

  subroutine test_d8h(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_d8h,at_d8h,xyz_d8h,sfsym)
    if (sfsym /= 'd8h') call test_failed(error,'expected d8h, got: '//sfsym)
  end subroutine test_d8h

  subroutine test_dinfh(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_dinfh,at_dinfh,xyz_dinfh,sfsym)
    if (sfsym /= 'din') call test_failed(error,'expected din, got: '//sfsym)
  end subroutine test_dinfh

  subroutine test_s4(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_s4,at_s4,xyz_s4,sfsym)
    if (sfsym /= 's4 ') call test_failed(error,'expected s4 , got: '//sfsym)
  end subroutine test_s4

  subroutine test_td(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_td,at_td,xyz_td,sfsym)
    if (sfsym /= 'td ') call test_failed(error,'expected td , got: '//sfsym)
  end subroutine test_td

  subroutine test_oh(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_oh,at_oh,xyz_oh,sfsym)
    if (sfsym /= 'oh ') call test_failed(error,'expected oh , got: '//sfsym)
  end subroutine test_oh

  subroutine test_ih(error)
    type(error_type),allocatable,intent(out) :: error
    character(len=3) :: sfsym
    call getsym(.false.,6,nat_ih,at_ih,xyz_ih,sfsym)
    if (sfsym /= 'ih ') call test_failed(error,'expected ih , got: '//sfsym)
  end subroutine test_ih

end module test_getsym
