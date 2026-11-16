# Search theory for search and rescue: the values of the reference
# implementation (rmoriebricklayer 1.1.1, computed once on these inputs and
# printed with 17 significant digits), closed forms, refusals and printing.

ref <- local({
  lr1 <- c(
    0.014831789952537533, 0.24001841429585147, 1, 0.9653392783061816,
    0.41605221397977077, 0.12583481900704341, 0.021287893166983318,
    0.0013439553511591651
  )
  lr2 <- c(
    0, 0.40000000000000002, 0.40000000000000002, 0.40000000000000002,
    0.40000000000000002, 0.40000000000000002, 0, 0
  )
  lr3 <- c(0, 0, 1, 1, 0, 0, 0, 0)
  lr4 <- c(0.65487412763980524, 0.31612805583085951, 0.071320237997985286)
  lr5 <- c(
    0.33653789842741721, 0.39999977067874248, 0.39997106066950366,
    0.00053995921265204139
  )
  sw1 <- 1.8729434674176879
  pod_random <- c(
    0, 0.095162581964040482, 0.39346934028736658, 0.63212055882855767,
    0.81731647594726531, 0.95021293163213605, 0.99966453737209748
  )
  pod_definite <- c(0, 0.10000000000000001, 0.5, 1, 1, 1, 1)
  pod_invcube <- c(
    0, 0.099738816310728007, 0.46911594893005937, 0.78990859455606288,
    0.9668806965043768, 0.99983004752482341, 1
  )
  pod_mbeta <- c(
    0, 0.10000000000000001, 0.5, 0.75, 0.90526771459310007, 0.984375,
    0.9999847412109375
  )
  psw <- c(
    0.96327939300563659, 0.91264875798198719, 0.6649980868389902,
    0.65255500109441444, 0.32388434665789001, 0.32428304071312092, 1,
    0.97017767030193813, 0.76666666666666672, 0.73971688956402015,
    0.32666666666666666, 0.3283690131922038, 0.96277333333333348,
    0.93597876680131586, 0.77333333333333332, 0.71054229092428134,
    0.32533333333333331, 0.33191276415089449
  )
  psw0 <- 0.66499809026450218
  scs <- c(
    0.77450000000000008, 0.33600000000000008, 0.23999999999999999,
    0.18000000000000002, 0.018499999999999999, 0.28381374722838154,
    0.26607538802660768, 0.088691796008869214, 0.13968957871396903,
    0.84000000000000008, 0.80000000000000004, 0.90000000000000002, 0.37
  )
  mp1 <- c(
    1.4604005398928015e-23, 4.6880059015098391e-19, 3.3629118881126611e-16,
    7.0645158385371882e-15, 7.0645158385371882e-15, 3.3629118881126606e-16,
    4.6880059015016481e-19, 1.4604005566922526e-23, 2.0267538661466795e-20,
    6.5060466809338187e-16, 4.6670721384718676e-13, 9.8041834394697589e-12,
    9.8041834394697589e-12, 4.6670721384718666e-13, 6.5060466809224512e-16,
    2.0267538894610729e-20, 1.0553553319456102e-17, 3.3877774550220735e-13,
    2.430201095545672e-10, 5.105157286754838e-09, 5.105157286754838e-09,
    2.4302010955456715e-10, 3.3877774550161539e-13, 1.0553553440856977e-17,
    2.0735465605753046e-15, 6.6562550803671971e-11, 4.7748231999596951e-08,
    1.0030537594983229e-06, 1.0030537594983229e-06, 4.7748231999596944e-08,
    6.6562550803555674e-11, 2.0735465844279691e-15, 1.5489034989192155e-13,
    4.9721076824139799e-09, 3.5667105343833091e-06, 7.4926384930559304e-05,
    7.4926384930559304e-05, 3.5667105343833082e-06, 4.972107682405292e-09,
    1.5489035167367444e-13, 4.4416931454039458e-12, 1.4258200479628666e-07,
    0.00010228031470820931, 0.0021486168156259053, 0.0021486168156259053,
    0.0001022803147082093, 1.4258200479603753e-07, 4.4416931964981532e-12,
    4.9465602441485556e-11, 1.5878865409379111e-06, 0.00113906053825026,
    0.023928403363662781, 0.023928403363662781, 0.0011390605382502598,
    1.5878865409351368e-06, 4.9465603010504095e-11, 2.1651029691633091e-10,
    6.9501586856158767e-06, 0.0049856531239859496, 0.10473430952607582,
    0.10473430952607582, 0.0049856531239859496, 6.9501586856037337e-06,
    2.1651029940691763e-10, 3.7605601785304667e-10, 1.2071707609220545e-05,
    0.0086595644036611165, 0.18191267544279871, 0.18191267544279871,
    0.0086595644036611148, 1.2071707609199454e-05, 3.7605602217893848e-10,
    2.6028296609892695e-10, 8.3552973845368657e-06, 0.0059936206339088116,
    0.12590882338637591, 0.12590882338637591, 0.0059936206339088107,
    8.3552973845222679e-06, 2.6028296909304459e-10, 7.1585889452854289e-11,
    2.2979659555894559e-06, 0.001648431591786411, 0.034628832025259522,
    0.034628832025259522, 0.0016484315917864108, 2.2979659555854409e-06,
    7.1585890276329507e-11, 7.7551725221384303e-12, 2.4894741927225383e-07,
    0.00017858088350870057, 0.0037514734907764997, 0.0037514734907764997,
    0.00017858088350870054, 2.4894741927181883e-07, 7.7551726113486428e-12,
    3.2707009558556895e-13, 1.0499219196571132e-08, 7.531549617008437e-06,
    0.00015821631172128499, 0.00015821631172128499, 7.5315496170084353e-06,
    1.0499219196552787e-08, 3.2707009934796014e-13, 5.3070720583728001e-15,
    1.7036137997604349e-10, 1.2220767678901619e-07, 2.5672336861355064e-06,
    2.5672336861355064e-06, 1.2220767678901616e-07, 1.7036137997574581e-10,
    5.3070721194217359e-15, 3.2796875612080288e-17, 1.0528066939211243e-12,
    7.5522433658448874e-10, 1.5865102818510489e-08, 1.5865102818510489e-08,
    7.5522433658448864e-10, 1.0528066939192848e-12, 3.2796875989353164e-17,
    7.6579344086057105e-20, 2.4582599581587429e-15, 1.7634174979815815e-12,
    3.7044356970756157e-11, 3.7044356970756157e-11, 1.7634174979815813e-12,
    2.4582599581544475e-15, 7.6579344966973622e-20
  )
  mp2 <- c(
    0.0014556365276158336, 0.0026087917615673626, 0.003108379871118274,
    0.0024732537445290403, 0.0013210859669765952, 0.00047671040464310698,
    0.00011707726903005351, 1.9736737915163359e-05, 0.0023104263189222063,
    0.0042118613003099764, 0.0051233185920777661, 0.004180712361648747,
    0.0023031311294576225, 0.00086294011600196842, 0.00022177311004961605,
    3.9448941310697594e-05, 0.0033324140280853402, 0.0061983881298269469,
    0.0077270041147292842, 0.0064972715773443501, 0.0037125516462895212,
    0.001453818231466232, 0.00039372723949300374, 7.4400798197444507e-05,
    0.0043790952453616767, 0.0083418951207790658, 0.010706606099786952,
    0.0093282704896334619, 0.0055642267067756808, 0.0022932328490553682,
    0.00065896642833449849, 0.00013303358364943028, 0.0052593904910977901,
    0.010306719786432876, 0.013693451996544603, 0.012440320634449805,
    0.0078002076387637433, 0.0034068999515801716, 0.001044887500691781,
    0.00022620506372114524, 0.00579496652120998, 0.011744489850936009,
    0.016252558347327795, 0.015502186307664792, 0.010288829356114004,
    0.0047916821517197993, 0.0015749937644447794, 0.00036604985741348918,
    0.0058839244017983025, 0.012407578481281694, 0.018006276096391136,
    0.01815897480956977, 0.012838225894553832, 0.0064038019537844956,
    0.0022594813923308989, 0.00056282595637243576, 0.0055337178157704205,
    0.012223493118775058, 0.018734172616897908, 0.020105463669011225,
    0.015214473573887824, 0.008145160220986096, 0.0030812889118480105,
    0.00081920670862759362, 0.0048482301162408297, 0.011297333729410724,
    0.018407638325544925, 0.021130769221654681, 0.017158525118569811,
    0.0098519331224079781, 0.0039808675374132562, 0.001122928241229022,
    0.0039808675374132579, 0.0098519331224079833, 0.017158525118569808,
    0.021130769221654674, 0.018407638325544914, 0.011297333729410712,
    0.0048482301162408236, 0.001441075525751489, 0.0030812889118480135,
    0.0081451602209860995, 0.01521447357388782, 0.020105463669011218,
    0.018734172616897887, 0.012223493118775046, 0.0055337178157704118,
    0.001721063850354608, 0.0022594813923309011, 0.0064038019537844982,
    0.012838225894553832, 0.01815897480956976, 0.018006276096391122,
    0.012407578481281687, 0.0058839244017982947, 0.001902119879962152,
    0.0015749937644447781, 0.0047916821517197923, 0.010288829356113992,
    0.015502186307664773, 0.016252558347327768, 0.011744489850935999,
    0.0057949665212099766, 0.0019356479479840779, 0.0010448875006917801,
    0.0034068999515801711, 0.0078002076387637363, 0.012440320634449787,
    0.013693451996544592, 0.010306719786432864, 0.0052593904910977866,
    0.0018058384404297553, 0.00065896642833449697, 0.0022932328490553656,
    0.0055642267067756712, 0.0093282704896334497, 0.010706606099786941,
    0.0083418951207790623, 0.0043790952453616775, 0.0015388843948119695,
    0.00039372723949300374, 0.0014538182314662331, 0.0037125516462895199,
    0.0064972715773443449, 0.007727004114729286, 0.0061983881298269513,
    0.0033324140280853445, 0.0011942200273645737
  )
  mp3 <- c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0.044642857142857137, 0.044642857142857137, 0.044642857142857137,
    0.008928571428571435, 0, 0, 0, 0, 0.044642857142857137,
    0.044642857142857137, 0.044642857142857137, 0.008928571428571435, 0, 0, 0,
    0, 0.044642857142857137, 0.044642857142857137, 0.044642857142857137,
    0.008928571428571435, 0, 0, 0, 0, 0.044642857142857137,
    0.044642857142857137, 0.044642857142857137, 0.008928571428571435, 0, 0, 0,
    0, 0.044642857142857137, 0.044642857142857137, 0.044642857142857137,
    0.008928571428571435, 0, 0, 0, 0, 0.044642857142857137,
    0.044642857142857137, 0.044642857142857137, 0.008928571428571435, 0, 0, 0,
    0, 0.044642857142857137, 0.044642857142857137, 0.044642857142857137,
    0.008928571428571435, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0
  )
  mp4 <- c(
    0.048346067491425775, 0.076584984509605214, 0.048346067491425782,
    0.096692134982851549, 0.15316996901921043, 0.096692134982851563,
    0.096692134982851549, 0.15316996901921043, 0.096692134982851563, 0, 0, 0
  )
  mp5 <- c(
    0.16937450979606611, 0.084687254898033054, 0.16937450979606611,
    0.16937450979606614, 0.084687254898033068, 0.16937450979606614
  )
  crc <- c(
    0, 0, 0, 0.83125705948441186, 0.40192908886352718, 0.92131773192356137,
    1.9454302308307641, 0.99557816332377636, 5.0462650440403198, 0, 0, 0,
    1.6625141189688237, 0.04756441971481129, 0.23032943298089034,
    3.8908604616615281, 0.5640729130594262, 1.26156626101008
  )
  opt2 <- c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.051992840185673052,
    0.051992840185673052, 0, 0, 0, 0, 0, 0, 0.18960158878893268,
    0.18960158878893268, 0, 0, 0, 0, 0, 0, 0.25840557102539269,
    0.25840557102539269, 0, 0, 0, 0, 0, 0, 0.25840557102539269,
    0.25840557102539269, 0, 0, 0, 0, 0, 0, 0.18960158878893268,
    0.18960158878893268, 0, 0, 0, 0, 0, 0, 0.051992840185672927,
    0.051992840185672927, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0.08148325486612211, 0.034194431042445066
  )
  opt15 <- c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0.076933674272789077, 0.076933674272789077, 0, 0, 0, 0, 0,
    0.084475806540037016, 0.35215741717544913, 0.35215741717544913,
    0.084475806540036891, 0, 0, 0, 0, 0.29089088520325285, 0.55857249583866497,
    0.55857249583866497, 0.29089088520325285, 0, 0, 0, 0, 0.42849963380651246,
    0.69618124444192464, 0.69618124444192464, 0.42849963380651235, 0, 0, 0, 0,
    0.49730361604297252, 0.76498522667838464, 0.76498522667838464,
    0.49730361604297235, 0, 0, 0, 0, 0.49730361604297252, 0.76498522667838464,
    0.76498522667838464, 0.49730361604297235, 0, 0, 0, 0, 0.42849963380651246,
    0.69618124444192464, 0.69618124444192464, 0.42849963380651235, 0, 0, 0, 0,
    0.29089088520325268, 0.55857249583866486, 0.55857249583866486,
    0.29089088520325257, 0, 0, 0, 0, 0.084475806540037016, 0.35215741717544913,
    0.35215741717544913, 0.084475806540036891, 0, 0, 0, 0, 0,
    0.076933674272789576, 0.076933674272789576, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.36933896028906665,
    0.015203647204722245
  )
  opt60 <- c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.12504312377300345, 0.12504312377300345, 0,
    0, 0, 0, 0, 0.2702125974320011, 0.53789420806741317, 0.53789420806741317,
    0.27021259743200088, 0, 0, 0, 0.07869519974361007, 0.61424809577026795,
    0.88192970640567991, 0.88192970640567991, 0.61424809577026773,
    0.078695199743610306, 0, 0, 0.35391894264627, 0.88947183867292789,
    1.15715344930834, 1.15715344930834, 0.88947183867292778,
    0.35391894264627044, 0, 0, 0.56033402130948595, 1.0958869173361436,
    1.3635685279715559, 1.3635685279715559, 1.0958869173361436,
    0.56033402130948617, 0, 0, 0.69794276991274551, 1.2334956659394034,
    1.5011772765748155, 1.5011772765748155, 1.2334956659394034,
    0.69794276991274584, 0, 0, 0.76674675214920562, 1.3022996481758633,
    1.5699812588112754, 1.5699812588112754, 1.3022996481758633,
    0.76674675214920596, 0, 0, 0.76674675214920562, 1.3022996481758633,
    1.5699812588112754, 1.5699812588112754, 1.3022996481758633,
    0.76674675214920596, 0, 0, 0.69794276991274551, 1.2334956659394034,
    1.5011772765748155, 1.5011772765748155, 1.2334956659394034,
    0.69794276991274584, 0, 0, 0.56033402130948573, 1.0958869173361436,
    1.3635685279715557, 1.3635685279715557, 1.0958869173361436,
    0.56033402130948606, 0, 0, 0.35391894264627, 0.88947183867292789,
    1.15715344930834, 1.15715344930834, 0.88947183867292778,
    0.35391894264627044, 0, 0, 0.078695199743610555, 0.61424809577026829,
    0.88192970640568036, 0.88192970640568036, 0.61424809577026829,
    0.078695199743610805, 0, 0, 0, 0.27021259743200088, 0.53789420806741306,
    0.53789420806741306, 0.27021259743200082, 0, 0, 0, 0, 0,
    0.12504312377300447, 0.12504312377300447, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0.72830983475943967, 0.0041935251103067013
  )
  opt3 <- c(
    1.90772034023008, 0.59656485559195949, 0.49571480417796249,
    1.90772034023008, 0.29828242779597974, 0.99142960835592497,
    0.85158165543626929, 0.25790827718134646, 0.62895413859067317,
    0.54026789701347122, 0.074209172281865368
  )
  sq1 <- c(
    1, 2, 4, 0.64321684927634259, 0.93595342672330695, 1.1183365957438716,
    2.5728673971053704, 3.7438137068932278, 4.4733463829754863,
    0.30376389417226796, 0.43816504002565243, 0.53392748708707882,
    0.60426186104311186, 0.28538525543333454, 0.39978355993755993,
    0.55114777792267788, 0.27941523356465869, 0.38366715332840928,
    0.16741879528618497, 0.12242998699863569, 0.20485043905449052,
    0.16741879528618497, 0.26935170135460584, 0.41902532612653132
  )
  sq2 <- c(
    0.93651343816339616, 1.2247448673296295, 0.53264495442026583,
    0.6103053266756463, 0.29098172233948655, 0.30515266474083896,
    0.2909817223394866, 0.50734053911755228
  )
  sq3 <- c(1.0275995216848992, 0.36639537509130282)
  environment()
})

test_that("every function equals the reference implementation", {
  xs <- c(-3, -0.7, 0, 0.2, 0.5, 1, 2.5, 10)
  tol <- 1e-10
  expect_equal(as.numeric(ltrlrng(xs, 1.3)), ref$lr1, tolerance = tol)
  expect_equal(
    as.numeric(ltrlrng(xs, 1, "mbeta", 0.4)), ref$lr2,
    tolerance = tol
  )
  expect_equal(
    as.numeric(ltrlrng(xs, 0.8, "definite")), ref$lr3,
    tolerance = tol
  )
  expect_equal(
    as.numeric(ltrlrng(c(0, 0.5, 1), 0.8, navsd = 0.3)), ref$lr4,
    tolerance = tol
  )
  expect_equal(
    as.numeric(ltrlrng(c(-1, 0, 0.3, 2), 1, "mbeta", 0.4, navsd = 0.25)),
    ref$lr5,
    tolerance = tol
  )
  xg <- seq(-5, 5, by = 0.05)
  p <- ltrlrng(xg, 2)
  expect_equal(swpwdth(xg, p), ref$sw1, tolerance = tol)
  expect_identical(swpwdth(p), swpwdth(xg, p))
  cv <- c(0, 0.1, 0.5, 1, 1.7, 3, 8)
  for (md in c("random", "definite", "invcube")) {
    expect_equal(
      as.numeric(srchpod(cv, md)), ref[[paste0("pod_", md)]],
      tolerance = tol
    )
  }
  expect_equal(
    as.numeric(srchpod(cv, "mbeta", m = 0.5)), ref$pod_mbeta,
    tolerance = tol
  )
  r <- NULL
  for (md in c("invcube", "definite", "mbeta")) {
    for (S in c(0.6, 1.3, 3)) {
      for (sd in c(0, 0.25)) {
        r <- c(
          r,
          pswppod(1, S, md, m = 0.8, navsd = sd, npts = 101, ntrk = 20)
        )
      }
    }
  }
  expect_equal(r, ref$psw, tolerance = tol)
  expect_equal(pswppod(1, 1.3), ref$psw0, tolerance = tol)
  poc <- c(0.4, 0.3, 0.2, 0.05)
  s <- srchscs(poc, rbind(c(0.8, 0.5, 0, 0.3), c(0.2, 0.6, 0.9, 0.1)))
  expect_equal(
    c(s$pos, s$posarea, s$pocupd, s$cumpod), ref$scs,
    tolerance = tol
  )
  xe <- seq(-4, 4, by = 0.5)
  ye <- seq(-3, 5, by = 1)
  expect_equal(
    as.vector(datumpmap(xe, ye, "point", c(0.3, 1), 0.5)$poc), ref$mp1,
    tolerance = tol
  )
  expect_equal(
    as.vector(
      datumpmap(xe, ye, "line", c(-2, -1), 1.5, c(3, 2), nodes = 32)$poc
    ),
    ref$mp2,
    tolerance = tol
  )
  expect_equal(
    as.vector(datumpmap(xe, ye, "area", c(-1, 0), 1, c(2.5, 3.2))$poc),
    ref$mp3,
    tolerance = tol
  )
  expect_equal(
    as.vector(
      datumpmap(
        0:4, c(-1.5, -0.5, 0.5, 1.5), "corridor", c(0.5, 0),
        end = c(3, 0)
      )$poc
    ),
    ref$mp4,
    tolerance = tol
  )
  expect_equal(
    as.vector(
      datumpmap(
        c(-1, 0, 1), c(-2, 0, 0.5, 2), "corridor", c(0, -1),
        sigma = 0.7, end = c(0, 1.5)
      )$poc
    ),
    ref$mp5,
    tolerance = tol
  )
  o <- NULL
  for (sg in c(0.5, 2)) {
    for (E in c(0, 1, 30)) {
      o <- c(o, unlist(optcrcsrch(sg, 1.5, E)[c("radius", "pos", "cntrdns")]))
    }
  }
  expect_equal(unname(o), ref$crc, tolerance = tol)
  m <- datumpmap(xe, ye, sigma = 1.5)
  for (E in c(2, 15, 60)) {
    a <- opteffort(m, W = 0.8, effort = E)
    expect_equal(
      c(a$effort, a$pos, a$lambda), ref[[paste0("opt", E)]],
      tolerance = tol
    )
  }
  a <- opteffort(c(0.5, 0.2, 0.1), c(1, 2, 0.5), W = 1, effort = 3)
  expect_equal(
    c(a$effort, a$coverage, a$pod, a$pos, a$lambda), ref$opt3,
    tolerance = tol
  )
  # optsqr() finds its factor by minimising, and the place of a minimum is
  # known only to about the square root of the machine precision, so the
  # last digits differ between platforms
  q <- optsqr(c(4, 4, 8), E = 2)
  expect_equal(
    unname(unlist(q[c(
      "cumeff", "factor", "side", "poc", "coverage", "pod", "pos", "cumpos"
    )])),
    ref$sq1,
    tolerance = 1e-6
  )
  q <- optsqr(c(2, 3), "mbeta", m = 0.5)
  expect_equal(
    unname(unlist(q[c("factor", "poc", "pos", "cumpos")])), ref$sq2,
    tolerance = 1e-6
  )
  q <- optsqr(4, "random")
  expect_equal(unname(unlist(q[c("factor", "pos")])), ref$sq3, tolerance = 1e-6)
})

test_that("closed forms and hand-checked values", {
  # inverse cube: 1 - exp(-1 / pi) half a sweep width off the track
  expect_equal(as.numeric(ltrlrng(0.5, 1)), 1 - exp(-1 / pi), tolerance = 1e-15)
  expect_equal(
    round(as.numeric(ltrlrng(c(0.5, 1, 2), 1)), 6),
    c(0.272623, 0.076494, 0.019698)
  )
  # its sweep width on [-L, L]: 2 (L p(L) + sqrt(pi a) erfc(sqrt(a) / L)),
  # a = W^2 / (4 pi)
  W <- 1.7
  L <- 40
  a <- W^2 / (4 * pi)
  xs <- seq(-L, L, length.out = 40001)
  exact <- 2 *
    (L * (1 - exp(-a / L^2)) +
       sqrt(pi * a) * 2 * stats::pnorm(-sqrt(2 * a) / L))
  expect_equal(swpwdth(ltrlrng(xs, W)), exact, tolerance = 1e-5)
  expect_equal(exact, W, tolerance = 0.02)
  # definite range and M-Beta curves sweep W
  xd <- seq(-3, 3, by = 0.001)
  expect_equal(swpwdth(ltrlrng(xd, 2, "definite")), 2, tolerance = 2e-3)
  expect_equal(swpwdth(ltrlrng(xd, 1, "mbeta", m = 0.5)), 1, tolerance = 2e-3)
  cv <- c(0, 0.5, 1, 2)
  expect_equal(as.numeric(srchpod(cv)), 1 - exp(-cv), tolerance = 1e-15)
  expect_equal(as.numeric(srchpod(cv, "definite")), c(0, 0.5, 1, 1))
  expect_equal(round(as.numeric(srchpod(1, "invcube")), 4), 0.7899)
  # parallel sweeps: inverse cube = erf form, definite range at S >= W is W / S
  expect_equal(
    pswppod(1, 1.3),
    as.numeric(srchpod(1 / 1.3, "invcube")),
    tolerance = 1e-12
  )
  expect_equal(pswppod(1, 2, "definite", npts = 2001), 0.5, tolerance = 1e-3)
  # an even number of points is raised to the next odd one
  expect_identical(pswppod(1, 1.3, npts = 100), pswppod(1, 1.3, npts = 101))
  expect_lt(
    pswppod(1, 1.3, navsd = 0.5, npts = 201),
    pswppod(1, 1.3, npts = 201) + 1e-12
  )
  # Bayes on three areas
  s <- srchscs(c(0.5, 0.3, 0.2), c(0.8, 0.5, 0))
  expect_equal(s$pos, 0.55)
  expect_equal(s$pocupd, c(0.1, 0.15, 0.2) / 0.45, tolerance = 1e-8)
  s2 <- srchscs(c(0.5, 0.5), rbind(c(0.5, 0), c(0.5, 0)))
  expect_equal(s2$cumpod, c(0.75, 0))
  expect_equal(s2$pocupd, c(0.2, 0.8))
  # a point datum on one cell [-1, 1]^2: (2 Phi(1) - 1)^2
  m1 <- datumpmap(c(-1, 1), c(-1, 1))
  expect_equal(m1$poc[1, 1], (2 * stats::pnorm(1) - 1)^2, tolerance = 1e-15)
  m <- datumpmap(c(-1, 0, 1), c(-1, 0, 1))
  expect_equal(round(c(m$poc[1, 1], m$outside), 6), c(0.116516, 0.533935))
  ma <- datumpmap(0:4, 0:2, "area", c(1, 0.5), end = c(3, 1.5))
  expect_equal(sum(ma$poc), 1, tolerance = 1e-8)
  expect_equal(
    ma$poc,
    rbind(c(0, 0.25, 0.25, 0), c(0, 0.25, 0.25, 0)),
    tolerance = 1e-8
  )
  # the segment's probability is that of its midpoint for a tiny line
  ml <- datumpmap(-3:3, -3:3, "line", c(0, 0), end = c(1e-9, 0))
  expect_equal(ml$poc, datumpmap(-3:3, -3:3)$poc, tolerance = 1e-8)
  mg <- datumpmap(0:3, 0:2, "general", wts = matrix(1:6, 2))
  expect_equal(mg$poc, matrix(1:6, 2) / 21, tolerance = 1e-15)
  # Koopman: u = 1 at effort pi sigma^2 / W gives POS 1 - 2 / e
  o <- optcrcsrch(1, 1, pi)
  expect_equal(o$pos, 1 - 2 / exp(1), tolerance = 1e-15)
  expect_equal(o$radius, sqrt(2), tolerance = 1e-15)
  expect_equal(round(optcrcsrch(1, 1, pi / 4)$pos, 6), 0.090204)
  # two equal cells share the effort; a fine grid approaches Koopman's
  # circular optimum
  expect_equal(
    opteffort(c(0.5, 0.5), c(1, 1), 1, 2)$effort,
    c(1, 1),
    tolerance = 1e-10
  )
  g <- seq(-6, 6, by = 0.4)
  fine <- opteffort(datumpmap(g, g), W = 1, effort = pi)
  expect_equal(fine$pos, 1 - 2 / exp(1), tolerance = 2e-2)
  expect_lt(fine$pos, 1 - 2 / exp(1))
  expect_equal(sum(fine$effort), pi, tolerance = 1e-10)
  expect_identical(dim(fine$effort), c(length(g) - 1L, length(g) - 1L))
  # an effort beyond the first bisection floor is still spent in full
  big <- opteffort(c(0.6, 0.4), c(1, 1), 1, 500)
  expect_equal(sum(big$effort), 500, tolerance = 1e-8)
  expect_equal(opteffort(c(0.6, 0.4), c(1, 1), 1, 0)$pos, 0, tolerance = 1e-8)
  # cell areas given as a matrix shape the allocation
  expect_identical(
    dim(opteffort(c(0.5, 0.5), matrix(1, 1, 2), 1, 1)$effort), c(1L, 2L)
  )
})

test_that("bad inputs are refused with the function named", {
  expect_error(ltrlrng(1, 1, "bogus"), "ltrlrng: `model` must be one of")
  expect_error(ltrlrng(1, -1), "ltrlrng: `W` must be positive")
  expect_error(ltrlrng(1, c(1, 2)), "`W` must be a single number")
  expect_error(ltrlrng("a", 1), "`x` must be numeric")
  expect_error(ltrlrng(factor(1), 1), "`x` must be numeric")
  expect_error(ltrlrng(TRUE, 1), "`x` must be numeric, not logical")
  expect_error(ltrlrng(NULL, 1), "`x` must be numeric, got NULL")
  expect_error(ltrlrng(numeric(), 1), "`x` must have at least one value")
  expect_error(ltrlrng(c(1, NA), 1), "`x` must be finite")
  expect_error(ltrlrng(NA, 1), "`x` must be finite")
  expect_error(ltrlrng(Inf, 1), "`x` must be finite")
  expect_error(ltrlrng(1, 1, "mbeta", m = 1.5), "`m` must be between 0 and 1")
  expect_error(
    swpwdth(c(2, 1), c(0, 1)),
    "swpwdth: `x` must hold at least two strictly increasing"
  )
  expect_error(swpwdth(matrix(1:4, 2), 1:4), "not a matrix")
  expect_error(swpwdth(1:3, c(0, 1)), "`p` must have one value per `x`")
  expect_error(swpwdth(1:2, c(0, 2)), "`p` must be between 0 and 1")
  expect_error(srchpod(-1), "srchpod: `coverage` must be at least 0")
  expect_error(pswppod(1, 0), "pswppod: `S` must be positive")
  expect_error(pswppod(1, 1, navsd = -1), "`navsd` must be at least 0")
  expect_error(pswppod(1, 1, npts = 3.5), "`npts` must be a whole number")
  expect_error(pswppod(1, 1, ntrk = 0), "`ntrk` must be between 1")
  expect_error(
    srchscs(c(0.8, 0.8), c(0.5, 0.5)),
    "srchscs: `poc` must sum to at most 1"
  )
  expect_error(
    srchscs(c(0, 0), c(0.5, 0.5)),
    "at least one cell with a positive"
  )
  expect_error(
    srchscs(c(0.5, 0.5), c(0.5, 0.5, 0.5)),
    "one value \\(or column\\) per area"
  )
  expect_error(
    srchscs(c(0.5, 0.5), c(0.5, -0.1)),
    "`pod` must be between 0 and 1"
  )
  expect_warning(
    {
      s <- srchscs(1, 1)
    },
    "POS = 1"
  )
  expect_true(is.na(s$pocupd))
  expect_error(datumpmap(1, 1:2), "datumpmap: `xedg` must hold at least two")
  expect_error(datumpmap(1:2, 1:2, center = 1), "`center` must be two numbers")
  expect_error(datumpmap(1:2, 1:2, "line"), "a line datum needs `end`")
  expect_error(
    datumpmap(1:2, 1:2, "line", end = 1),
    "`end` must be two numbers"
  )
  expect_error(
    datumpmap(1:2, 1:2, "area", end = c(0, 1)),
    "differ in x and in y"
  )
  expect_error(datumpmap(1:2, 1:2, sigma = 0), "`sigma` must be positive")
  expect_error(
    datumpmap(1:2, 1:2, "line", end = c(1, 1), nodes = 1),
    "`nodes` must be between 2 and 1000"
  )
  expect_error(datumpmap(1:2, 1:2, "volume"), "`datum` must be one of")
  expect_error(
    datumpmap(1:3, 1:3, "corridor", end = c(1, 1)),
    "must run along x or along y"
  )
  expect_error(
    datumpmap(1:3, 1:3, "corridor", end = c(0, 0)),
    "`center` and `end` to differ"
  )
  expect_error(datumpmap(1:3, 1:3, "general"), "needs `wts`")
  expect_error(
    datumpmap(1:3, 1:3, "general", wts = matrix(1, 3, 3)),
    "one row per y cell \\(2\\)"
  )
  expect_error(
    datumpmap(1:3, 1:3, "general", wts = matrix(0, 2, 2)),
    "at least one positive weight"
  )
  expect_error(optcrcsrch(1, 1, -1), "optcrcsrch: `effort` must be at least 0")
  expect_error(
    opteffort(c(0.5, 0.5), W = 1, effort = 1),
    "opteffort: `area` must be given"
  )
  expect_error(
    opteffort(c(0.5, 0.5), 1, W = 1, effort = 1),
    "`area` must have one value per cell"
  )
  expect_error(
    opteffort(c(0.5, 0.5), c(1, 0), W = 1, effort = 1),
    "`area` must be positive"
  )
  expect_error(
    opteffort(c(0.5, 0.5), c(1, 1), W = 1, effort = 1, tol = 0),
    "`tol` must be positive"
  )
  m <- "opteffort: `W` = %s is out of range for these cells"
  expect_error(
    opteffort(c(0.5, 0.5), c(1, 1), W = 5e-324, effort = 1),
    sprintf(m, "4.940656e-324"),
    fixed = TRUE
  )
  expect_error(
    opteffort(c(0.5, 0.5), c(1, 1e-300), W = 1e300, effort = 1),
    sprintf(m, "1e+300"),
    fixed = TRUE
  )
  # the kernels refuse what the wrappers would never pass them
  expect_error(
    .Call(C_clab_r_sarlr, 1, -1, 0L, 1, 0),
    "ltrlrng: a size or an index"
  )
})

test_that("srchpod() takes sweep width, effort and area for the coverage", {
  expect_equal(
    as.vector(srchpod(W = 0.6, effort = 12, area = 4, model = "invcube")),
    as.vector(srchpod(0.6 * 12 / 4, "invcube")),
    tolerance = 1e-12
  )
  expect_identical(
    attr(srchpod(W = 2, effort = c(1, 2), area = 4), "coverage"),
    c(0.5, 1)
  )
  m <- datumpmap(0:3, 0:2, center = c(1, 1))
  eff <- matrix(1:6, 2)
  pd <- srchpod(W = 0.5, effort = eff, area = m)
  expect_equal(attr(pd, "coverage"), 0.5 * as.vector(eff), tolerance = 1e-12)
  expect_error(srchpod(1, W = 0.6), "srchpod: give either `coverage`")
  expect_error(srchpod(W = 0.6, effort = 1), "srchpod: give `coverage`")
  expect_error(srchpod(), "srchpod: give `coverage`")
  expect_error(
    srchpod(W = 1, effort = 1:2, area = 1:3),
    "srchpod: `effort` and `area` must have the same length"
  )
  expect_error(srchpod(W = 1, effort = -1, area = 1), "srchpod: `effort`")
})

test_that("srchscs() takes a map and an allocation and returns the map", {
  m <- datumpmap(seq(-3, 3), seq(-3, 3), sigma = 1.5)
  a <- opteffort(m, W = 1, effort = 10)
  s <- srchscs(m, a)
  plain <- srchscs(as.vector(m$poc), as.vector(a$pod))
  expect_equal(s$pos, plain$pos, tolerance = 1e-12)
  expect_equal(as.vector(s$map$poc), plain$pocupd, tolerance = 1e-12)
  expect_s3_class(s$map, "datumpmap")
  expect_equal(sum(s$map$poc) + s$map$outside, 1, tolerance = 1e-12)
  expect_null(plain$map)
  # a POD matrix laid out like the map is one search
  expect_equal(srchscs(m, a$pod)$pos, s$pos, tolerance = 1e-12)
  b <- opteffort(s$map, W = 1, effort = 10)
  expect_s3_class(b, "opteffort")
  out <- capture.output(print(s))
  expect_true(any(grepl("updated in $map", out, fixed = TRUE)))
  expect_true(any(grepl("The 20 most likely of 36 cells", out, fixed = TRUE)))
  expect_error(srchscs(m, c(0.5, 0.5)), "srchscs: `pod` must have one value")
})

test_that("the help page's own search runs step by step", {
  own <- datumpmap(seq(-3, 7), seq(-2, 8), center = c(2, 3), sigma = 1.5)
  walked <- matrix(0, nrow(own$poc), ncol(own$poc))
  walked[5:6, 5:6] <- 12
  pd <- srchpod(W = 0.6, effort = walked, area = own, model = "invcube")
  expect_identical(sum(pd > 0), 4L)
  day1 <- srchscs(own, pd)
  expect_equal(day1$pos, sum(own$poc * as.vector(pd)), tolerance = 1e-12)
  expect_lt(sum(day1$map$poc[5:6, 5:6]), sum(own$poc[5:6, 5:6]))
  day2 <- opteffort(day1$map, W = 0.6, effort = 40)
  expect_equal(sum(day2$effort), 40, tolerance = 1e-8)
})

test_that("the optimal search squares are the ones the book reports", {
  o4 <- optsqr(4)
  expect_equal(round(o4$factor, 2), 1.02)
  expect_equal(
    round(c(o4$poc, o4$coverage, o4$pod, o4$pos), c(2, 2, 2, 3)),
    c(0.59, 0.96, 0.77, 0.458)
  )
  o2 <- optsqr(2)
  expect_equal(
    round(c(o2$factor, o2$poc, o2$coverage, o2$pod, o2$pos), 2),
    c(0.81, 0.44, 0.76, 0.66, 0.29)
  )
  o8 <- optsqr(8)
  expect_equal(
    round(c(o8$factor, o8$poc, o8$coverage, o8$pod, o8$pos), 2),
    c(1.29, 0.76, 1.21, 0.87, 0.66)
  )
  expect_named(
    o4,
    c(
      "search", "effort", "releff", "cumeff", "factor", "side", "poc",
      "coverage", "pod", "pos", "cumpos"
    )
  )
  expect_identical(optsqr(16, E = 2)$factor, o4$factor)
  expect_identical(optsqr(16, E = 2)$side, 2 * o4$factor * 2)
  s <- optsqr(c(4, 4))
  expect_equal(
    s$cumpos[2],
    1 - (1 - s$pos[1]) * (1 - s$pos[2]),
    tolerance = 1e-12
  )
  expect_gt(s$factor[2], s$factor[1])
  expect_lt(optsqr(4, "random")$pos, o4$pos)
  expect_equal(
    optsqr(4, "mbeta", m = 1)$pos,
    optsqr(4, "definite")$pos,
    tolerance = 1e-6
  )
  expect_error(optsqr(-1), "optsqr: `effort` must be positive")
  expect_error(optsqr(1, E = 0), "`E` must be positive")
})

test_that("srchfig() lists the book's figures", {
  tab <- srchfig()
  expect_identical(nrow(tab), 36L)
  expect_named(tab, c("fig", "caption", "how", "yours"))
  expect_setequal(unique(tab$how), c("model", "schematic", "approximation"))
  expect_identical(srchfig(5, 8)$caption, "Inverse cube")
  expect_identical(srchfig("6-1")$yours, "optsqr(effort = 2)")
  # each figure's own-numbers line that computes without drawing runs
  for (y in tab$yours[!grepl("^plot|^with", tab$yours)]) {
    expect_type(unlist(eval(str2lang(y))), "double")
  }
  expect_error(srchfig("9-9"), "srchfig: `fig` must be one of the book's")
  expect_error(srchfig(c("3-1", "3-2")), "`fig` must be a single")
  expect_error(srchfig(3.5, 1), "`fig` must be a whole number")
})

test_that("the search results print", {
  out <- c(
    capture.output(print(ltrlrng(seq(-3, 3, by = 0.25), 1, "mbeta", 0.5, 0.2))),
    capture.output(print(srchpod(seq(0, 3, by = 0.1), "mbeta", m = 0.5))),
    capture.output(print(srchscs(c(0.5, 0.5), rbind(c(0.5, 0), c(0.5, 0))))),
    capture.output(print(
      datumpmap(0:4, 0:2, "area", c(1, 0.5), end = c(3, 1.5))
    )),
    capture.output(print(optcrcsrch(1, 1, pi))),
    capture.output(print(
      opteffort(datumpmap(seq(-6, 6), seq(-6, 6)), W = 1, effort = 60)
    )),
    capture.output(print(optsqr(c(4, 2))))
  )
  expect_true(any(grepl("M-Beta model", out, fixed = TRUE)))
  expect_true(any(grepl("navigation error sd 0.2", out)))
  expect_true(any(grepl("and 5 more", out, fixed = TRUE)))
  expect_true(any(grepl("after 2 searches", out, fixed = TRUE)))
  expect_true(any(grepl("area datum", out, fixed = TRUE)))
  expect_true(any(grepl("radius 1.414", out)))
  expect_true(any(grepl("more cells", out, fixed = TRUE)))
  expect_true(any(grepl("Optimal square searches: 2", out, fixed = TRUE)))
  capture.output(r <- withVisible(print(srchpod(1))))
  expect_false(r$visible)
  expect_s3_class(r$value, "srchpod")
})

test_that("datumpmap() checks `end` whatever the datum", {
  xe <- seq(-3, 3, 1)
  for (e in list("a", sum, TRUE, list(1, 2))) {
    expect_error(datumpmap(xe, xe, "point", end = e), "^datumpmap: `end`")
  }
})
