# Scaling and ideal points (src/svt.cpp). Each fit is checked against the
# values of the rmorie function it ports, computed once on the same small
# input; then the behaviour carried over from that function's tests. The
# iterative fits (SMACOF, unfolding, EM IRT) stop at a convergence
# tolerance, and the last digits of their iterates follow the order of
# floating-point sums, so they are compared at 1e-9; the rest at 1e-10.

test_that("votagr() equals the reference agreement matrix", {
  set.seed(1)
  v <- matrix(rbinom(30, 1, 0.5), 5)
  v[2, 3] <- NA
  r <- votagr(v)
  # rmorie 1.5.0: R$sptag(), the same input,
  # its values printed with 17 significant digits
  ref <- c(1, 0.59999999999999998, 0.5, 0.5, 0.33333333333333331,
        0.59999999999999998, 1, 0.59999999999999998, 0.40000000000000002,
        0.80000000000000004, 0.5, 0.59999999999999998, 1,
        0.33333333333333331, 0.5, 0.5, 0.40000000000000002,
        0.33333333333333331, 1, 0.16666666666666666, 0.33333333333333331,
        0.80000000000000004, 0.5, 0.16666666666666666, 1)
  expect_equal(unname(c(r$agree)), c(ref), tolerance = 1e-10)
  ref <- 0.47333333333333333
  expect_equal(unname(c(r$mean)), c(ref), tolerance = 1e-10)
})


test_that("amscal() equals the reference Aldrich-McKelvey fit", {
  set.seed(2)
  tr <- c(-1.2, -0.4, 0.1, 0.9, 1.5)
  Z <- t(vapply(1:12, function(i) {
    4 + rnorm(1, 0, 0.3) + runif(1, 0.6, 1.4) * tr + rnorm(5, 0, 0.3)
  }, numeric(5)))
  Z[2, 3] <- NA
  r <- amscal(Z)
  # rmorie 1.5.0: R$morie_spatial_voting_aldrich_mckelvey(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-1.3968015311364701, -0.45576260986428757, -0.004343515119547359,
        0.68642693052414538, 1.1704807255961596)
  expect_equal(unname(c(r$zhat)), c(ref), tolerance = 1e-10)
  ref <- c(3.9818675053048178, 4.1642301973753497, 3.8618934197277142,
        3.8342990500177208, 4.5484643854903286, 3.8536561007526453,
        4.1079499252499394, 4.09562704340303, 4.0964582332536414,
        3.8351169155428777, 4.003858747746583, 4.1488743037932343)
  expect_equal(unname(c(r$alpha)), c(ref), tolerance = 1e-10)
  ref <- c(1.1103787041988711, 1.1696842713418745, 1.5660408231527345,
        0.74185197666821434, 0.70269181199822639, 1.0320456821491473,
        0.8332978965919039, 1.2634291105775803, 1.5342437616247422,
        1.2869976327734034, 0.38809750227701056, 1.3441422435984427)
  expect_equal(unname(c(r$beta)), c(ref), tolerance = 1e-10)
  ref <- c(1.0271059666710176, 1.081963918862505, 1.448595751546836,
        0.6862168634370196, 0.6499935113173626, 0.95464752161042687,
        0.77080480593462664, 1.1686783734530928, 1.4191833074009716,
        1.190479376733693, 0.35899216972684594, 1.2433384333056021)
  expect_equal(unname(c(r$weights)), c(ref), tolerance = 1e-10)
})


test_that("bbscal() equals the reference Basic Space fit", {
  set.seed(3)
  X <- matrix(rnorm(24), 12) %*% matrix(rnorm(12), 2) +
    matrix(rnorm(72, 0, 0.3), 12)
  X[c(5, 30, 44)] <- NA
  r <- bbscal(X, 2, 5)
  # rmorie 1.5.0: R$morie_spatial_voting_blackbox(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-0.70998337454914395, 0.34946019597299244, 0.78270624351505913,
        -0.59111861623757944, -0.087195100955407168, -0.23893973491552098,
        1.1377034151260779, 1.0623092586779368, -0.6977344825313494,
        0.6508593279941699, -0.081079950573413537, -1.5769871815238217,
        0.067857539659783458, 0.270279154638007, 0.01564631920960162,
        0.41395160849782847, -0.54346628458363955, -0.31129406652879166,
        0.82711687999328132, -0.33428094533054248, 0.65065231578663785,
        -0.95324265324782276, 0.26080152994556449, -0.3640213980399073)
  expect_equal(unname(c(colsgn(r$ideal, ref))), c(ref), tolerance = 1e-10)
  ref <- c(-1.092912986400159, 1.7545954901064482, -0.83277773174660374,
        1.0573839509113865, 1.0994902102241524, 0.43151556980633976,
        0.27328108670824508, -0.47323096605039405, -1.1851001823101246,
        -0.047356143614779869, -0.25508093156693762, 1.0952303278770372)
  expect_equal(unname(c(colsgn(r$weights, ref))), c(ref), tolerance = 1e-10)
  ref <- c(55.946332588245276, 8.8204895435381658)
  expect_equal(unname(c(r$eigen)), c(ref), tolerance = 1e-10)
  ref <- 0.98943751712465788
  expect_equal(unname(c(r$fit)), c(ref), tolerance = 1e-10)
  ref <- c(0.35579191905573526, -0.70737693110774391, 0.51937947262822171,
        -0.36762141193096931, -0.46719044780073854, -0.1063540634357983)
  expect_equal(unname(c(r$cmeans)), c(ref), tolerance = 1e-10)
})


test_that("optcls() equals the Optimal Classification it states", {
  # rmorie 1.5.0 (morie_spatial_voting_optimal_classification) puts every
  # legislator on the unit sphere, which in one dimension leaves two
  # positions, and splits each vote at the midpoint of its two sides, so
  # its fit can be worse than voting with the majority. crimlab keeps the
  # mean of the directions and splits each vote where it misclassifies
  # least; this is that algorithm written out in R.
  oc <- function(v, maxit, nstart, seed) {
    set.seed(seed)
    n <- nrow(v)
    best <- Inf
    for (r in seq_len(nstart)) {
      x <- rnorm(n)
      nr <- numeric(ncol(v))
      for (k in seq_len(maxit)) {
        for (j in seq_len(ncol(v))) {
          y <- v[, j]
          d <- mean(x[which(y == 1)]) - mean(x[which(y == 0)])
          if (!is.na(d) && d != 0) nr[j] <- sign(d)
        }
        s <- ifelse(v == 1, 1, -1) * rep(nr, each = n)
        x <- rowMeans(s, na.rm = TRUE)
      }
      e <- 0
      for (j in seq_len(ncol(v))) {
        ok <- !is.na(v[, j])
        p <- x[ok] * nr[j]
        y <- v[ok, j][order(p)]
        p <- sort(p)
        m <- length(p)
        up <- vapply(0:m, function(k) {
          sum(y[seq_len(k)]) + sum(y[seq_len(m) > k] == 0)
        }, 0)
        up <- up[c(TRUE, p[-m] < p[-1], TRUE)]
        e <- e + min(up, m - up)
      }
      if (e < best) {
        best <- e
        fit <- list(ideal = x, normals = nr, errors = e)
      }
    }
    fit
  }
  set.seed(4)
  v <- matrix(rbinom(80, 1, plogis(outer(rnorm(8), rnorm(10)) * 3)), 8)
  v[2, 3] <- NA
  r <- optcls(v, 1, 10, 2, 7)
  ref <- oc(v, 10, 2, 7)
  expect_equal(unname(c(r$ideal)), ref$ideal, tolerance = 1e-12)
  expect_equal(unname(c(r$normals)), ref$normals, tolerance = 1e-12)
  expect_identical(r$errors, as.integer(ref$errors))
  # the errors of voting with the majority, as rmorie counts them
  expect_identical(r$null, 25L)
  expect_equal(r$pre, (25 - ref$errors) / 25, tolerance = testthat_tolerance())
})


test_that("dblctr() equals the reference double centring", {
  set.seed(5)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  r <- dblctr(D)
  # rmorie 1.5.0: R$morie_spatial_voting_double_centering(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.95541009322337256, -1.0342532175759751, 1.2740661301037903,
        -0.085262105230570789, -1.9568549977194867, 0.84689409719886966,
        -1.0342532175759753, 1.9543003340665854, -1.6686724195212987,
        -0.14735843483915731, 1.4554741714993467, -0.55949043362950024,
        1.2740661301037903, -1.6686724195212985, 1.7993877534954024,
        -0.030588322447425657, -2.3796440989621286, 1.0054509573316599,
        -0.085262105230570803, -0.14735843483915745, -0.030588322447425615,
        0.076418602316878692, 0.36495309405716242, -0.17816283385688664,
        -1.9568549977194867, 1.4554741714993467, -2.3796440989621281,
        0.36495309405716214, 4.534405507357695, -2.0183336762325896,
        0.84689409719886999, -0.55949043362950002, 1.0054509573316599,
        -0.17816283385688664, -2.0183336762325896, 0.90364188918844701)
  expect_equal(unname(c(r$B)), c(ref), tolerance = 1e-10)
})


test_that("cmdsf() equals the reference classical scaling", {
  set.seed(5)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  r <- cmdsf(D)
  # rmorie 1.5.0: R$morie_spatial_voting_classical_mds(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.97550941473972774, -0.99845742163080309, 1.2809099300892759,
        -0.10358051211826237, -2.0437277891939911, 0.88934637811405337,
        -0.06157495414148504, 0.97845956035850223, -0.39831834566598728,
        -0.25629998015254962, -0.59798179824634734, 0.33571551784786735)
  expect_equal(unname(c(colsgn(r$coords, ref))), c(ref), tolerance = 1e-10)
  ref <- c(8.5677552691358194, 1.6558089105125635)
  expect_equal(unname(c(r$eigen)), c(ref), tolerance = 1e-10)
  ref <- 2.0651770504185459e-16
  expect_equal(unname(c(r$stress)), c(ref), tolerance = 1e-10)
  ref <- 0.99999999999999978
  expect_equal(unname(c(r$fit)), c(ref), tolerance = 1e-10)
})


test_that("smacf() equals the reference SMACOF fit", {
  set.seed(5)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  r <- smacf(D, maxit = 40)
  # rmorie 1.5.0: R$morie_spatial_voting_smacof(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.39774166480659567, -1.264250114939238, 0.81736884706550483,
        0.20993433832860894, -0.15585841350872284, -0.0049363217527489192,
        0.89224077624019138, -0.59423369426109751, 1.0647243563465856,
        -0.18817814746627251, -2.1239131825945585, 0.94935989173515145)
  expect_equal(unname(c(r$coords)), c(ref), tolerance = 1e-10)
  ref <- 4.9009244786191755e-05
  expect_equal(unname(c(r$stress)), c(ref), tolerance = 1e-10)
  ref <- 40
  expect_equal(unname(c(r$iter)), c(ref), tolerance = 1e-10)
})


test_that("smacf() with weights equals the reference fit", {
  set.seed(5)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  w <- matrix(seq(0.5, 1.5, length.out = 36), 6)
  w <- (w + t(w)) / 2
  r <- smacf(D, maxit = 20, w = w)
  # rmorie 1.5.0: R$morie_spatial_voting_smacof(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.34197650741053448, -1.2148261713840094, 0.71871419355552391,
        0.27318491929774102, -0.051675657312573293, -0.067373791567216873,
        0.90817639245799042, -0.67062933473735753, 1.1422193334692721,
        -0.18398332917436114, -2.1308768054940388, 0.93509374347849361)
  expect_equal(unname(c(r$coords)), c(ref), tolerance = 1e-09)
  ref <- 0.0049825778718769189
  expect_equal(unname(c(r$stress)), c(ref), tolerance = 1e-09)
})


test_that("nmdsf() equals the reference nonmetric scaling", {
  set.seed(5)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  r <- nmdsf(D, maxit = 30)
  # rmorie 1.5.0: R$morie_spatial_voting_nonmetric_mds(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.2601351727808689, -0.49178271496170345, 0.19758197763432603,
        0.18215178267922996, -0.13841158499730769, -0.0096746331354138574,
        0.24495088143877694, -0.18035844987741778, 0.51540608260810483,
        -0.1590401599574213, -0.75081826603076285, 0.32985991181871993)
  expect_equal(unname(c(r$coords)), c(ref), tolerance = 1e-10)
  ref <- 0.0021260064168143024
  expect_equal(unname(c(r$stress)), c(ref), tolerance = 1e-10)
  ref <- 30
  expect_equal(unname(c(r$iter)), c(ref), tolerance = 1e-10)
})


test_that("mdsfit() equals the reference Mardia fit", {
  ev <- c(4, 2, -1, 0.5)
  r <- mdsfit(ev)
  # rmorie 1.5.0: R$morie_spatial_voting_mds_fit_stats(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.53333333333333333, 0.26666666666666666, 0, 0.066666666666666666)
  expect_equal(unname(c(r$fit)), c(ref), tolerance = 1e-10)
  ref <- c(0.53333333333333333, 0.80000000000000004, 0.80000000000000004,
        0.8666666666666667)
  expect_equal(unname(c(r$cumfit)), c(ref), tolerance = 1e-10)
})


test_that("unfstr() equals the reference stress", {
  set.seed(6)
  xr <- matrix(rnorm(8), 4)
  xs <- matrix(rnorm(6), 3)
  Du <- matrix(runif(12), 4)
  Du[2, 2] <- NA
  r <- unfstr(xr, xs, Du)
  # rmorie 1.5.0: R$morie_spatial_voting_unfolding_stress(), the same input,
  # its values printed with 17 significant digits
  ref <- 26.180309032504869
  expect_equal(unname(c(r$stress)), c(ref), tolerance = 1e-10)
})


test_that("mlsmu() equals the reference unfolding", {
  set.seed(7)
  Dm <- matrix(runif(32), 8)
  r <- mlsmu(Dm, maxit = 30, nstart = 2)
  # rmorie 1.5.0: R$morie_spatial_voting_mlsmu6(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.19702564349045115, -0.0062469271168970946, 0.096512604335521945,
        0.020641500146212601, 0.098037554800090659, -0.13943926112969732,
        -0.20841070910955412, -0.058120405416127861, 0.15694714025682976,
        -0.17459432442651268, -0.047834469627150163, -0.21877841468661205,
        -0.31887826023902455, 0.059024973391177654, -0.031390152106662912,
        0.57550350743795498)
  expect_equal(unname(c(r$xr)), c(ref), tolerance = 1e-09)
  ref <- c(0.084707694176971604, -0.0068616903819012504, 0.066676849083221065,
        -0.14452285287829142, -0.28742910419470968, 0.19576589040742981,
        -0.28392620251301376, 0.37558941630029358)
  expect_equal(unname(c(r$xs)), c(ref), tolerance = 1e-09)
  ref <- 2.9461106199472185
  expect_equal(unname(c(r$stress)), c(ref), tolerance = 1e-09)
  ref <- 30
  expect_equal(unname(c(r$iter)), c(ref), tolerance = 1e-09)
})


test_that("smunf() equals the reference SMACOF unfolding", {
  set.seed(7)
  Dm <- matrix(runif(32), 8)
  r <- smunf(Dm, maxit = 30)
  # rmorie 1.5.0: R$morie_spatial_voting_smacof_unfolding(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.61597656290590397, 0.21212535350108197, 0.43290310737559617,
        0.25248895097851348, 0.59852224317032199, 0.35209456593758037,
        0.56275284837813211, 0.50715552028624966, 0.27857900131946739,
        0.095909378347548868, 0.19563373027191805, 0.11409598187469779,
        0.27037580426620234, 0.15926243577354032, 0.25432160416584604,
        0.22940337057109289)
  expect_equal(unname(c(r$xr)), c(ref), tolerance = 1e-10)
  ref <- c(-0.44643643664619437, -0.27742855218340223, -0.53638691223485024,
        -0.50675767520224302, -0.20197221275773347, -0.12534552621461398,
        -0.24244926035238795, -0.22902365397042143)
  expect_equal(unname(c(r$xs)), c(ref), tolerance = 1e-10)
  ref <- 9.7518435687469065
  expect_equal(unname(c(r$stress)), c(ref), tolerance = 1e-10)
})


test_that("ipsrec() passes the respondent coordinates through", {
  set.seed(6)
  xr <- matrix(rnorm(8), 4)
  xs <- matrix(rnorm(6), 3)
  Du <- matrix(runif(12), 4)
  Du[2, 2] <- NA
  r <- ipsrec(xr, xs)
  # rmorie 1.5.0: R$morie_spatial_voting_ideal_point_recovery(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.26960598203752734, -0.62998541403927499, 0.8686598276522649,
        1.7271955171152407, 0.024187641767722562, 0.36802517699251047,
        -1.3092042982085272, 0.73862193071140592)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-10)
})


test_that("nomutl() equals the reference utilities", {
  set.seed(8)
  x <- matrix(rnorm(8), 4)
  zy <- matrix(rnorm(6), 3)
  zn <- matrix(rnorm(6), 3)
  r <- nomutl(x, zy, zn, 10, c(1, 0.5))
  # rmorie 1.5.0: R$morie_spatial_voting_nominate_utility(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.13477395289271824, 0.005891477150315183, 0.37939088874552129,
        0.38215476684984623, 8.6787642187215681, 3.4556136887740712,
        9.4916583474379159, 7.514076167585559, 4.7552586702102779,
        2.3309364017913743, 8.1719437676400517, 9.7322060493220253)
  expect_equal(unname(c(r$uyea)), c(ref), tolerance = 1e-10)
  ref <- c(8.7043881538846843, 7.403294154383147, 8.6657726229780003,
        7.5907553449098462, 7.0481083377340736, 1.8015629955405075,
        6.4753836109265448, 4.0255349882182712, 2.4878925266711933,
        7.9704125596910664, 1.423362367304789, 1.0906621272739767)
  expect_equal(unname(c(r$unay)), c(ref), tolerance = 1e-10)
  ref <- c(0.00018974984304811368, 0.00061246709957542039,
        0.00025186089499537292, 0.00073964479295696498,
        0.83625946808660201, 0.83943776194015551, 0.95330397324490623,
        0.97035996672670877, 0.90613801273102312, 0.0035421386890378179,
        0.99882883143516787, 0.99982341717363377)
  expect_equal(unname(c(r$prob)), c(ref), tolerance = 1e-10)
})


test_that("nompr() equals the reference probability", {
  x <- c(0.1, 0.2)
  r <- nompr(x, c(0.3, 0), c(-0.3, 0.1))
  # rmorie 1.5.0: R$morie_spatial_voting_nominate_vote_prob(), the same input,
  # its values printed with 17 significant digits
  ref <- 0.65343152359269663
  expect_equal(unname(c(r$prob)), c(ref), tolerance = 1e-10)
})


test_that("nomll() equals the reference log-likelihood", {
  set.seed(8)
  x <- matrix(rnorm(8), 4)
  zy <- matrix(rnorm(6), 3)
  zn <- matrix(rnorm(6), 3)
  vv <- matrix(c(1, 0, 1, 1, 0, 0, 1, NA, 1, 0, 1, 0), 4)
  r <- nomll(vv, x, zy, zn)
  # rmorie 1.5.0: R$morie_spatial_voting_nominate_loglik(), the same input,
  # its values printed with 17 significant digits
  ref <- -52.341511479387705
  expect_equal(unname(c(r$loglik)), c(ref), tolerance = 1e-10)
  ref <- 0.36363636363636365
  expect_equal(unname(c(r$gmp)), c(ref), tolerance = 1e-10)
  ref <- 4
  expect_equal(unname(c(r$correct)), c(ref), tolerance = 1e-10)
  ref <- 11
  expect_equal(unname(c(r$total)), c(ref), tolerance = 1e-10)
})


test_that("procrs() equals the reference rotation", {
  set.seed(9)
  A <- matrix(rnorm(12), 6)
  B <- A %*% matrix(c(0.6, -0.8, 0.8, 0.6), 2) + 0.05 * rnorm(12)
  r <- procrs(A, B)
  # rmorie 1.5.0: R$morie_spatial_voting_procrustes(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-1.4211558016837851, -0.50141448839648439, 0.096909991106336191,
        0.10289122493023745, -0.74367700335552733, -0.37670535748095813,
        0.1476950152922111, -0.64038571606243244, -0.25259652701660146,
        -0.43055797366152626, 1.1455946628832254, -1.2102900249143684)
  expect_equal(unname(c(r$rotated)), c(ref), tolerance = 1e-10)
  ref <- c(0.618980089151414, -0.7854066776098273, 0.78540667760982752,
        0.61898008915141423)
  expect_equal(unname(c(r$rotation)), c(ref), tolerance = 1e-10)
  ref <- 4.8134542624450383
  expect_equal(unname(c(r$scale)), c(ref), tolerance = 1e-10)
  ref <- 0.002791611832415023
  expect_equal(unname(c(r$mse)), c(ref), tolerance = 1e-10)
})


test_that("irtlik() equals the reference likelihood", {
  vv <- matrix(c(1, 0, 1, 1, 0, 0, 1, NA, 1, 0, 1, 0), 4)
  x <- matrix(c(-1, 0.5, 0, 1.2))
  b <- matrix(c(1, -0.5, 2))
  r <- irtlik(vv, x, c(0, 0.2, -0.1), b)
  # rmorie 1.5.0: R$morie_spatial_voting_bayesian_irt_likelihood(),
  # the same input, its values printed with 17 significant digits
  ref <- -17.301215088491542
  expect_equal(unname(c(r$loglik)), c(ref), tolerance = 1e-10)
  ref <- c(0.15865525393145705, 0.69146246127401301, 0.5, 0.88493032977829178,
        0.61791142218895267, 0.32635522028792008, 0.42074029056089696, 0,
        0.028716559816001803, 0.86433393905361733, 0.53982783727702899,
        0.99379033467422384)
  expect_equal(unname(c(r$prob)), c(ref), tolerance = 1e-10)
  ref <- 0.27272727272727271
  expect_equal(unname(c(r$accuracy)), c(ref), tolerance = 1e-10)
})


test_that("irtpst() equals the reference summaries", {
  set.seed(10)
  ch <- array(rnorm(20 * 4 * 2), c(20, 4, 2))
  r <- irtpst(ch)
  # rmorie 1.5.0: R$morie_spatial_voting_bayesian_irt_posterior(),
  # the same input, its values printed with 17 significant digits
  ref <- c(0.12536087849772037, -0.51377702800764369, 0.37353480222527613,
        0.014881347284647281, 0.20313093251790099, 0.24294228030738496,
        -0.23058909721969109, -0.21548411560559486)
  expect_equal(unname(c(r$mean)), c(ref), tolerance = 1e-10)
  ref <- c(0.77880878451931512, 0.86273716437226933, 0.64629386561041635,
        0.97356724916922166, 0.7893632673064801, 0.8160793525542489,
        0.84533565353169859, 0.97150101146684631)
  expect_equal(unname(c(r$sd)), c(ref), tolerance = 1e-10)
  ref <- c(-1.361583034631856, -1.4409010519037322, -0.58922116009508063,
        -1.2683564916991679, -1.3283753333465036, -1.238316833549189,
        -1.4240362025275681, -1.4455782379908069)
  expect_equal(unname(c(r$lo)), c(ref), tolerance = 1e-10)
  ref <- c(1.2793202019249099, 1.1333395180834769, 1.2989923744428444,
        1.2932203633591519, 1.4069010156861252, 1.2215292445147985,
        1.0888066317893998, 1.2577774903995906)
  expect_equal(unname(c(r$hi)), c(ref), tolerance = 1e-10)
})


test_that("ordocl() equals the reference ordered classification", {
  set.seed(11)
  Y <- matrix(sample(1:4, 40, replace = TRUE), 10)
  Y[3, 2] <- NA
  r <- ordocl(Y, 1, 20)
  # rmorie 1.5.0: R$morie_spatial_voting_ordered_oc(), the same input,
  # its values printed with 17 significant digits
  ref <- c(1.3709584471466685, -0.56469817139608869, 0.3631284113373392,
        0.63286260496104041, 0.40426832314099903, -0.10612451609148403,
        1.5115219974389389, -0.094659038413097557, 2.0184237138770418,
        -0.062714099052420993)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-10)
  # one vector of cutpoints per item, as in rmorie
  ref <- c(1.0019105260530865, -0.10039177725221385, -1, -0.10039177725224689,
    -0.5185654640506463, -0.5185654640506463, 0.1003917772522185,
    0.3354113437425833, -0.51856546404915971, -0.33438631009553998)
  expect_equal(unlist(r$cuts), ref, tolerance = 1e-10)
  ref <- 0.28205128205128205
  expect_equal(unname(c(r$correct)), c(ref), tolerance = 1e-10)
  ref <- 3
  expect_equal(unname(c(r$iter)), c(ref), tolerance = 1e-10)
})


test_that("anchv() equals the reference vignette correction", {
  set.seed(12)
  y <- sample(1:5, 10, TRUE)
  V <- matrix(sample(1:5, 30, TRUE), 10)
  V[4, 2] <- NA
  V[5, 1:2] <- NA
  y[6] <- NA
  r <- anchv(y, V)
  # rmorie 1.5.0: R$morie_spatial_voting_anchoring_vignettes(), the same input,
  # its values printed with 17 significant digits
  ref <- c(2, 2, 0, 4, 4, NA, 2, 3, 2, 4)
  expect_equal(unname(c(r$scores)), c(ref), tolerance = 1e-10)
  ref <- c(2, 1, 4, 2, 1, 2, 2, 1, 2, 3, 2, 1, 4, 2, 2.333333333333333, 2, 2, 1,
        2, 3, 4, 4, 5, 5, 3.6666666666666665, 5, 4, 3, 4, 3, 5, 4, 5, 5, 5,
        5, 4, 4, 5, 4)
  expect_equal(unname(c(r$thresh)), c(ref), tolerance = 1e-10)
  ref <- c(0.94280904158206336, 0.87770745147251106, 0.75030857849485033,
        0.5163977794943222)
  expect_equal(unname(c(r$dif)), c(ref), tolerance = 1e-10)
  ref <- c(1, 2, 3)
  expect_equal(unname(c(r$order)), c(ref), tolerance = 1e-10)
})


test_that("indscl() equals the reference INDSCAL fit", {
  set.seed(13)
  D1 <- as.matrix(dist(matrix(rnorm(10), 5)))
  D2 <- as.matrix(dist(matrix(rnorm(10), 5)))
  r <- indscl(list(D1, D2), maxit = 20)
  # rmorie 1.5.0: R$morie_spatial_voting_indscal(), the same input,
  # its values printed with 17 significant digits
  ref <- c(8.2028999136934179e-34, 8.1656051096257012e-37,
        3.6943023677271783e-38, 5.178855824760028e-39,
        1.9218139428718008e-40, -1.1435603671161978e-41,
        -2.9505574644837055e-43, -2.2514067942700399e-44,
        -1.2853643008711673e-45, -2.769819702276587e-48)
  expect_equal(unname(c(colsgn(r$config, ref))), c(ref), tolerance = 1e-10)
  ref <- c(1082151.4044401974, 1070329.6815013324, 1082151.4044401974,
        1070329.6815013324)
  expect_equal(unname(c(r$weights)), c(ref), tolerance = 1e-10)
  ref <- 93.249833256970476
  expect_equal(unname(c(r$stress)), c(ref), tolerance = 1e-10)
})


test_that("nrmvec() equals the reference normal vector", {
  set.seed(14)
  x <- matrix(rnorm(16), 8)
  y <- rnorm(8)
  y[3] <- NA
  r <- nrmvec(x, y)
  # rmorie 1.5.0: R$morie_spatial_voting_normal_vectors(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.40050387876246457, -0.9162950633372533)
  expect_equal(unname(c(r$normal)), c(ref), tolerance = 1e-10)
  ref <- -66.390317849774235
  expect_equal(unname(c(r$degrees)), c(ref), tolerance = 1e-10)
  ref <- 0.088897631298404001
  expect_equal(unname(c(r$r2)), c(ref), tolerance = 1e-10)
  ref <- c(0.47944637203201823, 0.27283115757210058, -0.62419830634442686)
  expect_equal(unname(c(r$coef)), c(ref), tolerance = 1e-10)
})


test_that("cutlns() equals the reference cutting lines", {
  nv <- rbind(c(0.6, 0.8), c(1, 0), c(-0.3, 0.2))
  r <- cutlns(nv, c(0.1, -0.2, 0))
  # rmorie 1.5.0: R$morie_spatial_voting_cutting_lines(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-1, -0.20000000000000001, -1, 0.87499999999999989, -10,
        -1.4999999999999998, 1, -0.20000000000000001, 1, -0.625, 10,
        1.4999999999999998)
  expect_equal(unname(c(r$ends)), c(ref), tolerance = 1e-10)
  ref <- c(53.13010235415598, 0, 146.30993247402023)
  expect_equal(unname(c(r$angles)), c(ref), tolerance = 1e-10)
  ref <- c(0, -0.20000000000000001, 0, 0.12499999999999994, 0, 0)
  expect_equal(unname(c(r$mids)), c(ref), tolerance = 1e-10)
})


test_that("dwnom() equals the reference fit", {
  set.seed(4)
  v <- matrix(rbinom(80, 1, plogis(outer(rnorm(8), rnorm(10)) * 3)), 8)
  v[2, 3] <- NA
  r <- dwnom(v, maxit = 10)
  # rmorie 1.5.0: R$morie_spatial_voting_dw_nominate(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.39814836668668474, 0.26350131947302319, -0.19496557619466176,
        -0.094965772263708309, -0.28776161180205173, -0.19260235617712054,
        0.40029453015737015, 0.10548000631226913, 0.4439700884372425,
        0.32047856291968496, 0.39278527118357337, 0.39258724651120186,
        0.33787243085299257, 0.37097480838591818, 0.35203118589463345,
        0.31614193324419826)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-10)
  ref <- c(-0.9783122842910561, -0.9991225330728668, -0.99716112360037812,
        -0.95128931589520449, -0.96758146693762925, -0.97636583913412223,
        -0.99999803930953446, -0.99969479895542757, -0.9783122842910561,
        0.96758146693762925, 0.20713540113948636, 0.041882739954043478,
        -0.075297367683283864, -0.3082995904370196, -0.25255911157355876,
        -0.21612438125283667, 0.0019802467237147966, 0.02470443161595132,
        0.20713540113948636, 0.25255911157355876)
  expect_equal(unname(c(r$normals)), c(ref), tolerance = 1e-10)
  ref <- c(-0.038259824345575164, -0.058139746586905613, -0.18482382099038586,
        -0.13392520212548933, -0.31692802242298856, -0.12611102685676417,
        0.023664863703928559, -0.064465948727279038, -0.038259824345575164,
        0.31692802242298856)
  expect_equal(unname(c(r$cuts)), c(ref), tolerance = 1e-10)
  ref <- -41.595827535466299
  expect_equal(unname(c(r$loglik)), c(ref), tolerance = 1e-10)
  ref <- 0.759493670886076
  expect_equal(unname(c(r$gmp)), c(ref), tolerance = 1e-10)
})


test_that("nombt() equals the reference bootstrap", {
  set.seed(4)
  v <- matrix(rbinom(80, 1, plogis(outer(rnorm(8), rnorm(10)) * 3)), 8)
  v[2, 3] <- NA
  f <- dwnom(v, maxit = 5)
  r <- nombt(v, f$ideal, f$normals, f$cuts, 4, 3)
  # rmorie 1.5.0: R$morie_spatial_voting_nominate_bootstrap(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.099310462063456506, 0.064584289606832532, 0.046184588370997315,
        0.062073240841644624, 0.084713983470872931, 0.036674157733610685,
        0.10970965573270321, 0.13484323113299243, 0.13738266010299674,
        0.07570779998272327, 0.053782537086450463, 0.12669676866598553,
        0.094307149226298681, 0.13104376439640339, 0.065472462913275506,
        0.057820495654432702)
  expect_equal(unname(c(r$se)), c(ref), tolerance = 1e-10)
  ref <- c(0.47868570004828143, 0.31388506463787491, -0.11310691772773283,
        -0.035958100769631102, -0.096542177027432605, -0.1480247136011873,
        0.50440611516192868, 0.25215377411148882, 0.49976175591707972,
        0.23191776856446383, 0.46106053261045499, 0.43324117443297899,
        0.33149533906796624, 0.38781082424549129, 0.24779021606754464,
        0.1891243129209717)
  expect_equal(unname(c(r$mean)), c(ref), tolerance = 1e-10)
})


test_that("irtem() equals the reference EM fit", {
  set.seed(4)
  v <- matrix(rbinom(80, 1, plogis(outer(rnorm(8), rnorm(10)) * 3)), 8)
  v[2, 3] <- NA
  r <- irtem(v)
  # rmorie 1.5.0: R$morie_spatial_voting_em_irt(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.26483709096739982, 0.69308154521629173, -0.83752523916727939,
        -0.31692647471499974, -0.78350110323011024, -0.35118098148520183,
        0.85941688393090232, 0.40336088465359538)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-09)
  ref <- c(-3.8418224262068823, -4.5780422791938369, -2.6051568238860798,
        -0.38257059334989896, -0.51776204412030236, -0.3598108830319045,
        -4.4170951250150035, -0.95145883494800787, -3.8418224262068823,
        0.51776204412030236)
  expect_equal(unname(c(r$disc)), c(ref), tolerance = 1e-09)
  ref <- c(1.1218044124558195, -0.027977505321687038, 0.86764148011386744,
        -0.68991384553852297, 1.1909173629682641, -0.32641686314713636,
        -2.4398369134692781, -0.026621910646079855, 1.1218044124558195,
        -1.1909173629682641)
  expect_equal(unname(c(r$diff)), c(ref), tolerance = 1e-09)
  ref <- -25.079902609614237
  expect_equal(unname(c(r$loglik)), c(ref), tolerance = 1e-09)
  ref <- 116
  expect_equal(unname(c(r$iter)), c(ref), tolerance = 1e-09)
})


test_that("irtem() in two dimensions equals the reference fit", {
  set.seed(4)
  v <- matrix(rbinom(80, 1, plogis(outer(rnorm(8), rnorm(10)) * 3)), 8)
  v[2, 3] <- NA
  r <- irtem(v, 2, 30, conv = "abs")
  # rmorie 1.5.0: R$morie_spatial_voting_em_irt(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.34844575630714486, 0.57865278067735526, -0.78626522545299937,
        -0.35224040917945632, -0.53079303705405667, -0.35778035146370224,
        0.77834509442534727, 0.26625430259885369, 1.017927689127708,
        -0.27569960335755583, 0.18845225915804137, 0.18855810751831778,
        -0.50304740890012689, 0.0010127434122069547, -0.026430916614526474,
        -0.64451500090088687)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-09)
  ref <- -11.73278034120076
  expect_equal(unname(c(r$loglik)), c(ref), tolerance = 1e-09)
})


test_that("scboot() equals the reference bootstrap", {
  set.seed(2)
  tr <- c(-1.2, -0.4, 0.1, 0.9, 1.5)
  Z <- t(vapply(1:12, function(i) {
    4 + rnorm(1, 0, 0.3) + runif(1, 0.6, 1.4) * tr + rnorm(5, 0, 0.3)
  }, numeric(5)))
  r <- scboot(Z, nboot = 5, seed = 5)
  # rmorie 1.5.0: R$morie_spatial_voting_nonparametric_bootstrap(),
  # the same input, its values printed with 17 significant digits
  ref <- c(0.055653564556787959, 0.12997557285678746, 0.092370380124841625,
        0.039820359295141203, 0.027931849613045383)
  expect_equal(unname(c(r$se)), c(ref), tolerance = 1e-10)
  ref <- c(-1.3505684348454483, -0.49032091207949713, -0.045240687729130848,
        0.68806595044717644, 1.1980640842068999)
  expect_equal(unname(c(r$mean)), c(ref), tolerance = 1e-10)
  ref <- c(-1.4162506060697497, -0.67014437484045897, -0.15743085239004803,
        0.64952619042260906, 1.168914381891909)
  expect_equal(unname(c(r$lo)), c(ref), tolerance = 1e-10)
  ref <- c(-1.2778189305708971, -0.37373264761905767, 0.062738666422013478,
        0.74497619267000559, 1.2382666323092966)
  expect_equal(unname(c(r$hi)), c(ref), tolerance = 1e-10)
})


test_that("wrdfsh() equals the reference Wordfish fit", {
  set.seed(15)
  dtm <- matrix(rpois(80, 5), 8)
  r <- wrdfsh(dtm, 20)
  # rmorie 1.5.0: R$morie_spatial_voting_wordfish(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-0.11621625173092627, 1.4091497138587468, -1.5889446022482421,
        -0.17502380729718459, 1.1597596660507157, -0.18566729430337917,
        0.43380491876176747, -0.93686234309149785)
  expect_equal(unname(c(r$pos)), c(ref), tolerance = 1e-10)
  ref <- c(-0.18930761136255975, 0.18244043950713507, 0.37279242520118838,
        0.41045819283631629, -0.28327607429299989, -0.17082389414646329,
        -0.21534471983946415, -0.057618477885617497, 0.10218411359253651,
        -0.16290119258587599)
  expect_equal(unname(c(r$beta)), c(ref), tolerance = 1e-10)
  ref <- c(-2.2577353590785663, -2.2906397823110192, -2.3755162380573687,
        -2.5603138919916266, -2.5947245807410786, -2.3880414433548447,
        -2.2839061440433119, -2.180171116464602, -2.1040303085846608,
        -2.5134837230122078)
  expect_equal(unname(c(r$alpha)), c(ref), tolerance = 1e-10)
  ref <- 330.71355399127162
  expect_equal(unname(c(r$loglik)), c(ref), tolerance = 1e-10)
})


# Behaviour of the scaling functions, carried over from the tests of the
# methods they port: recovery on simulated data, invariants, refusals.

test_that("amscal() gives z-scored stimuli in the true order", {
  set.seed(1)
  truth <- c(-1.2, -0.4, 0.1, 0.9, 1.5)
  Z <- t(vapply(1:60, function(i) {
    4 + rnorm(1, 0, 0.3) + runif(1, 0.6, 1.4) * truth + rnorm(5, 0, 0.3)
  }, numeric(5)))
  f <- amscal(Z)
  expect_s3_class(f, c("amscal", "svfit"))
  expect_equal(mean(f$zhat), 0, tolerance = 1e-12)
  expect_equal(sd(f$zhat), 1, tolerance = 1e-12)
  expect_identical(order(f$zhat), order(truth))
  expect_gt(cor(f$zhat, truth), 0.99)
  # the stimuli minimise the Aldrich-McKelvey loss over unit, centred z
  loss <- function(z) {
    z <- (z - mean(z)) / sqrt(sum((z - mean(z))^2))
    sum(apply(Z, 1, function(r) sum(lm.fit(cbind(1, r), z)$residuals^2)))
  }
  l0 <- loss(f$zhat)
  for (k in 1:5) {
    e <- replace(numeric(5), k, 1e-3)
    expect_gt(loss(f$zhat + e), l0 - 1e-12)
  }
  expect_error(amscal(matrix(1:4, 4)), "at least two stimuli")
  # one row with a gap, one placing both stimuli at the same point
  expect_error(amscal(matrix(c(1, 1, NA, 1), 2)), "at least one respondent")
})

test_that("bbscal() gives the dimensions asked, sparse rows unscaled", {
  set.seed(1)
  X <- matrix(rnorm(60), 30, 2) %*% matrix(rnorm(16), 2, 8) +
    matrix(rnorm(240, 0, 0.3), 30, 8)
  X[1, 1:6] <- NA
  f <- bbscal(X, ndim = 2)
  expect_identical(dim(f$ideal), c(30L, 2L))
  expect_true(all(is.na(f$ideal[1, ])))
  expect_false(anyNA(f$ideal[-1, ]))
  expect_gt(f$fit, 0.9)
  expect_error(bbscal(matrix(1:3)), "two issues")
  expect_error(bbscal(X, minscale = 99, ndim = 40), "Too few respondents")
})

test_that("classical scaling recovers a configuration from its distances", {
  set.seed(1)
  X <- matrix(rnorm(20), 10)
  D <- as.matrix(dist(X))
  f <- cmdsf(D)
  expect_equal(unname(as.matrix(dist(f$coords))), unname(D), tolerance = 1e-10)
  expect_equal(f$stress, 0, tolerance = 1e-10)
  B <- dblctr(D)$B
  expect_equal(B, t(B), tolerance = testthat_tolerance())
  ev <- eigen(B, symmetric = TRUE)$values
  expect_equal(mdsfit(ev)$fit[1:2], pmax(ev[1:2], 0) / sum(abs(ev)),
               tolerance = testthat_tolerance())
  expect_equal(sum(mdsfit(ev)$fit[1:2]), f$fit,
               tolerance = testthat_tolerance())
  expect_length(mdsfit(c(0, 0))$fit, 0L)
  expect_error(cmdsf(matrix(1:6, 2)), "square")
  expect_error(cmdsf(D, ndim = 11), "at most")
})

test_that("SMACOF stress does not rise with more iterations", {
  set.seed(2)
  D <- as.matrix(dist(matrix(rnorm(30), 10)))
  s <- vapply(c(1, 5, 20, 80), function(m) smacf(D, maxit = m)$stress, 1)
  expect_true(all(diff(s) <= 1e-12))
  expect_true(all(is.finite(s) & s >= 0))
  f <- smacf(D, init = cmdsf(D)$coords)
  expect_lt(f$stress, s[1])
  expect_error(smacf(D, w = matrix(1, 2, 2)), "`w` must be")
  expect_error(smacf(D, init = matrix(0, 3, 2)), "`init` must have")
  g <- nmdsf(D, maxit = 50)
  expect_true(is.finite(g$stress))
  expect_gte(g$stress, 0)
})

test_that("unfolding gives finite stress and passes ideal points on", {
  set.seed(3)
  D <- matrix(runif(40), 10)
  expect_gte(unfstr(matrix(rnorm(20), 10), matrix(rnorm(8), 4), D)$stress, 0)
  u <- smunf(D, maxit = 30)
  expect_true(is.finite(u$stress))
  expect_identical(dim(u$xs), c(4L, 2L))
  expect_identical(ipsrec(u$xr, u$xs)$ideal, u$xr)
  m <- mlsmu(D, maxit = 30, nstart = 2)
  expect_true(is.finite(m$stress))
  expect_error(unfstr(matrix(0, 2, 2), matrix(0, 3, 2), D), "must have a row")
})

test_that("NOMINATE utilities, probabilities and likelihood", {
  set.seed(4)
  x <- matrix(rnorm(8), 4)
  zy <- matrix(rnorm(6), 3)
  zn <- matrix(rnorm(6), 3)
  u <- nomutl(x, zy, zn)
  expect_identical(dim(u$prob), c(4L, 3L))
  expect_true(all(u$prob > 0 & u$prob < 1))
  expect_equal(u$udiff, u$uyea - u$unay, tolerance = testthat_tolerance())
  expect_equal(nompr(c(0, 0), c(1, 0), c(-1, 0))$prob, 0.5)
  expect_equal(nompr(x[1, ], zy[1, ], zn[1, ])$prob, u$prob[1, 1],
               tolerance = testthat_tolerance())
  v <- matrix(rbinom(12, 1, 0.5), 4)
  l <- nomll(v, x, zy, zn)
  expect_true(is.finite(l$loglik))
  expect_gte(l$gmp, 0)
  expect_lte(l$gmp, 1)
  expect_identical(l$total, 12L)
  expect_error(nompr(c(0, 0), 1, 1), "same length")
  expect_error(nomutl(x, zy, zn, w = 1:3), "one weight per dimension")
})

test_that("procrs() recovers an orthogonal rotation", {
  set.seed(5)
  A <- matrix(rnorm(20), 10)
  Q <- matrix(c(cos(1), sin(1), -sin(1), cos(1)), 2)
  f <- procrs(A, A %*% Q)
  expect_equal(f$rotation, Q, tolerance = 1e-12)
  expect_equal(f$mse, 0, tolerance = 1e-20)
  expect_equal(det(f$rotation), 1, tolerance = testthat_tolerance())
  expect_error(procrs(A, A[1:3, ]), "same dimensions")
})

test_that("Wordfish separates two blocs of documents", {
  set.seed(6)
  lam <- rep(c(9, 1, 5), each = 10)
  dtm <- rbind(t(replicate(5, rpois(30, lam))),
               t(replicate(5, rpois(30, rev(lam)))))
  f <- wrdfsh(dtm, maxit = 200)
  expect_true(all(sign(f$pos[1:5]) == -sign(f$pos[6:10])))
  expect_equal(mean(f$pos), 0, tolerance = 1e-10)
  expect_error(wrdfsh(matrix(-1, 2, 2)), "counts from 0")
})

test_that("geometry: normal vectors and cutting lines", {
  set.seed(7)
  x <- matrix(rnorm(40), 20)
  y <- 2 * x[, 1] - x[, 2] + 1
  f <- nrmvec(x, y)
  expect_equal(f$r2, 1, tolerance = testthat_tolerance())
  expect_equal(f$normal, c(2, -1) / sqrt(5), tolerance = testthat_tolerance())
  expect_equal(f$coef, c(1, 2, -1), tolerance = testthat_tolerance())
  expect_equal(f$degrees, f$radians * 180 / pi,
               tolerance = testthat_tolerance())
  g <- cutlns(rbind(c(0, 1), c(1, 0), c(0, 0)), c(0.5, 0.2, 0))
  expect_equal(g$ends[1, ], c(-1, 0.5, 1, 0.5))
  expect_equal(g$ends[2, ], c(0.2, -10, 0.2, 10))
  expect_equal(g$ends[3, ], rep(0, 4), tolerance = testthat_tolerance())
  expect_equal(g$angles, c(90, 0, 0), tolerance = testthat_tolerance())
  expect_error(cutlns(matrix(1, 2, 3), 1:2), "two columns")
})

test_that("IRT likelihood and posterior summaries", {
  v <- rbind(c(1, 0), c(0, 1))
  f <- irtlik(v, matrix(c(1, -1)), c(0, 0), matrix(c(1, -1)))
  expect_equal(f$prob, matrix(pnorm(c(1, -1, -1, 1)), 2),
               tolerance = testthat_tolerance())
  expect_equal(f$loglik, 4 * log(pnorm(1)), tolerance = testthat_tolerance())
  expect_identical(f$accuracy, 1)
  set.seed(8)
  ch <- array(rnorm(200 * 4), c(200, 4, 1))
  p <- irtpst(ch, std = FALSE)
  expect_equal(p$mean[, 1], colMeans(ch[, , 1]),
               tolerance = testthat_tolerance())
  expect_equal(p$lo[, 1], unname(apply(ch[, , 1], 2, quantile, 0.025)),
               tolerance = testthat_tolerance())
  q <- irtpst(ch)
  expect_equal(mean(q$mean), 0, tolerance = 1e-12)
  expect_error(irtpst(1:3), "numeric array")
})

test_that("Optimal Classification beats the null model on clean votes", {
  set.seed(9)
  x <- sort(rnorm(20))
  v <- outer(x, rnorm(30), ">") * 1
  f <- optcls(v, maxit = 20, nstart = 2)
  expect_identical(dim(f$ideal), c(20L, 1L))
  expect_gt(f$pre, 0.8)
  expect_lte(f$errors, f$null)
  g <- ordocl(matrix(sample(1:4, 60, TRUE), 15), ndim = 1, maxit = 20)
  expect_length(g$cuts, 4)
  expect_gte(g$correct, 0)
  expect_lte(g$correct, 1)
})

test_that("anchoring vignettes recode each rating by its own thresholds", {
  f <- anchv(c(2, 4, NA), rbind(c(1, 3, 5), c(1, 3, 5), c(2, 2, NA)), ncat = 3)
  expect_equal(f$thresh[1, ], c(1, 3), tolerance = testthat_tolerance())
  expect_equal(f$scores, c(1, 2, NA), tolerance = testthat_tolerance())
  expect_identical(f$order, c(1L, 2L, 3L))
  expect_error(anchv(1:2, matrix(1, 3, 2)), "one row per rating")
})

test_that("INDSCAL weights follow each individual's stretch", {
  set.seed(10)
  X <- matrix(rnorm(16), 8)
  D1 <- as.matrix(dist(X))
  D2 <- as.matrix(dist(X %*% diag(c(2, 0.5))))
  f <- indscl(list(D1, D2), maxit = 50)
  expect_identical(dim(f$weights), c(2L, 2L))
  expect_true(all(f$weights >= 0.01))
  expect_identical(f$nind, 2L)
  expect_error(indscl(list(D1, D1[1:3, 1:3])), "one size")
  expect_error(indscl(1:3), "list of square")
})

test_that("DW-NOMINATE fit, bootstrap and EM ideal points", {
  set.seed(11)
  x <- rnorm(20)
  v <- matrix(rbinom(600, 1, pnorm(outer(x, rnorm(30)) * 2)), 20)
  f <- dwnom(v, maxit = 20)
  expect_identical(dim(f$ideal), c(20L, 2L))
  expect_true(all(sqrt(rowSums(f$ideal^2)) <= 1 + 1e-12))
  b <- nombt(v, f$ideal, f$normals, f$cuts, nboot = 5)
  expect_true(all(b$se >= 0))
  e <- irtem(v)
  expect_true(e$conv)
  expect_gt(abs(cor(e$ideal[, 1], x)), 0.9)
  # A restart from the fit stays at it: each of three more steps changes
  # the measure of convergence by less than tol = 1e-6. Under the default
  # measure (1 minus the correlation) the points may still move by 4e-3,
  # as rmorie's do on this input.
  s <- irtem(v, start = list(e$diff, e$disc, e$ideal), maxit = 3)
  expect_lt(1 - min(cor(s$ideal, e$ideal)), 3e-6)
  a <- irtem(v, conv = "abs", maxit = 5000)
  expect_true(a$conv)
  s <- irtem(v, start = list(a$diff, a$disc, a$ideal), maxit = 3,
    conv = "abs")
  expect_lt(max(abs(s$ideal - a$ideal)), 3e-6)
  expect_error(irtem(v, xprior = 1), "mean and a positive variance")
  expect_error(irtem(v, start = list(1)), "`start` must")
})

test_that("irtem() reproduces the EM fit of emIRT::binIRT", {
  V <- em_votes()
  r <- irtem(V)
  # emIRT::binIRT with the same starts, priors and control: 83 iterations
  expect_identical(r$iter, 83L)
  expect_true(r$conv)
  t <- irtem(V, maxit = 20000, tol = 1e-13, conv = "abs")
  expect_lt(abs(t$ideal[1, 1] + 0.5391190086), 1e-9)
  expect_lt(max(abs(t$diff[1:2] - c(0.6150891743, -0.4241460078))), 1e-9)
  expect_lt(abs(t$disc[1, 1] + 6.5721230364), 1e-9)
  expect_lt(abs(t$loglik + 267.3624873117), 1e-9)
})

test_that("the two-dimensional EM fixed point is a posterior mode", {
  V <- em_votes()
  r <- irtem(V, ndim = 2, maxit = 20000, tol = 1e-12, conv = "abs")
  lp <- function(a, b, x) {
    M <- matrix(a, nrow(V), ncol(V), byrow = TRUE) + x %*% t(b)
    sum(pnorm(ifelse(V == 1, M, -M), log.p = TRUE)[!is.na(V)]) -
      0.5 * sum(x^2) - 0.5 * (sum(a^2) + sum(b^2)) / 25
  }
  h <- 1e-5
  e <- matrix(0, 30, 2)
  e[1, 2] <- h
  g <- (lp(r$diff, r$disc, r$ideal + e) -
    lp(r$diff, r$disc, r$ideal - e)) / (2 * h)
  expect_lt(abs(g), 1e-5)
})

test_that("bootstrap of scaling positions", {
  set.seed(12)
  Z <- matrix(rnorm(100), 20) + rep(1:5, each = 20)
  f <- scboot(Z, nboot = 10)
  expect_length(f$se, 5)
  expect_true(all(f$lo <= f$hi))
  g <- scboot(Z, "blackbox_t", nboot = 3)
  expect_length(g$se, 20)
  expect_error(scboot(Z, "x"), "must be \"am\"")
})

test_that("print, format and summary describe every part", {
  f <- cmdsf(as.matrix(dist(1:5)), 1)
  expect_output(print(f), "cmdsf:")
  expect_length(format(f), length(f) + 1L)
  s <- summary(f)
  expect_s3_class(s, "summary.svfit")
  expect_identical(s$parts$part, names(f))
  expect_output(print(s), "parts")
  expect_match(format(votagr(matrix(1:4, 2)))[2], "2 x 2 array")
})

test_that("the refusals name the argument", {
  expect_error(votagr("a"), "numeric matrix")
  expect_error(dblctr(matrix(c(1, NA, NA, 1), 2)), "missing values")
  expect_error(optcls(matrix(1, 2, 2), ndim = 0), "whole number")
  expect_error(irtpst(array(1, c(2, 2, 1)), std = NA), "TRUE or FALSE")
  expect_error(mlsmu(matrix(1, 2, 2), tol = -1), "at least 0")
})

test_that("amscal() leaves a one-placement respondent unscaled", {
  set.seed(2)
  tr <- c(-1.2, -0.4, 0.1, 0.9, 1.5)
  Z <- t(vapply(1:12, function(i) {
    4 + rnorm(1, 0, 0.3) + runif(1, 0.6, 1.4) * tr + rnorm(5, 0, 0.3)
  }, numeric(5)))
  # one placement cannot give an intercept and a slope: 0 and 1
  z1 <- Z
  z1[3, -2] <- NA
  a <- amscal(z1)
  expect_identical(c(a$alpha[3], a$beta[3]), c(0, 1))
  # two stimuli placed alike by everyone get one position, and a
  # respondent who placed only those two has nothing to regress on
  z2 <- cbind(Z, Z[, 5])
  z2[4, ] <- NA
  z2[4, 5:6] <- c(3, 4)
  expect_error(
    amscal(z2),
    paste("amscal: a respondent's placements could not be regressed on the",
          "stimulus positions."),
    fixed = TRUE
  )
})

test_that("each scaling kernel refuses arguments it cannot use", {
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  set.seed(7)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  # issue 10 is answered only by respondents with too few answers
  X <- matrix(sample(1:7, 200, TRUE), 20)
  X[1:3, 10] <- 5
  X[4:20, 10] <- NA
  X[1:3, 1:5] <- NA
  ef(bbscal(X, minscale = 8),
     paste("bbscal: Every issue needs at least one response among the",
           "scaled respondents."))
  ef(smacf(D, w = "a"), "smacf: `w` must be a numeric matrix.")
  ef(nmdsf(D, maxit = 0),
     "nmdsf: `maxit` must be one whole number of at least 1.")
  ef(unfstr(matrix(1:4, 2), matrix(1:4, 2), matrix(1, 2, 2), w = "a"),
     "unfstr: `w` must be a numeric matrix.")
  # a respondent with no dissimilarity at all: every start is NaN
  dm <- matrix(runif(12), 3)
  dm[2, ] <- NA
  ef(mlsmu(dm, ndim = 1), "mlsmu: no start reached a finite stress.")
  ef(smunf(D, maxit = 0),
     "smunf: `maxit` must be one whole number of at least 1.")
  x <- matrix(rnorm(6), 3, 2)
  zy <- matrix(rnorm(4), 2)
  zn <- matrix(rnorm(4), 2)
  ef(nomutl(x, zy, zn[1, , drop = FALSE]),
     "nomutl: `zy` and `zn` must have the same rows, and the columns of `x`.")
  ef(nomll(matrix(1, 2, 2), x, zy, zn),
     paste("nomll: `votes` must have a row per row of `x` and a column per",
           "row of `zy`."))
  ef(irtlik(matrix(1, 2, 2), "a", 1, 1),
     "irtlik: `x` must be a numeric matrix.")
  ef(irtlik(matrix(1, 2, 2), matrix(0, 3, 1), 1:2, matrix(0, 2, 1)),
     paste("irtlik: `x` must have a row per legislator, `alpha` and `beta`",
           "one entry (row) per vote, and `beta` the columns of `x`."))
  ef(irtpst(array(c(1, NA, 3, 4), c(2, 2))),
     "irtpst: `chain` must hold finite numbers.")
  ef(indscl(list(D, D), maxit = 0),
     "indscl: `maxit` must be one whole number of at least 1.")
  dn <- D
  dn[1, 2] <- NA
  ef(indscl(list(D, dn)),
     "indscl: `D` must hold finite numbers with no missing values.")
  ef(indscl(list(D, D), ndim = 7),
     "indscl: `ndim` must be at most the number of stimuli, 6.")
  ef(nrmvec(x, 1:2), "nrmvec: `y` must have one value per row of `x`, 3.")
  ef(nrmvec(matrix(1, 3, 1), c(1, 2, NA)),
     paste("nrmvec: `x` with a constant is rank deficient on the rows where",
           "`y` is known."))
  ef(nombt(matrix(1, 2, 2), "a", matrix(0, 2, 1), 1:2),
     "nombt: `x` must be a numeric matrix.")
  ef(nombt(matrix(1, 2, 2), matrix(0, 3, 1), matrix(0, 2, 1), 1:2),
     paste("nombt: `x` must have a row per legislator, `normals` and `cuts`",
           "a row (value) per vote, and `normals` the columns of `x`."))
  set.seed(2)
  ef(scboot(matrix(rnorm(20), 4), nboot = 0),
     "scboot: `nboot` must be one whole number of at least 1.")
  ef(wrdfsh(matrix(1:3, 1)), "wrdfsh: `dtm` must have two documents at least.")
})

test_that("BUG: an unbounded number is not refused as 'at least -inf'", {
  x <- matrix(rnorm(6), 3, 2)
  zy <- matrix(rnorm(4), 2)
  zn <- matrix(rnorm(4), 2)
  e <- expect_error(nomutl(x, zy, zn, beta = "a"), "`beta` must be one number",
                    fixed = TRUE)
  expect_false(grepl("-inf", conditionMessage(e), fixed = TRUE))
  e <- expect_error(nompr(1:2, 1:2, 1:2, beta = "a"))
  expect_false(grepl("-inf", conditionMessage(e), fixed = TRUE))
})

test_that("indscl() takes an array of matrices as a list of them", {
  set.seed(7)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  a <- array(0, c(2, 6, 6))
  a[1, , ] <- D
  a[2, , ] <- D * 1.1
  expect_identical(unclass(indscl(a, ndim = 2)),
                   unclass(indscl(list(D, D * 1.1), ndim = 2)))
})

test_that("irtem() refuses its start, its convergence and singular steps", {
  set.seed(4)
  v <- matrix(rbinom(80, 1, 0.5), 8)
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  ef(irtem(v, start = list(1, 2, 3)),
     paste("irtem: `start` must hold 10 intercepts, 10 by 1 slopes and 8 by",
           "1 ideal points."))
  ef(irtem(v, maxit = 0),
     "irtem: `maxit` must be one whole number of at least 1.")
  ef(irtem(v, conv = "x"), "irtem: `conv` must be \"cor\" or \"abs\".")
  # identical ideal points and a prior of variance 1e300: intercept and
  # slope cannot be told apart
  ef(irtem(v, bprior = c(0, 1e300),
           start = list(rep(0, 10), rep(0, 10), rep(0, 8))),
     "irtem: the M-step of the roll-call parameters is singular.")
  # two dimensions from one roll call with a prior of variance 1e300
  ef(irtem(v[, 1, drop = FALSE], ndim = 2, xprior = c(0, 1e300)),
     "irtem: the M-step of the ideal points is singular.")
})

test_that("the binding refuses arguments of a type no kernel takes", {
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  ef(votagr(factor(c("a", "b"))), "votagr: `x` must not be a factor.")
  ef(votagr(1i), "votagr: `x` must be numbers, strings or a list of them.")
  ef(dynirt(matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1), 3),
            c(NA, "1", "2", "2")),
     "dynirt: `period` must not hold NA.")
  ef(nompr("a", 1, 1), "nompr: `x` must be a numeric vector.")
  ef(nompr(c(1, Inf), 1:2, 1:2),
     "nompr: `x` must hold finite numbers with no missing values.")
  ef(nrmvec(matrix(1:3), c(1, Inf, 2)),
     "nrmvec: `y` must hold finite numbers or NA.")
})

test_that("the printed parts name NULL, lists and short numbers", {
  s <- structure(list(a = NULL, b = list(1, 2), c = c(1.23456, 2)),
                 class = c("x", "svfit"))
  expect_identical(capture.output(print(s)),
                   c("x:", "  a NULL", "  b list of 2", "  c 1.235 2.000"))
})

test_that("dimensions and replicates past their bounds are refused at once", {
  # ndim, nboot or ncat times the points past 2^26 values: refused before
  # any configuration is drawn
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  nd <- "`ndim` times the points scaled must be at most 67108864."
  V <- matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1), 4)
  D <- as.matrix(dist(1:4))
  Y <- matrix(c(1, 2, 3, 1, 2, 3, 3, 2, 1, 1, 3, 2), 4)
  ef(optcls(V, ndim = 1e9), paste("optcls:", nd))
  ef(smacf(D, ndim = 1e9), paste("smacf:", nd))
  ef(nmdsf(D, ndim = 1e9), paste("nmdsf:", nd))
  ef(mlsmu(D[, 1:3], ndim = 1e9), paste("mlsmu:", nd))
  ef(smunf(D[, 1:3], ndim = 1e9), paste("smunf:", nd))
  ef(ordocl(Y, ndim = 1e9), paste("ordocl:", nd))
  ef(indscl(list(D, D), ndim = 1e9), paste("indscl:", nd))
  ef(dwnom(V, ndim = 1e9), paste("dwnom:", nd))
  ef(irtem(V, ndim = 1e9), paste("irtem:", nd))
  ef(scboot(Y, nboot = 1e9),
     "scboot: `nboot` times the positions must be at most 67108864.")
  ef(anchv(c(1, 2, 3, 2), Y, ncat = 1e9),
     "anchv: `ncat` times the ratings must be at most 67108864.")
})

test_that("the vote scalings refuse codes other than 1, 0 and NA", {
  set.seed(1)
  v <- matrix(rbinom(200, 1, 0.5), 20, 10)
  v[3, 4] <- 7
  na <- matrix(NA_real_, 6, 4)
  cd <- "`votes` must hold 1 (yea), 0 (nay) or NA (missing)."
  no <- paste("`votes` must hold at least one vote, 1 (yea) or 0 (nay),",
              "not only NA (missing).")
  for (f in c("irtem", "optcls", "dwnom", "anomin")) {
    expect_error(do.call(f, list(votes = v, ndim = 1)), paste0(f, ": ", cd),
                 fixed = TRUE)
  }
  v[3, 4] <- -1
  expect_error(irtem(v), paste("irtem:", cd), fixed = TRUE)
  for (f in c("irtem", "optcls", "dwnom", "cjrirt")) {
    expect_error(do.call(f, list(votes = na)), paste0(f, ": ", no),
                 fixed = TRUE)
  }
  expect_error(ordocl(na, ndim = 1),
               "ordocl: `Y` must hold at least one response, not only NA",
               fixed = TRUE)
  # a vector is refused at once, not classified as one item of as many
  # categories as values
  t0 <- proc.time()[["elapsed"]]
  expect_error(ordocl(rnorm(1e5)),
               "ordocl: `Y` must be a numeric matrix, respondents by items.",
               fixed = TRUE)
  expect_lt(proc.time()[["elapsed"]] - t0, 5)
})

test_that("a date or a time is refused where a count goes", {
  set.seed(4)
  v <- matrix(rbinom(48, 1, 0.5), 8)
  expect_error(cjrirt(votes = v, nsamp = 3, burn = Sys.Date(), chains = FALSE),
               "^cjrirt: `burn` must be numbers, not dates or times")
  expect_error(cjrirt(votes = v, nsamp = Sys.time(), chains = FALSE),
               "^cjrirt: `nsamp` must be numbers, not dates or times")
})
