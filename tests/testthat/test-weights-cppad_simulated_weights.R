set.seed(1234)
m <- ppi_egmodel(1000, maxden = 4)
#simulate weights
set.seed(134)
vw <- virtualweights(m$sample)
acut = 0.1

test_that("w = rep(1, nrow(Y)) is near the result as if w omitted", {
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, m$p), m$theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, m$p), 1:length(m$theta), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out_constant <- cppadest(smoppi, m$theta * 0 + 1, m$sample, control = list(tol = 1E-12), w = rep(1, nrow(m$sample)))
  out_ommit <- cppadest(smoppi, m$theta * 0 + 1, m$sample, control = list(tol = 1E-12))

  expect_equal(out_ommit$par, out_constant$par)
  expect_equal(out_ommit$value, out_constant$value)
})

test_that("smobj, smobjgrad, smobjhess matches for simulated weights and constant weights", {
  intheta <- cdabyppi:::ppi_paramvec(m$p)
  tapes <- buildsmotape("sphere", "ppi",
                        m$sample[1, ], intheta,
                        weightname = "minsq",
                        acut = acut)
  smo_sim <- smobj(tapes$smotape, m$theta, vw$newY)
  smo_direct <- smobj(tapes$smotape, m$theta, m$sample, w=vw$w)
  expect_equal(smo_sim, smo_direct)

  smograd_sim <- smobjgrad(tapes$smotape, m$theta, vw$newY)
  smograd_direct <- smobjgrad(tapes$smotape, m$theta, m$sample, w=vw$w)
  expect_equal(smograd_sim, smograd_direct)

  smohess_sim <- smobjhess(tapes$smotape, m$theta, vw$newY)
  smohess_direct <- smobjhess(tapes$smotape, m$theta, m$sample, w=vw$w)
  expect_equal(smohess_sim, smohess_direct)
})

test_that("cppadest() for ppi with minsq match itself and estimatorall1", {
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, m$p), m$theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, m$p), 1:length(m$theta), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out_sim <- cppadest(smoppi, m$theta * 0 + 1, vw$newY, control = list(tol = 1E-12))
  out_dir <- cppadest(smoppi, m$theta * 0 + 1, m$sample, control = list(tol = 1E-12), w = vw$w)
  expect_equal(out_sim[!(names(out_sim) %in% c("counts", "SE"))], 
     out_dir[!(names(out_sim) %in% c("counts", "SE"))],
     tolerance = 1E-3)

  directestimate <- estimatorall1(m$sample, acut, w = vw$w)

  expect_lt(out_dir$value,
            smobj(smoppi, directestimate$estimator1, m$sample, w = vw$w) + 1E-5 * abs(out_dir$value))

  expect_lt_v(abs(out_dir$par - directestimate$estimator1) / out_dir$SE, 1E-3) #proxy for optimisation flatness
  expect_lt_v(abs(out_dir$par - m$theta) / out_dir$SE, 3)
})

