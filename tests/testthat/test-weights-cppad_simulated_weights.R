set.seed(1234)
m <- ppi_egmodel(1000, maxden = 4)
#simulate weights
set.seed(134)
vw <- virtualweights(m$sample)
acut = 0.1

test_that("cppad_closed() w = rep(1, nrow(Y)) is near the result as if w omitted", {
  tapes <- buildsmotape("sphere", "ppi",
                        ytape = rep(1/p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        weightname = "minsq", acut = acut)
  out_constant <- cppad_closed(tapes$smotape, m$sample, w = rep(1, nrow(m$sample)))
  out_ommit <- cppad_closed(tapes$smotape, m$sample)

  expect_equal(out_ommit$est, out_constant$est)
  expect_equal(out_ommit$value, out_constant$value)
})

test_that("cppad_search() w = rep(1, nrow(Y)) is near the result as if w omitted", {
  tapes <- buildsmotape("sphere", "ppi",
                        ytape = rep(1/p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        weightname = "minsq", acut = acut)
  out_constant <- cppad_search(tapes$smotape, m$theta * 0 + 1, m$sample, control = list(tol = 1E-12), w = rep(1, nrow(m$sample)))
  out_ommit <- cppad_search(tapes$smotape, m$theta * 0 + 1, m$sample, control = list(tol = 1E-12))

  expect_equal(out_ommit$est, out_constant$est)
  expect_equal(out_ommit$value, out_constant$value)
})

test_that("tape_eval_wsum() matches for simulated weights and constant weights", {
  intheta <- ppi_paramvec(m$p)
  tapes <- buildsmotape("sphere", "ppi",
                        m$sample[1, ], intheta,
                        weightname = "minsq",
                        acut = acut)
  smo_u <- tapeSwap(tapes$smotape)
  smo_sim <- tape_eval_wsum(smo_u, vw$newY, m$theta)
  smo_hardcoded <- tape_eval_wsum(smo_u, m$sample, m$theta, w=vw$w)
  expect_equal(smo_sim, smo_hardcoded)

  # compare results to old smobj  
  expect_equal(smobj(smofun = tapes$smotape, theta = m$theta, utabl = vw$newY), smo_sim / nrow(vw$newY))
  expect_equal(smobj(smofun = tapes$smotape, theta = m$theta, utabl = m$sample, w=vw$w), smo_hardcoded/sum(vw$w))

  # compare results to manual calculation
  smo_sim_manual <- sum(vapply(1:nrow(vw$newY), function(i){pForward0(tapes$smotape$ptr, m$theta, vw$newY[i, ])}, FUN.VALUE = 1.3))
  smo_dir_manual_v <- vapply(1:nrow(m$sample), function(i){pForward0(tapes$smotape$ptr, m$theta, m$sample[i, ])}, FUN.VALUE = 1.3)
  smo_dir_manual <- sum(smo_dir_manual_v*vw$w)
  expect_equal(smo_sim_manual, smo_sim)
  expect_equal(smo_sim_manual, smo_dir_manual)
})

test_that("tape_eval_wsum() matches for simulated weights and constant weights with boundary data", {
  intheta <- ppi_paramvec(m$p)
  tapes <- buildsmotape("sphere", "ppi",
                        m$sample[1, ], intheta,
                        weightname = "minsq",
                        acut = acut)
  smo_u <- tapeSwap(tapes$smotape)
  Y <- m$sample
  isbdry <- simplex_isboundary(Y, 1E-2)
  Yapproxcentres <- Y 
  Yapproxcentres[!isbdry, ] <- NA
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(Y[isbdry, , drop = FALSE], shiftsize = 1E-3)
  origYcentres <- Yapproxcentres
  
  Y <- vw$newY 
  isbdry <- simplex_isboundary(Y, 1E-2)
  Yapproxcentres <- Y 
  Yapproxcentres[!isbdry, ] <- NA
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(Y[isbdry, , drop = FALSE], shiftsize = 1E-3)
  simYcentres <- Yapproxcentres

  smo_sim <- tape_eval_wsum(smo_u, vw$newY, pmat = m$theta, xcentres = simYcentres)
  smo_hardcoded <- tape_eval_wsum(smo_u, m$sample, pmat = m$theta, xcentres = origYcentres, w=vw$w)
  expect_equal(smo_sim, smo_hardcoded)
})

test_that("cppad_search() for ppi with minsq matches itself", {
  tapes <- buildsmotape("sphere", "ppi",
               ytape = rep(1/m$p, m$p),
               usertheta = ppi_paramvec(m$p),
               weightname = "minsq",
               acut = acut)

  out_sim <- cppad_search(tapes$smotape, m$theta *0 + 1, vw$newY, control = list(tol = 1E-12))
  out_dir <- cppad_search(tapes$smotape, m$theta *0 + 1, m$sample, control = list(tol = 1E-12, maxit = 1000), w = vw$w)
  expect_equal(out_sim[!(names(out_sim) %in% c("counts", "SE"))], 
     out_dir[!(names(out_sim) %in% c("counts", "SE"))],
     tolerance = 1E-3)

  expect_equal(out_dir[c("est", "value", "sqgradsize")], out_sim[c("est", "value", "sqgradsize")], tolerance = 1E-4)
})


