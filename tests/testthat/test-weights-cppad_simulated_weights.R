set.seed(1234)
m <- rppi_egmodel(1000, maxden = 4)
#simulate weights
set.seed(134)
vw <- virtualweights(m$sample)
acut = 0.1

test_that("cppad_closed() w = rep(1, nrow(Y)) is near the result as if w omitted", {
  p <- ncol(m$sample)
  tapes <- tape_smd("sim","sqrt", "sph", "ppi",
                        ytape = rep(1/p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        bdryw = "minsq", acut = acut)
  out_constant <- cppad_closed(tapes$smdtape, m$sample, w = rep(1, nrow(m$sample)))
  out_ommit <- cppad_closed(tapes$smdtape, m$sample)

  expect_equal(out_ommit$est, out_constant$est)
  expect_equal(out_ommit$value, out_constant$value)
})

test_that("evaltape_wsum() matches for simulated weights and constant weights", {
  intheta <- ppi_paramvec(m$p)
  tapes <- tape_smd("sim","sqrt", "sph", "ppi",
                        m$sample[1, ], intheta,
                        bdryw = "minsq",
                        acut = acut)
  smd_u <- tape_swap(tapes$smdtape)
  smd_sim <- evaltape_wsum(smd_u, vw$newY, m$theta)
  smd_hardcoded <- evaltape_wsum(smd_u, m$sample, m$theta, w=vw$w)
  expect_equal(smd_sim, smd_hardcoded)

  # compare results to manual calculation
  smd_sim_manual <- sum(vapply(1:nrow(vw$newY), function(i){tapes$smdtape$eval(m$theta, vw$newY[i, ])}, FUN.VALUE = 1.3))
  smd_dir_manual_v <- vapply(1:nrow(m$sample), function(i){tapes$smdtape$eval(m$theta, m$sample[i, ])}, FUN.VALUE = 1.3)
  smd_dir_manual <- sum(smd_dir_manual_v*vw$w)
  expect_equal(smd_sim_manual, smd_sim)
  expect_equal(smd_sim_manual, smd_dir_manual)
})

test_that("evaltape_wsum() matches for simulated weights and constant weights with boundary data", {
  intheta <- ppi_paramvec(m$p)
  tapes <- tape_smd("sim","sqrt", "sph", "ppi",
                        m$sample[1, ], intheta,
                        bdryw = "minsq",
                        acut = acut)
  smd_u <- tape_swap(tapes$smdtape)
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

  smd_sim <- evaltape_wsum(smd_u, vw$newY, pmat = m$theta, xcentres = simYcentres)
  smd_hardcoded <- evaltape_wsum(smd_u, m$sample, pmat = m$theta, xcentres = origYcentres, w=vw$w)
  expect_equal(smd_sim, smd_hardcoded)
})

test_that("cppad_search() for ppi with minsq matches cppad_closed()", {
  tapes <- tape_smd("sim","sqrt", "sph", "ppi",
               ytape = rep(1/m$p, m$p),
               usertheta = ppi_paramvec(m$p),
               bdryw = "minsq",
               acut = acut)

  out_closed <- cppad_closed(tapes$smdtape, vw$newY)
  suppressWarnings({out_dir <- cppad_search(tapes$smdtape, m$theta *0.9, m$sample, control = list(tol = 1E-8, maxit = 500), w = vw$w)})
  expect_equal(as.vector(out_dir$est), 
     out_closed$est,
     tolerance = 1E-1)
})


