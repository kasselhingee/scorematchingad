#helper function: given a sample, simulate integer weights, then return weights and a version where each measurement is replicated by the weight size
virtualweights <- function(Y, sizefactor = 1.5){
  ind <- sample(1:nrow(Y), ceiling(sizefactor*nrow(Y)), replace = TRUE)
  weights <- rep(0, nrow(Y))
  weights[as.numeric(names(table(ind)))] <- table(ind)
  newsample <- Y[ind, ]
  return(list(
    newY = newsample,
    w = weights
  ))
}

set.seed(1234)
m <- sec2_3model(1000, maxden = 4)
#simulate weights
set.seed(134)
vw <- virtualweights(m$sample)
acut = 0.1

test_that("smobj, smobjgrad, smobjhess matches for simulated weights", {
  intheta <- cdabyppi:::ppi_cppad_thetaprocessor(m$p)
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

  smoSE_sim <- smestSE(tapes$smotape, m$theta, vw$newY)
  smoSE_direct <- smestSE(tapes$smotape, m$theta, m$sample, w=vw$w)
  expect_equal(smoSE_sim, smoSE_direct, tolerance = 1E-3)
})

test_that("smest() for ppi with minsq match itself and estimatorall1", {
  psphere <- pmanifold("sphere")
  pppi <- ptapell(rep(0.1, m$p), m$theta, llname = "ppi", psphere, fixedtheta = rep(FALSE, length(m$theta)), verbose = FALSE)
  smoppi <- ptapesmo(rep(0.1, m$p), 1:length(m$theta), pll = pppi, pman = psphere, "minsq", acut = acut, verbose = FALSE) #tape of the score function

  out_sim <- smest(smoppi, m$theta * 0 + 1, vw$newY, control = list(tol = 1E-15))
  out_dir <- smest(smoppi, m$theta * 0 + 1, m$sample, control = list(tol = 1E-15), w = vw$w)
  expect_equal(out_sim[names(out_sim) != "counts"], out_dir[names(out_sim) != "counts"], tolerance = 1E-3)

  directestimate <- estimatorall1(m$sample, acut, w = vw$w)

  expect_lt(out_dir$value,
            smobj(smoppi, directestimate$estimator1, m$sample, w = vw$w) + 1E-5 * abs(out_dir$value))

  cdabyppi:::expect_lt_v(abs(out_dir$par - directestimate$estimator1) / out_dir$SE, 1) #proxy for optimisation flatness
  cdabyppi:::expect_lt_v(abs(out_dir$par - m$theta) / out_dir$SE, 3)
})
