test_that("Solution without boundary considerations for PPI has zero gradient and matches numerical minimum", {
  set.seed(13411)
  mod <- ppi_egmodel(100)
  Y <- mod$sample

  Ralr <- pmanifold("Ralr")
  ppitape <- tapell(llname = "ppi",
                  xtape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)), 
                  pmanifoldtransform = Ralr)
  smotape <- tapesmo(lltape = ppitape,
                      pmanifoldtransform = Ralr,
                      divweight = "ones",
                      verbose = FALSE)

  estobj <- cppad_closed(smotape, Y)

  grads <- t(apply(Y, MARGIN = 1, function(x) pJacobian(smotape, estobj$est, x))) 
  totalgrad <- colSums(grads)
  expect_lt(sum(totalgrad^2), 1E-20)

  numericalmin <- ppi(Y, paramvec = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)), trans = "alr", method = "closed")
  expect_equal(numericalmin$est$paramvec, c(estobj$est, tail(mod$beta, 1)), ignore_attr = TRUE)
  expect_equal(numericalmin$SE$paramvec, c(estobj$SE, 0))
})

test_that("Closed-from solution with boundary points matches hard-coded version", {
  mnongamma <- ppi_egmodel(1)
  theta <- ppi_paramvec(beta = c(-0.95, -0.9, 0.5), AL = mnongamma$AL, bL = 0)
  set.seed(1234)
  Ycts <- rppi(1000, paramvec = theta)
  dsample <- round(Ycts * 100)/ 100
  dsample[, 3] <- 1 - rowSums(dsample[, 1:2])

  isbdry <- simplex_isboundary(dsample)
  Yapproxcentres <- dsample
  Yapproxcentres[!isbdry, ] <- NA 
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(dsample[isbdry, , drop = FALSE])

  Ralr <- pmanifold("Ralr")
  ppitape <- tapell(llname = "ppi",
                  xtape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3, bL = 0, betap = tail(theta, 1)), 
                  pmanifoldtransform = Ralr)
  smotape <- tapesmo(lltape = ppitape,
                      pmanifoldtransform = Ralr,
                      divweight = "ones",
                      verbose = FALSE)

  estobj <- cppad_closed(smotape, Y = dsample, Yapproxcentres, approxorder = 10)

  est_hardcode <- ppi(dsample, paramvec = ppi_paramvec(p = 3, bL = 0, betap = tail(theta, 1)),
      trans = "alr", method = "hardcoded")
  expect_equal(est_hardcode$est$paramvec, t_fu2t(estobj$est, attr(ppitape, "usertheta")))

})

test_that("Closed-form solution with all boundary points and alr matches hardcoded", {
  set.seed(123)
  m <- ppi_egmodel(100)
  #add some zeroes
  pushtozero <- function(x){
    whichmin <- which.min(x)
    x[whichmin] <- 0
    x <- x / sum(x) #normalise
    return(x)
  }
  newsample <- t(apply(m$sample, MARGIN = 1, pushtozero))
  mean(apply(newsample, 1, min) == 0) #100% have a zero

  hardcoded <- ppi(Y = newsample,
                   paramvec = ppi_paramvec(bL = 0, p = m$p, betap = tail(m$beta, 1)),
                   method = "hardcoded",
                   trans = "alr")
  cppad <- ppi(Y = newsample,
                   paramvec = ppi_paramvec(bL = 0, p = m$p, betap = tail(m$beta, 1)),
                   method = "closed",
                   trans = "alr")
})

