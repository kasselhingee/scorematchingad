test_that("custom_uld() can generate a working tape", {
  utape = rep(0.2, 5)
  dyntape = rep(-0.1, 5)
  dirichlet <- tape_uld(system.file("demo_custom_uld.cpp", package = "scorematchingad"),
                        utape, dyntape)
  expect_equal(dirichlet$tape$name, "dirich")
  expect_equal(dirichlet$tape$xtape, utape)
  expect_equal(dirichlet$tape$dyntape, dyntape)
  expect_equal(dirichlet$tape$size_order, 0)
  expect_equal(dirichlet$tape$domain, 5)

  # try further taping
  maninfo <- manifoldtransform("sim", "sqrt", "sph")
  dirichwrtsph <- reembed(dirichlet$tape, maninfo$tran)
  dirichwrtsph <- fixdynamic(dirichwrtsph, dirichwrtsph$dyntape, c(rep(FALSE, 4), TRUE))
  hardwired <- tapell(ll = "dirichlet",
                  ytape = rep(0.2, 5),
                  usertheta = c(NA, NA, NA, NA, dyntape[5]), 
                  tranobj = maninfo$tran) 
  newu <- c(0.1, 0.1, 0.2, 0.2, 0.4)
  newbeta <- c(-0.9, -0.5, -0.1, -0.1, -0.1)
  expect_equal(dirichwrtsph$eval(sqrt(newu), newbeta[-5]), hardwired$eval(sqrt(newu), newbeta[-5]))
  expect_equal(dirichwrtsph$Jac(sqrt(newu), newbeta[-5]), hardwired$Jac(sqrt(newu), newbeta[-5]))
})


if (!requireNamespace("RcppEigen", quietly = TRUE)){skip("Need RcppEigen package")}
test_that("customll_test() returns TRUE or FALSE", {
  skip_on_cran()
  suppressWarnings({out <- customll_test()})
  expect_true(out %in% c(TRUE, FALSE))
})

# warning: for interactive testing load_all() with install doesn't include scorematchingad.h properly

test_that("customll errors correctly with wrong signature", {
  expect_error({dirll <- customll("
a1type dirichlet(const vecd &u, const veca1 &beta) {
        size_t d  = u.size();
        a1type y(0.);  // initialize summation
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
}")})
})

test_that("a tape_uld() gets all the way to the correct score matching estimate", {
  set.seed(13411)
  mod <- rppi_egmodel(100)
  Y <- mod$sample
  ll <- tape_uld("a1type customppi(const veca1 &u, const veca1 &theta) {
          a1type y;
          y = ll::ll_ppi(u, theta);
          return y;
  }", c(0.2, 0.3, 0.5), mod$theta)

  tapes_custom <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = ll$tape,
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
     bdryw = "ones",
     verbose = FALSE)
  tapes <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
     bdryw = "ones",
     verbose = FALSE)

  est_custom <- cppad_closed(tapes_custom$smdtape, Y) 
  est <- cppad_closed(tapes$smdtape, Y) 
  expect_equal(est_custom$est, est$est)
  expect_equal(est_custom$covar, est$covar)
})

