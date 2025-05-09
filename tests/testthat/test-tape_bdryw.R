test_that("tape_bdryw_inbuilt() works", {
  tape <- tape_bdryw_inbuilt("minsq", rep(0.2, 5), 0.1)
  expect_equal(tape$xtape, rep(0.2, 5))

  expect_equal(tape$forward(0, rep(0.2, 5)), 0.1^2)
  expect_equal(tape$forward(0, c(0.01, 0.19, 0.2, 0.2, 0.2)), 0.01^2)

  expect_equal(tape$Jacobian(rep(0.2, 5)), rep(0, 5))
  expect_equal(tape$Jacobian(c(0.01, 0.19, 0.2, 0.2, 0.2)), c(2*0.01, rep(0, 4)))
})


# warning: for interactive testing load_all() with install doesn't include scorematchingad.h properly
test_that("tape_bdryw() can generate a working tape", {
  out <- tape_bdryw(
    "a1type myminsq(const veca1 &x){
    veca1 xsq = x.array().square();
    a1type minval(0.1 * 0.1);
    for(size_t i=0;i<x.size();i++){
      minval = CppAD::CondExpLe(xsq[i], minval, xsq[i], minval);
    }
    return(minval);
    }", 
       rep(0.2, 5))
  expect_equal(out$fun(c(0.01, 0.19, 0.2, 0.2, 0.2)), 0.0001)
  expect_equal(out$tape$forward(0, c(0.01, 0.19, 0.2, 0.2, 0.2)), 0.0001)
  expect_equal(out$tape$Jacobian(c(0.01, 0.19, 0.2, 0.2, 0.2)), c(2*0.01, rep(0, 4)))
  expect_equal(out$tape$Jacobian(c(0.2, 0.2, 0.2, 0.2, 0.2)), rep(0, 5))
  expect_equal(out$tape$name, "myminsq")
})

test_that("a tape_bdryw() gets all the way to the correct score matching estimate", {
  set.seed(13411)
  mod <- rppi_egmodel(100)
  Y <- mod$sample
  tapes <- tape_smd(
     start = "sim",
     tran = "sqrt",
     end = "sph",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
     bdryw = "minsq",
     acut = 0.1,
     verbose = FALSE)
  
  myminsq <- tape_bdryw(
    "a1type myminsq(const veca1 &x){
    veca1 xsq = x.array().square();
    a1type minval(0.1 * 0.1);
    for(size_t i=0;i<x.size();i++){
      minval = CppAD::CondExpLe(xsq[i], minval, xsq[i], minval);
    }
    return(minval);
    }", 
    rep(0.2, 3))
  
  tapes_custom <- tape_smd(
    start = "sim",
    tran = "sqrt",
    end = "sph",
    ll = "ppi",
    ytape = c(0.2, 0.3, 0.5),
    usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
    bdryw = myminsq$tape,
    verbose = FALSE)

  est_custom <- cppad_closed(tapes_custom$smdtape, Y) 
  est <- cppad_closed(tapes$smdtape, Y) 
  expect_equal(est_custom$est, est$est)
  expect_equal(est_custom$covar, est$covar)
})

