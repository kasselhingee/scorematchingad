
test_that("cppad-based Score2 estimate leads to a match for large number of observations", {
  smofun <- ptapesmo_simplex(c(1,1,1,3,3,3), 3)
  beta = c(-0.3, -0.1, 3)
  n = 1000
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  smobj <- function(beta, utabl){
    sc_perpt <- lapply(1:n, function(i){
      scobj <- psmo(smofun, utabl[i, ], beta)
      return(scobj)
    })
    scmo <- mean(unlist(sc_perpt))
    return(scmo)
  }

  out <- optim(par = beta*0,
               fn = function(beta){smobj(beta,utabl)},
               # gr = function(beta){smobj(beta, utabl)[["grad"]]},
               method = "BFGS")

  expect_equal(beta, out$par, tolerance = 1E-3, ignore_attr = TRUE)
}
