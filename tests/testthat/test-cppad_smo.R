
test_that("cppad-based Score2 function leads to a match with direct estimators for Dirichlet distribution", {
  beta = c(-0.3, -0.1, 3)
  n = 10
  set.seed(134)
  utabl <- MCMCpack::rdirichlet(n, beta+1)

  smobj <- function(beta, utabl){
    ztabl <- sqrt(utabl)
    sc_perpt <- lapply(1:n, function(i){
      scobj <- smo_n_grad(ztabl[i, ], beta)
      args = paste(c(ztabl[i, ], beta), sep = " ", collapse = " ")
      grad = NA #as.numeric(unlist(strsplit(trimws(out[grep("Gradient:", out) + 1]), " +")))
      hess = NA #as.numeric(unlist(strsplit(trimws(out[grep("Hessian:", out) + 1:length(beta)]), " +")))
      return(list(grad = grad,
                  hess = hess,
                  scmo = scobj))
    })
    grad <- NA #colMeans(do.call(rbind, lapply(sc_perpt, "[[", "grad")))
    hess <- NA #apply(simplify2array(lapply(sc_perpt, "[[", "hess")), MARGIN = c(1,2), FUN = mean)
    scmo <- mean(vapply(sc_perpt, "[[", "scmo", FUN.VALUE = -0.1))
    return(list(
     grad = grad,
     hess = hess,
     scmo = scmo))
  }

# There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = beta*0,
        fn = function(beta){smobj(beta,utabl)[["scmo"]]},
        # gr = function(beta){smobj(beta, utabl)[["grad"]]},
        method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimator2_dir(utabl, 1)
  expect_equal(directestimate, out$par, tolerance = 1E-3, ignore_attr = TRUE)
})

test_that("ADFun XPtr for computing values", {
  smofun <- ptapesmo(c(1,1,1,3,3,3), 3) #tape of the score function
  beta = c(-0.3, -0.1, 3)
  n = 10
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

  smobjgrad <- function(beta, utabl){
    sc_perpt <- lapply(1:n, function(i){
      scobj <- psmograd(smofun, utabl[i, ], beta)
      return(scobj)
    })
    scmo <- colMeans(do.call(rbind, sc_perpt))
    return(scmo)
  }

  # There are better optimisers than below: John Nash at https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/)
  out <- optim(par = beta*0,
               fn = function(beta){smobj(beta,utabl)},
               gr = function(beta){smobjgrad(beta, utabl)},
               method = "BFGS")

  # memoisation could be used to avoid calling the smobj function again for gradient computation
  directestimate <- estimator2_dir(utabl, 1)
  expect_equal(directestimate, out$par, tolerance = 1E-3, ignore_attr = TRUE)
})
