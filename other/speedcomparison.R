# compare speed

beta = c(-0.3, -0.1, 3)
n = 1000
set.seed(134)
utabl <- MCMCpack::rdirichlet(n, beta+1)

smfunA <- function(utable, beta){
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

smfunB <- function(beta, utabl){
  smofun <- ptapesmo(c(1,1,1,3,3,3), 3)
  ztabl <- sqrt(utabl)
  sc_perpt <- lapply(1:n, function(i){
    scobj <- psmo_n_grad(smofun, ztabl[i, ], beta)
    return(scobj)
  })
  scmo <- mean(unlist(sc_perpt))
  return(scmo)
}

system.time(smfunA(beta, utabl), gcFirst = TRUE)
system.time(smfunB(beta, utabl), gcFirst = TRUE)
