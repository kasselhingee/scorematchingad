a <- smobj_sum_b(smofun = smotape, utabl= Y, w = w, theta = theta)
b <- tape_eval_wsum(smotape, xmat = theta, pmat = Y, w = w)
d <- tape_eval_wsum(smofun_u, xmat = Y, pmat = theta, w = w)

smobjgrad_sum_b(smofun = smotape, utabl= Y, w = w, theta = theta)
tape_eval_wsum(Jsmofun_u, xmat = Y, pmat = theta, w = w)

  smoobjold <- function(atheta){
    a <- smobj_sum_b(smofun = smotape, utabl = Y, w = w, theta = atheta)
    attr(a, "normaliser") <- NULL
    return(a)
  }

  smoobjnew <- function(atheta){
    tape_eval_wsum(smotape, pmat = Y, xmat = atheta, w = w)
  }

out <- Rcgmin::Rcgmin(par = theta,
                        fn = smoobj,
                        gr = smograd,
                        control = list(trace = 3, maxit = 10))

out2 <- Rcgmin::Rcgmin(par = theta,
                        fn = smoobjold,
                        gr = smograd,
                        control = control)

out3 <- Rcgmin::Rcgmin(par = theta,
                        fn = smoobjnew,
                        gr = smograd,
                        control = list(trace = 3, maxit = 10))

# looks like out2 and out3 were identical, but `out` is different
# the only difference is the fn component
# noticeably after 10 evaluations of the gradient


