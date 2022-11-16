test_that("PPI model with sqrt transformation, minsq divergence weight, acut of 0.1 is quadratic", {
  sqrtman <- pmanifold("sphere")
  ppitape <- tapell(llname = "ppi",
                    xtape = c(0.2, 0.3, 0.5),
                    usertheta = ppi_paramvec(p = 3), 
                    pmanifoldtransform = sqrtman)
  ppismotape <- tapesmo(lltape = ppitape,
                        pmanifoldtransform = sqrtman,
                        divweight = "minsq",
                        acut = 0.1,
                        verbose = FALSE)
  ppismotapeJ <- pTapeJacobian(ppismotape, attr(ppismotape, "xtape"), attr(ppismotape, "dyntape"))
  pForward0(ppismotapeJ, 
    ppi_paramvec(p = 3, AL=0, bL=0, beta=c(-0.1,-0.1,0.5)),
    c(0.1, 0.1, 0.8)) 
  
  objs <- buildsmotape("sphere", "ppi", c(0.2, 0.3, 0.5),
                       weightname = "minsq", acut = 0.1,
                       ppi_paramvec(p = 3))
  pForward0(objs$smotape,
    ppi_paramvec(p = 3, AL=0, bL=0, beta=c(-0.1,-0.1,0.5)),
    c(0.1, 0.1, 0.8))
  pJacobian(objs$smotape,
            ppi_paramvec(p = 3, AL=0, bL=0, beta=c(-0.1,-0.1,0.5)),
            c(0.1, 0.1, 0.8))
})