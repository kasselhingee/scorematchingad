test_that("Jacobian, Hessian, GradOffset, Swap and LogJacDet all produce new tapes", {
maninfo <- manifoldtransform("sim", "sqrt", "sph")
ppitape <- tapell(ll = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran) 
expect_s4_class(tapeJacobian(ppitape), "Rcpp_ADFun")
expect_s4_class(tapeHessian(ppitape), "Rcpp_ADFun")
expect_s4_class(tapeGradOffset(ppitape), "Rcpp_ADFun")
expect_condition(tapeLogJacDet(ppitape), class = "Rcpp::exception", regexp = "equal")
expect_s4_class(tapeLogJacDet(tapeJacobian(ppitape)), "Rcpp_ADFun")
expect_s4_class(tapeSwap(ppitape), "Rcpp_ADFun")
})

