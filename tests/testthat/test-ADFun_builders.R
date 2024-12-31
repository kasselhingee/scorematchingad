test_that("Jacobian, Hessian, GradOffset, Swap and LogJacDet all produce new tapes", {
maninfo <- manifoldtransform("sim", "sqrt", "sph")
ppitape <- tapell(ll = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran) 
expect_s4_class(tape_Jacobian(ppitape), "Rcpp_ADFun")
expect_s4_class(tape_Hessian(ppitape), "Rcpp_ADFun")
expect_s4_class(tape_gradoffset(ppitape), "Rcpp_ADFun")
expect_condition(tape_logJacdet(ppitape), class = "Rcpp::exception", regexp = "equal")
expect_s4_class(tape_logJacdet(tape_Jacobian(ppitape)), "Rcpp_ADFun")
expect_s4_class(tape_swap(ppitape), "Rcpp_ADFun")
})

