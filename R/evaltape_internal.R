#' @name evaltapes_internal
#' @title Advanced: Evaluate `CppAD` Tapes via their Pointer
#' @description The recommended method for evaluating tapes is [`evaltape()`].
#' Internally, `evaltape()` and other methods are using the methods documented here.
#' There methods access the tapes using `Rcpp::XPtr` objects and perform evaluations a single point at a time. 
#' @family tape evaluators
NULL

