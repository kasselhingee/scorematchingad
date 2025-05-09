#ifndef tapell_h
#define tapell_h

# include <RcppEigen.h>
# include <Rcpp.h>
# include <cppad/cppad.hpp>
# include "scorematchingad_forward.h"
# include "likelihoods/likelihoods.hpp"
# include "utils/PrintFor.hpp"
# include "utils/pADFun.h"

//' @noRd
//' @title Get an XPtr to a named log-density function in source code of package
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to a `llPtr` object of the log-density function. Since `llPtr` is itself a pointer object, we have an XPtr pointing to a pointer that points to a function.
// [[Rcpp::export]]
Rcpp::XPtr<llPtr> getllptr(std::string llname);

//' @rdname tape_uld
//' @name tape_uld
//' @param name Name of an inbuilt function. See details.
//' @details
//' For `tape_uld_inbuilt()`, currently available unnormalised log-density functions are:
//'
//' ```{r, results = "asis", echo = FALSE}
//' cat(paste(" +", llnames), sep = "\n")
//' ```
//' @export
// [[Rcpp::export]]
pADFun tape_uld_inbuilt(std::string name, veca1 x, veca1 theta);

#endif
