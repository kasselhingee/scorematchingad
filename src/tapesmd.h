#ifndef tapesmd_h
#define tapesmd_h

# include "scorematchingad_forward.h"
# include <Rcpp.h>
# include <cppad/cppad.hpp>
# include "tapebdryweight.h"
# include "manifoldtransforms/bdryweights.hpp"
# include "manifoldtransforms/manifolds.hpp"
# include "utils/PrintFor.hpp"
# include <utils/pADFun.h>

//' @noRd
//' @title The score matching objective calculator.
//' @param xbetain a concatenated vector of sqrt(x) and beta
//' @param n The dimension of x.
//' @param manifoldname The name of the manifold to transform to
//' @param weightname The name of the weight function to use
//' @param acut The constraint a_c in the weight function
//' @return An RCpp::XPtr object pointing to the ADFun
// [[Rcpp::export]]
pADFun tapesmd(pADFun & uldtape, //uld wrt end manifold of `tran` ie M
               transform<a1type> &tran,
               manifold<a1type> &M,
               std::string weightname, //the weight function h^2 - name of a a1type (*h2fun)(const veca1 &, const double &)
               const double & acut, //the acut constraint for the weight functions
               bool verbose
               );

#endif
