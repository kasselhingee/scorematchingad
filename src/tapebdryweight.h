#ifndef tapebdryweight
#define tapebdryweight

# include "scorematchingad_forward.h"
# include "utils/pADFun.h"

CppAD::ADFun<double> tapeh2(veca1 z,
                            a1type (*h2fun)(const veca1 &, const double &),
                            const double & acut);

//' @rdname tape_bdryw
//' @name tape_bdryw
//' @param name Name of an inbuilt function. See details.
//' @param acut The \eqn{a_c} threshold used by the "minsq" and "prodsq" built-in functions. See details.
//' @details
//' For `tape_bdryw_inbuilt()`, currently available functions are:
//' * The function "ones" applies no weights and should be used whenever the manifold does not have a boundary.
//' \deqn{h(x)^2 = 1.}
//' * The function "minsq" is the minima-based boundary weight function \insertCite{@Equation 12, @scealy2023sc}{scorematchingad}
//' \deqn{h(x)^2 = \min(x_1^2, x_2^2, ..., x_p^2, a_c^2).}{h(x)^2 = min(x1^2, x2^2, ..., xp^2, a_c^2),}
//' where \eqn{x_j}{xj} is the jth component of x. 
//' * The function "prodsq" is the product-based \insertCite{@Equation 9, @scealy2023sc}{scorematchingad}
//' \deqn{h(x)^2 = \min(\prod_{j=1}^{p} x_j^2, a_c^2).}{h(x)^2 = min(x1^2 * x2^2 * ... * xp^2, a_c^2),}
//' where \eqn{x_j}{xj} is the jth component of x.
//'
//' The "minsq" and "prodsq" functions are useful when the manifold is the positive orthant, the p-dimensional unit sphere restricted to the positive orthant, or the unit simplex.
//' \insertCite{@scealy2023sc}{scorematchingad} prove that both "minsq" and "prodsq" can be used for score matching the PPI model on the simplex or positive orthant of the sphere.
//' @export
// [[Rcpp::export]]
pADFun tape_bdryw_inbuilt(std::string name, veca1 x, const double & acut);
#endif
