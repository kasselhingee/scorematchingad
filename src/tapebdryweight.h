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
//' @details
//' For `tape_bdryw_inbuilt()`, currently available functions are:
//'
//' ```{r, results = "asis", echo = FALSE}
//' cat(paste(" +", bdrywnames), sep = "\n")
//' ```
//' See [`ppi()`] for details on prodsq and minsq. "ones" is the function that always returns `1`, so is the function to use for manifolds without boundary.
//' @export
// [[Rcpp::export]]
pADFun tape_bdryw_inbuilt(std::string name, veca1 x, const double & acut);
#endif
