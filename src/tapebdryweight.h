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
//' For `tape_bdryw_inbuilt()`, currently available unnormalised log-density functions are:
//'
//' ```{r, results = "asis", echo = FALSE}
//' cat(paste(" +", bdrywnames), sep = "\n")
//' ```
//' @export
// [[Rcpp::export]]
pADFun tape_bdryw_inbuilt(std::string name, veca1 x, const double & acut);
#endif
