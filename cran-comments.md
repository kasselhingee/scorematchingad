Dear CRAN Team,

Please find attached my an R package for Hyvarinen score matching estimators by automatic differentiation. There are a few quirks that you might spot when looking through the package:

0. The RcppEigen package is listed in the Enhances section of DESCRIPTION because it is needed for calling Rcpp::cppFunction() (via RcppXPtrUtils::cppXPtr) with an RcppEigen dependency.

1. The feature involving customll() appears to only work on linux + gcc combinations. I have a provided a function for users to test if the feature works, and I hope to expand the functionality in the future.

2. Although the source tar ball is only 700KB, the installed size of the package is large due to incorporating the full CppAD source code. I have looked into avoiding this but opted otherwise because:
  + it will be very rare for CppAD to be installed on users computers, and it looks hard to check
  + the CppAD version inside the TMB R package is unsuitable for my purposes: the CppAD version is from 2015 (which misses at least one crucial feature I'm using) and contains modifications bespoke to TMB.


