# Version 0.1.4
Custom boundary weight functions enabled through `tape_bdryw()`, the results of which can be passed to `tape_smd()`.

# Version 0.1.3
Bug fixes: export of `evaltape_range()` and `keeprange()`.
Added a new tape builder that averages range `avgrange()`.

# Version 0.1.2
Exposing new functions for modifying the arguments of Rcpp_adfun objects: fixdynamic() and fixindependent().

# Version 0.1.1
No longer using RcppXPtrUtils.
Small fix to creation of XPtr to ADFun in movetoXPtr()

# Version 0.1.0
Major version number to reflect that functionality of the package is now pretty good.

# Version 0.0.71
Cleaning up of help.

# Version 0.0.70
+ Looking for `cmake` in the `/Applications/CMake.app/Contents/bin` means installation works on mac.
+ Help for `Rcpp_ADFun` objects.
+ Removed old `customll()` associated functions/code
+ More stable test of `cppad_search()` due to differing convergence of `optimx::Rcgmin()` on differing OS.

# Version 0.0.69
__Major improvement and breaking change:__ taping of custom log-density functions occurs at the *same* time as compilation. This works across all OS (previous `customll()` only worked on some cases of unix) and means the major object passed between functions is a 'tape' object. Use `tape_uld()` for this functionality.
 + New function `tape_uld()`
 + New class wrapping `ADFun` is from an Rcpp module and contains a special object (`pADFun` class) that points to an `ADFun`.
   + This new class has many of the properties available of `ADFun`.
   + It took a lot of experimentation to work out how to get expose this, and then further expose it to external R packages.

Converted all cases of std::exit to Rcpp::stop. See `cout2Rcout.sh` for script.

# Version 0.0.68
CppAD errors handled when compiled in debug mode

# Version 0.0.67
changed source code to use && or || for boolean objects

# Version 0.0.66
+ Deleted many unneeded files from original CppAD source.
+ Rearranged calling of wrapas.hpp. Now much of the types are declared in `scorematchingad_forward.h`. Then, following RcppEigen, `scorematchingad.h` includes `scorematchingad_forward.h`, then `Rcpp.h` then `wrapas.hpp`

# Version 0.0.65
Moved to a recent version of `CppAD` (version 20240000.5), which contains some macros defined in `./inst/cppad/include/cppad/configure.hpp` via a call to `cmake`. Thus this package builds via a configure script for installation. The `CppAD` source code is now in `./inst/cppad`. After `cmake` the `./inst/cppad/include/cppad` directory is moved to `./inst/include/cppad/`. One `cpp` file is also moved to `./src/`

Also moved to using `rlang` instead of `ellipsis` to `check_dots_used()`

# Version 0.0.64

URLs to CppAD and Eigen help pages, data locations and even paper locations are breaking frequently enough for CRAN testing to pick up. So I've removed a few urls in this version.

# Version 0.0.63

No change.

# Version 0.0.62

+ Bug for M1Mac test not floating point but actually because sign eigenvectors of Bingham A matrix are undetermined.

# Version 0.0.61

+ Bug fix to passing `paramvec_start` to `ppi_robust()`
+ More generous limit for equality for a M1Mac equality error
+ Faster examples

# Version 0.0.60
Again, tiny polish of description.

# Version 0.0.59
Small polish of DESCRIPTION and examples in response to CRAN review.

# Version 0.0.58
Small clean up of package version number and spelling in help.

