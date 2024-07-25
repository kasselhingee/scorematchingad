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

