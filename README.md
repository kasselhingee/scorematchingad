
# scorematchingad

<!-- badges: start -->
[![R-CMD-check](https://github.com/kasselhingee/scorematchingad/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kasselhingee/scorematchingad/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of `scorematchingad` is to enable fast implementation of score matching estimators through the use of automatic differentiation in the CppAD library. 
Such implementation is best done by either contributing to this package or creating a new package that links to this package. On linux with the `gcc` compiler it is possible to create estimators for new models interactively using `customll()` (*I am pondering how it is that this feature only works on linux with gcc*).

See the file `DESCRIPTION` for a slightly longer description, and `./R/scorematchingad-package.R` (equivalently `help(scorematchingad, scorematchingad)` from within `R`) for an even longer description. The built-in help for `R` packages is well populated. 

## Installation

You can install the development version of scorematchingad from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kasselhingee/scorematchingad")
```

## Example

Some models are already incorporated into `scorematchingad`. Below is an example of estimating the Polynomially-Tilted Pairwise Interaction model (Scealy and Wood, 2023) for compositional data:

``` r
library(scorematchingad)
model <- rppi_egmodel(100)
estalr <- ppi(model$sample,
              paramvec = ppi_paramvec(betap = -0.5, p = ncol(model$sample)),
              trans = "alr")
```

This is an example of obtaining a tape of the score matching discrepancy of a custom likelihood for compositional data, which most naturally lies on the simplex:

``` r
myll <- customll("a1type dirichlet(const veca1 &u, const veca1 &beta) {
  size_t d  = u.size();
  a1type y(0.);  // initialize summation at 0
  for(size_t i = 0; i < d; i++)
  {   y   += beta[i] * log(u[i]);
  }
  return y;
}")

tapes <- buildsmdtape("sim", "identity", "sim", 
 myll, rep(1/3, 3), rep(NA, 3), 
 bdryw="minsq", acut = 0.01)
```

## Using scorematchingad in Other Packages [draft notes]

+ Avoid including the implementations of `wrap` and `as` for `veca1`, `mata1` etc except in `RcppExports.cpp`. This makes sure that the speciailsations definitions are not duplicated for each `.cpp` file in your `./src` directory. In practise you can get the `scorematchingad` types by including just the `_forward.h` header file.

