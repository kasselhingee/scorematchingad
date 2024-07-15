This update is purely to update the internals of the package. It has been checked on windows, mac and linux on old, current and devel versions of R.

Currently parts of the external library `CppAD` receive compilation warnings like "use of bitwise '&' with boolean operands" on mac. The developer of `CppAD` has chosen bitwise operation due to speed (see <https://cppad.readthedocs.io/latest/cmake.html#clang>). So I am hoping CRAN is comfortable with leaving the bitwise operations in the source code?

Summary of changes:

+ Upgraded to a recent version of the external library `CppAD`.
+ Moved to an install method via a `configure` that is more faithful to the recommended installation methods for `CppAD`.
+ Replaced use of `ellipsis` with `rlang`

