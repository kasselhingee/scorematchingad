This update is purely to update the internals of the package. It has been checked on windows, mac and linux on old, current and devel versions of R.

I've fixed compilation warnings like "use of bitwise '&' with boolean operands" on mac.

Summary of changes:

+ Upgraded to a recent version of the external library `CppAD`.
+ Moved to an install method via a `configure` that is more faithful to the recommended installation methods for `CppAD`.
+ Replaced use of `ellipsis` with `rlang`

