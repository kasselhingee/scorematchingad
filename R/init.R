.onLoad <- function(libname, pkgname){
  # When compiled with debugging, CppAD's errors can be turned in R errors by this initiation
  set_cppad_error_handler()
}
