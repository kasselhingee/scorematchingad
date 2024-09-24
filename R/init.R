.onLoad <- function(libname, pkgname){
  # When compiled with debugging, CppAD's errors can be turned in R errors by this initiation
  set_cppad_error_handler()

  # for now needed to make cppad_module load fully (no lazy load),
  # which is needed to have C++ functions return the pADFun class
  # without Error in .getClassesFromCache(Class) : 
  # class should be either a character-string name or a class definition
  # and with the manifolds module working
  cppad_module$empty()
}
