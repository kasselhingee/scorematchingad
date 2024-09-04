// Code for overriding CppAD's errors with Rcpp::stop()

#include <scorematchingad_forward.h>
#include <Rcpp.h> //for Rcpp::stop()

// A function that handles errors in the format required by CppAD
void Rcpphandler(
   bool known       ,
   int  line        ,
   const char *file ,
   const char *exp  ,
   const char *msg  )
{  // error handler must not return, so throw an exception
   std::ostringstream oss;
   oss << "CppAD Error: " << msg << "\n";
   oss << "Expression: " << exp << "\n";
   oss << "File: " << file << " Line: " << line << "\n";
   if (known) {
       oss << "This error is known and expected.\n";
   }
   Rcpp::stop(oss.str());
}

// function to set the global error handler
// [[Rcpp::export]]
void set_cppad_error_handler() {
    CppAD::ErrorHandler::ErrorHandler info(Rcpphandler);
}


