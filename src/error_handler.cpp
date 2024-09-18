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
   if (known) {
       oss << "CppAD error from a known source: ";
   } else {
       oss << "CppAD Error: " ;
   }
   oss << msg << "\n";
   oss << "Expression: " << exp << "\n";
   oss << "File: " << file << " Line: " << line << "\n";
   Rcpp::stop(oss.str());
}

// function to set the global error handler
// [[Rcpp::export]]
void set_cppad_error_handler() {
  static CppAD::ErrorHandler RcppErrors(Rcpphandler);
}

// the following should use Rcpphandler if the global initiation worked
// [[Rcpp::export]]
void test_Rcpphandler(){
      CppAD::ErrorHandler::Call(
         true     , // reason for the error is known
         __LINE__ , // current source code line number
         __FILE__ , // current source code file name
         "1 > 0"  , // an intentional error condition
         "Testing ErrorHandler"     // reason for error
      );
}

