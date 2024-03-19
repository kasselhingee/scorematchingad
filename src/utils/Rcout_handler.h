#ifndef Rcout_handler
#define Rcout_handler

# include <cppad/utility/error_handler.hpp>
# include <cstring>

// for replacing cppad errors with Rcpp::Rcout

void Rcout_handler(const bool known, const int  line, const char *file , const char *exp, const char *msg);

#endif
