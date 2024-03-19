#ifndef Rcout_handler
#define Rcout_handler

// for replacing cppad errors with Rcpp::Rcout

void Rcout_handler(
   bool known       ,
   int  line        ,
   const char *file ,
   const char *exp  ,
   const char *msg  );

#endif
