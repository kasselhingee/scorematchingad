# include Rcout_handler.h

void myhandler(
   bool known       ,
   int  line        ,
   const char *file ,
   const char *exp  ,
   const char *msg  )
{  // error handler must not return, so throw an exception
   Rcpp::Rcout << "Rcout: " << line;
}

