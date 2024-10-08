# to do it for everything:

# replace std::cout and std::cerr: 
sed -i "s/std::cout/Rcpp::Rcout/g" `find inst/cppad/ -type f \( -iname "*.hpp" \)`
sed -i "s/std::cerr/Rcpp::Rcerr/g" `find inst/cppad/ -type f \( -iname "*.hpp" \)`
sed -i "s/ cout/ Rcout/g" `find inst/cppad/ -type f \( -iname "*.hpp" \)`
sed -i "s/ cerr/ Rcerr/g" `find inst/cppad/ -type f \( -iname "*.hpp" \)`

# replace std::exit(1) with Rcpp::stop
sed -i 's/std::exit(1)/Rcpp::stop("An error occurred")/g' `find inst/cppad/ -type f \( -iname "*.hpp" \)`

# add include RcppCommon to any file that uses Rcpp and doesn't already have RcppCommon in it
sed -i "0,/^# include.*/s//# include <RcppCommon.h>\n&/" `grep -L "<.*RcppCommon" $(grep -l "Rcpp" inst/cppad/* -r)`


