# to do it for everything:
# add include RcppCommon
sed -i "0,/^# include.*/s//# include <RcppCommon.h>\n&/" inst/include/cppad/*.hpp

# replace things: 
sed -i "s/std::cout/Rcpp::Rcout/g" `find inst/include/cppad/ -type f \( -iname "*.hpp" -o -iname "*.omh" \)`
sed -i "s/std::cerr/Rcpp::Rcerr/g" `find inst/include/cppad/ -type f \( -iname "*.hpp" -o -iname "*.omh" \)`
sed -i "s/ cout/ Rcout/g" `find inst/include/cppad/ -type f \( -iname "*.hpp" -o -iname "*.omh" \)`
sed -i "s/ cerr/ Rcerr/g" `find inst/include/cppad/ -type f \( -iname "*.hpp" -o -iname "*.omh" \)`

# add include RcppCommon to any file that uses Rcpp
sed -i "0,/^# include.*/s//# include <RcppCommon.h>\n&/" `grep -l "Rcpp" inst/include/cppad/utility/*.hpp`

