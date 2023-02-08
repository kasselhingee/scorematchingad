# add include RcppCommon to cppad.hpp
sed -i "/^# include.*/i # include <RcppCommon.h>" inst/include/cppad/cppad.hpp
# then use it instead of std::cout
sed -i "s/std::cout/Rcpp::Rcout/g" inst/include/cppad/core/*.hpp inst/include/cppad/core/forward/*.omh
sed -i "s/std::cerr/Rcpp::Rcerr/g" inst/include/cppad/core/*.hpp
sed -i "s/std::cerr/Rcpp::Rcerr/g" inst/include/cppad/*/*.hpp
sed -i "s/ cerr/ Rcerr/g" inst/include/cppad/*/*.hpp
# the above is sufficient for avoiding check notes on 20230208, I suspect the rest of CPPAD isn't used yet

# to do it for everything:
sed -i "0,/^# include.*/ //i # include <RcppCommon.h>" inst/include/cppad/*.hpp
sed -i "s/std::cout/Rcpp::Rcout/g" `find inst/include/cppad/ -type f \( -iname "*.hpp" -o -iname "*.omh" \)`
sed -i "s/std::cerr/Rcpp::Rcerr/g" `find inst/include/cppad/ -type f \( -iname "*.hpp" -o -iname "*.omh" \)`
sed -i "s/ cout/ Rcout/g" `find inst/include/cppad/ -type f \( -iname "*.hpp" -o -iname "*.omh" \)`
sed -i "s/ cerr/ Rcerr/g" `find inst/include/cppad/ -type f \( -iname "*.hpp" -o -iname "*.omh" \)`


