# include <iostream>        // standard input/output
# include <vector>          // standard vector
# include <cppad/example/cppad_eigen.hpp>  //load eigen
# include <cppad/cppad.hpp> // the CppAD package
# include <Rcpp.h>
using namespace Rcpp;

typedef std::vector<double> svecd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecd; //a vector of a1type values
typedef CppAD::AD<double> a1type;   // for first (outer) level of taping
typedef Eigen::Matrix<a1type, Eigen::Dynamic, 1> veca1; //a vector of a1type values
typedef Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> mata1;//a matrix of a1types
typedef CppAD::AD<a1type> a2type;  // for second (inner) level of taping
typedef Eigen::Matrix<a2type, Eigen::Dynamic, 1> veca2;
