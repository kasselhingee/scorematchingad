# include <iostream>        // standard input/output
# include <vector>          // standard vector
# include <cppad/example/cppad_eigen.hpp>  //load eigen
# include <cppad/cppad.hpp> // the CppAD package
# include <Rcpp.h>
using namespace Rcpp;

typedef std::vector<double> svecd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecd; //a vector of double values
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matd;//a matrix of double
typedef CppAD::AD<double> a1type;   // for first (outer) level of taping
typedef Eigen::Matrix<a1type, Eigen::Dynamic, 1> veca1; //a vector of a1type values
typedef Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> mata1;//a matrix of a1types
typedef CppAD::AD<a1type> a2type;  // for second (inner) level of taping
typedef Eigen::Matrix<a2type, Eigen::Dynamic, 1> veca2;

template <typename T>
struct manifold {
  Eigen::Matrix<T, Eigen::Dynamic, 1> (*toM)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &); //map from simplex to manifold
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> (*Pmatfun)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &); //projection matrix for manifold
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> (*dPmatfun)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &, const int &); //elementwise derivative of projection matrix for manifold
  Eigen::Matrix<T, Eigen::Dynamic, 1> (*fromM)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &); //transformation from manifold to simplex
  T (*logdetJfromM)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &); //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
};
