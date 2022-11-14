# ifndef scorecompdir_TYPES
# define scorecompdir_TYPES

# include <RcppCommon.h>
# include <iostream>        // standard input/output
# include <vector>          // standard vector
# include <cppad/example/cppad_eigen.hpp>  //load eigen
# include <cppad/cppad.hpp> // the CppAD package
# include <cppad/utility/index_sort.hpp> //for index sorting - for Rivest model


using namespace Rcpp;

typedef std::vector<double> svecd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecd; //a vector of double values
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matd;//a matrix of double
typedef CppAD::AD<double> a1type;   // for first (outer) level of taping
typedef Eigen::Matrix<a1type, Eigen::Dynamic, 1> veca1; //a vector of a1type values
typedef Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> mata1;//a matrix of a1types
typedef CppAD::AD<a1type> a2type;  // for second (inner) level of taping
typedef Eigen::Matrix<a2type, Eigen::Dynamic, 1> veca2;

typedef CppAD::ADFun<double> ADFun_double; //required for RCPP_MODULE exporting. See RcppAnnoy for an example

// template <typename T>
// struct manifold {
//   Eigen::Matrix<T, Eigen::Dynamic, 1> (*toM)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &); //map from simplex to manifold
//   Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> (*Pmatfun)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &); //projection matrix for manifold
//   Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> (*dPmatfun)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &, const int &); //elementwise derivative of projection matrix for manifold
//   Eigen::Matrix<T, Eigen::Dynamic, 1> (*fromM)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &); //transformation from manifold to simplex
//   T (*logdetJfromM)(const Eigen::Matrix<T, Eigen::Dynamic, 1> &); //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
// };

template <typename T>
struct manifold { //exactly like a class, but with default public members https://stackoverflow.com/questions/24196885/can-c-struct-have-member-functions
  //make these members 'pure virtual' means this is an 'abstract class'
  virtual Eigen::Matrix<T, Eigen::Dynamic, 1> toM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &) = 0; //map from simplex to manifold
  virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &) = 0; //projection matrix for manifold
  virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dPmatfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &, const int &) = 0; //elementwise derivative of projection matrix for manifold
  virtual Eigen::Matrix<T, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &) = 0; //transformation from manifold to simplex
  virtual T logdetJfromM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &) = 0; //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
  //for taping will need to pass copies - so that the coefficients of the tape are not updated by other calls part way through
  virtual ~manifold(){}; //destructor
  manifold(){};
};

namespace Rcpp {
  template <> veca1 as( SEXP );
  template <> SEXP wrap(const veca1&);
}

# endif

