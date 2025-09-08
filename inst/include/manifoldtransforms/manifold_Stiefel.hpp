#ifndef mantrans_manifold_Stiefel
#define mantrans_manifold_Stiefel

// code for various tools for the Stiefel manifold with Frobenius norm as metric and embedded in Euclidean space using vec()
#include <RcppEigen.h>
namespace mantran {
template <typename Type>
struct Stiefel : public manifold<Type> {
  int nrow;
  int ncol;
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> commutation_mat;
  ~Stiefel(){};
  Stiefel(int nrow_, int ncol_) : nrow(nrow_), ncol(ncol_) {
      commutation_mat = build_commutation_mat(nrow_, ncol_);
  };

  std::string name() const {
    std::string out = "Stiefel";
    return(out);
  }

  // for internal use: the commutation matrix builder
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> build_commutation_mat(int m, int p) {
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> K =
        Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Zero(m*p, m*p);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            K(i*p + j, j*m + i) = Type(1);
    return K;
  }
  
  // for internal use: the Kronecker product
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>
  kronecker(const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> & A, const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> & B) {
      Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> 
          K( A.rows() * B.rows(), A.cols() * B.cols());
      for (int i = 0; i < A.rows(); i++)
          for (int j = 0; j < A.cols(); j++)
              K.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i,j) * B;
      return K;
  }


  // manifold tangent-plane projection matrix P. Shoud match result of R function Stiefel_projmat(invvec(x, m))
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
    // build needed identity matrices
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Im = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(this->nrow, this->nrow);
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Ip = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(this->ncol, this->ncol);

    // build matrix out of x (invvec)
    Eigen::Map<const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> Xmat(x.data(), this->nrow, this->ncol);

    // first term of proj matrix
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> term1 =  kronecker(Ip, Im - Type(0.5) * (Xmat * Xmat.transpose()));
    
    // second term or proj matrix
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> term2 =  Type(0.5) * kronecker(Xmat.transpose(), Xmat) * this->commutation_mat;
    
    return term1 - term2;
  }

  //partial derivative of the tangent-plane projection matrix
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> dPmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const int &d) override {
    //temporary dPmatfun to test compiling
    return Pmatfun(x);
  }

};
}
#endif
