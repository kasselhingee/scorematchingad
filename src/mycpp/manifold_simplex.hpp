#ifndef mycpp_manifold_simplex
#define mycpp_manifold_simplex
#include <RcppEigen.h>
namespace mantran {
template <typename T>
struct simplex : public manifold<T> {
  ~simplex(){};
  simplex(){};

  Eigen::Matrix<T, Eigen::Dynamic, 1> toM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x) override {
    Eigen::Matrix<T, Eigen::Dynamic, 1> out(x.size());
    out = x;
    return(out);
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x) override {
    Eigen::Matrix<T, Eigen::Dynamic, 1> out(x.size());
    out = x;
    return(out);
  }


  // manifold tangent-plane projection matrix P (for isometric(?) embeddings this is closely related to the manifold metric
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x) override {
    int n = x.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);
    Eigen::Matrix<T, Eigen::Dynamic, 1> ones(n);
    ones.setOnes();
    double nd = n;
    Pmat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n) - (ones*ones.transpose()/nd);
    return(Pmat);
  }

  //partial derivative of the tangent-plane projection matrix
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dPmatfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x, const int &d) override {
    int n = x.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bvx(n, n);
    bvx.setZero();
    return(bvx);
  }
};
}
#endif
