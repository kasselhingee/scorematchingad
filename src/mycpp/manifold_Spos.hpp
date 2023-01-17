#ifndef mycpp_manifold_Spos
#define mycpp_manifold_Spos
// code for various tools for the positive quadrant of the sphere
#include <RcppEigen.h>
namespace mantran {
template <typename Type>
struct Spos : public manifold<Type> {
  ~Spos(){};
  Spos(){};

  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.cwiseSqrt();
     return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.cwiseProduct(x);
     return(out);
  }

  Type logdetJfromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &z) override {
     Type out;
     out = z.array().log().sum() + 0.6931472 * z.size(); //final number here is log(2)
     return(out);
  }


  // manifold tangent-plane projection matrix P (for isometric(?) embeddings this is closely related to the manifold metric
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);
    Pmat = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n) - x*x.transpose();
    return(Pmat);
  }

  //partial derivative of the tangent-plane projection matrix
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> dPmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const int &d) override {
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, 1> basisvec(n);
    basisvec.setZero();
    basisvec(d) = 1;
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bvx(n, n);
    bvx = -basisvec * x.transpose();
    bvx += bvx.transpose().eval(); //eval() means the tranposition happens in a temporary location
    return(bvx);
  }


};
}
#endif
