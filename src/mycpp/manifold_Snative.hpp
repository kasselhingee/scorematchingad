#ifndef mycpp_manifold_Snative
#define mycpp_manifold_Snative

// code for various tools for the sphere, without transformation required
#include <RcppEigen.h>
namespace mantran {
template <typename Type>
struct Snative : public manifold<Type> {
  ~Snative(){};
  Snative(){};

  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x;
     return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x;
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
