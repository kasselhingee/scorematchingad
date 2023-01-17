#ifndef mycpp_manifold_Ralr
#define mycpp_manifold_Ralr
// code for various tools for the additive log ratio transform
#include <RcppEigen.h>
namespace mantran {//names space for manifold-transformation pair (triplets}
template <typename Type>
struct Ralr : public manifold<Type> {
  ~Ralr(){};
  Ralr(){};

  std::string name() const {
    std::string out = "Ralr";
    return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size() - 1);
     out = x.block(0,0,x.size() - 1, 1) / x[x.size() - 1];
     out = out.array().log();
     return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size() + 1);
     Type one_on_u_p;
     one_on_u_p = x.array().exp().sum() + 1.;
     out << x.array().exp(), 1.;
     out /= one_on_u_p;
     return(out);
  }

  Type logdetJfromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &z) override {
    Eigen::Matrix<Type, Eigen::Dynamic, 1> u(z.size() + 1);
    u = fromM(z);
    Type out;
    out = u.array().log().sum();
    return(out);
  }


  // manifold tangent-plane projection matrix P (for isometric(?) embeddings this is closely related to the manifold metric
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);
    Pmat = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n);
    return(Pmat);
  }

  //partial derivative of the tangent-plane projection matrix
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> dPmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const int &d) override {
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bvx(n, n);
    bvx.setZero();
    return(bvx);
  }
};
}//namesspace mantran for manifold-transformation pair (triplets}

#endif
