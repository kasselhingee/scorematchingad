#ifndef mycpp_manifold_Hclr
#define mycpp_manifold_Hclr
// code for various tools for the additive log ratio transform
#include <RcppEigen.h>
#include <cppad/cppad.hpp> //for CppAD::log
namespace mantran {//names space for manifold-transformation pair (triplets}
template <typename Type>
struct Hclr : public manifold<Type> {
  ~Hclr(){};
  Hclr(){};

  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.array().log(); //log all elements of x
     Eigen::Matrix<Type, 1, 1> sumlog; //use a matrix so that -= is a known operation
     sumlog << out.mean(); //sum logged values - mean would work just as well, but sum has fewer operations (except maybe when dimensions are very large?)
     out -= sumlog; //take the sumlog away from each element
     
     //make toM non-generate in enclosing R^p space by adding a vector perpendicular to the space
     Eigen::Matrix<Type, Eigen::Dynamic, 1> unitones(x.size());
     unitones.setConstant(1.0).normalize();
     out += unitones * (unitones.transpose() * x - 1/sqrt(x.size()));
     return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.array().exp(); //exp all elements of x
     Type sumexp = out.sum(); 
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out2(x.size());
     out2 = out / sumexp; //normalise by sum
     
     //make fromM non-degenerate in enclosing R^p space by adding a vector perpendicular to the space
     Eigen::Matrix<Type, Eigen::Dynamic, 1> unitones(x.size());
     unitones.setConstant(1.0).normalize();
     out2 += unitones * (unitones.transpose() * x);


     return(out2);
  }
  
  Type logdetJfromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &z) override {
    Eigen::Matrix<Type, Eigen::Dynamic, 1> u(z.size());
    u = fromM(z);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> unotp(u.size() - 1);
    unotp = u.head(u.size() - 1);
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> jac(u.size() - 1, u.size() - 1);
    jac = unotp.asDiagonal() + unot * (u.tail(1) - unotp.transpose());
    Type out;
    out = CppAD::log(jac.determinant());
    return(out);
  }


  //Pmat is the same as for simplex (both planes with normal of (1,1,1,....1)
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
}//namesspace mantran for manifold-transformation pair (triplets}

#endif

