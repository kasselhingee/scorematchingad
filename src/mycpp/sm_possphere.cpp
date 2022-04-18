// code for various tools for the positive quadrant of the sphere
namespace Spos { // begin the empty namespace

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> toS(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     // for (int i=0; i<x.size(); i++){
     //   out[i] = CppAD::sqrt(x[i]);
     // }
     out = x.cwiseSqrt();
     return(out);
  }

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromS(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.cwiseProduct(x);
     return(out);
  }

  template <class Type>
  Type logdetJ_fromS(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &z){
     Type out;
     out = z.array().log().sum();
     return(out);
  }


  // manifold tangent-plane projection matrix P (for isometric(?) embeddings this is closely related to the manifold metric
  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmat_S(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);
    Pmat = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n) - x*x.transpose();
    return(Pmat);
  }

  //partial derivative of the tangent-plane projection matrix
  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> dPmat_S(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const int &d){
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, 1> basisvec(n);
    basisvec.setZero();
    basisvec(d) = 1;
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bvx(n, n);
    bvx = -basisvec * x.transpose();
    bvx += bvx.transpose().eval(); //eval() means the tranposition happens in a temporary location
    return(bvx);
  }


////////////////////APPROX HELPERS/////////////////////
  //automatically choose approximation centre
  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> taylorapprox_bdry(
		  CppAD::ADFun<Type> &f,
		  const size_t order,
		  const Eigen::Matrix<Type, Eigen::Dynamic, 1> xbeta,
		  double shiftsize=1E-5){
     Eigen::Matrix<Type, Eigen::Dynamic, 1> x(xbeta.size() / 2);
     x << xbeta.block(0,0,x.size(), 1);
     Eigen::Matrix<Type, Eigen::Dynamic, 1> shiftdir(x.size());
     shiftdir.setOnes();
     shiftdir *= -1 * shiftsize;
     shiftdir = Pmat_S(x) * shiftdir;

     Eigen::Matrix<Type, Eigen::Dynamic, 1> centre(xbeta.size());
     centre = xbeta;
     centre.block(0,0,x.size(), 1) << x + shiftdir;
     std::cout << "Approximation centre is:" << shiftdir.transpose() << std::endl;
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(0);
     out = taylorapprox(f, centre, order, xbeta);
     return(out);
  }

}
