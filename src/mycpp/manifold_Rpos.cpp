// code for various tools for the positive quadrant of the Reals
namespace Rpos { // begin the empty namespace

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
    out = -x.array().log();
    return(out);
  }

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
    out = Eigen::exp(-x.array());
    return(out);
  }

  template <class Type>
  Type logdetJ_fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &z){
    // du/dz = [-exp(-z), ditto, ditto,] , ...
    // log(abs(det(du/dz))) = log(abs(-exp(-z) * -exp(-z) * ....))
    // = log(abs(-1^p)) + log(exp(-z)) + ...
    // = 0 - z - z - ...
    Type out;
    out = -z.array().sum();
    return(out);
  }


  // manifold tangent-plane projection matrix P (for isometric(?) embeddings this is closely related to the manifold metric
  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmat_M(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);
    Pmat = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n);
    return(Pmat);
  }

  //partial derivative of the tangent-plane projection matrix
  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> dPmat_M(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const int &d){
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bvx(n, n);
    bvx.setZero();
    return(bvx);
  }


}
