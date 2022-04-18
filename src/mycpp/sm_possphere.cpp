// code for various tools for the positive quadrant of the sphere
namespace { // begin the empty namespace

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> toS(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
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


  //////////////////////////////////////////
  // weight function and grad(h^2) functions
  // squared without constraint
  template <class Type>
  Type prodsq(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
    Type out;
    out = x.array().square().prod();
    return(out);
  }

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> gradprodsq(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
    size_t n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(n);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> avoidone(n-1);
    Type prodx;
    prodx = 2 * x.array().prod();
    for (int i=0; i < n; i++){
	avoidone << x.head(i), x.tail(n-i-1);
	out[i] = prodx * avoidone.prod();
    }
    return(out);
  }

////////////////////APPROX HELPERS/////////////////////
  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> taylorapprox(
		  CppAD::ADFun<Type> &f,
		  const Eigen::Matrix<Type, Eigen::Dynamic, 1> centre,
		  const size_t order,
		  const Eigen::Matrix<Type, Eigen::Dynamic, 1> x){
    //In CppAD speak consider the input function X(t) to be
    //X(t) = centre + t*(x - centre). So X(0) = centre, X(1) = x.
    //First derivative of X at 0, is x - centre
    //Higher order derivative are all 0.
    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(1);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> diff(x.size());
    out.setZero();
    out += f.Forward(0, centre); //zeroth order - constant component of taylor expansion
    if (order >= 1){
      diff = x - centre; // for some reason forward can't take the lazy evaluation version of x - centre direclty. (x - centre).eval() also works
      out += f.Forward(1, diff); //now out[0] is evaluation of a linear approximation of f
    }
    if (order >= 2){
      for (int i=2; i<=order; i++){
	diff.setZero();
        out += f.Forward(i, diff); //now out[0] is evaluation of a quadratic approximation of f
      }
    }
    return(out);
  }

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
