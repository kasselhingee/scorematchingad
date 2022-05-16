// code for various tools for the positive quadrant of the sphere
template <typename Type>
struct Spos : public manifold<Type> {
  ~Spos(){};
  Spos(){};

  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     // for (int i=0; i<x.size(); i++){
     //   out[i] = CppAD::sqrt(x[i]);
     // }
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


////////////////////APPROX HELPERS/////////////////////
  // function that returns whether point is close to boundary or not, for this particular manifold
  bool close2bdry(const Eigen::Matrix<double, Eigen::Dynamic, 1> x, double threshold){
    double minval;
    minval = x.array().abs().minCoeff(); //this is not quite the distance on the manifold, but close enough
    bool out;
    if (minval < threshold){
      out = true;
    } else {
      out = false;
    }
    return(out);
  }

  //function that produces a new location given a location too close to boundary
  Eigen::Matrix<Type, Eigen::Dynamic, 1> approxcentre(const Eigen::Matrix<Type, Eigen::Dynamic, 1> x,
                                                      const double shiftsize=1E-5) {
    Eigen::Matrix<Type, Eigen::Dynamic, 1> shiftdir(x.size());
    shiftdir.setOnes();
    shiftdir *= -1 * shiftsize;
    shiftdir = Pmat_S(x) * shiftdir; //shift in tangent to manifold at x

    Eigen::Matrix<Type, Eigen::Dynamic, 1> centre(x.size());
    centre = x + shiftdir;
    centre = centre / (centre.array().square().sum().sqrt()); //project onto the manifold by normalising
    return(centre);
  }

  //automatically choose approximation centre
  Eigen::Matrix<Type, Eigen::Dynamic, 1> taylorapprox_bdry(
		  CppAD::ADFun<Type> &f,
		  //above is the smo function - whos arguments are the model parameters, but need to differente wrt the values of a measurement
		  const size_t order,
		  const Eigen::Matrix<Type, Eigen::Dynamic, 1> xbeta,
		  double shiftsize=1E-5) {
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

};
