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

  T logdetJfromM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &z) override {
    T out;
    out = 0.;
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

  ////////////////////APPROX HELPERS/////////////////////
  bool close2bdry(const Eigen::Matrix<double, Eigen::Dynamic, 1> x, double threshold) {
    stop("Not implemented yet");
    return(true);
  };

  Eigen::Matrix<T, Eigen::Dynamic, 1> approxcentre(const Eigen::Matrix<T, Eigen::Dynamic, 1> x, const double shiftsize=1E-5) {//do nothing by default
    stop("Not approxcentre() not yet implemented for this manifold.");
    Eigen::Matrix<T, Eigen::Dynamic, 1> out(x.size());
    out = x;
    return(out);
  };

  //automatically choose approximation centre
  Eigen::Matrix<T, Eigen::Dynamic, 1> taylorapprox_bdry(
		  CppAD::ADFun<T> &f,
		  const size_t order,
		  const Eigen::Matrix<T, Eigen::Dynamic, 1> xbeta,
		  double shiftsize=1E-5){
     Eigen::Matrix<T, Eigen::Dynamic, 1> x(xbeta.size() / 2);
     x << xbeta.block(0,0,x.size(), 1);
     Eigen::Matrix<T, Eigen::Dynamic, 1> shiftdir(x.size());
     shiftdir.setOnes();
     shiftdir *= -1 * shiftsize;
     shiftdir = Pmat_M(x) * shiftdir;

     Eigen::Matrix<T, Eigen::Dynamic, 1> centre(xbeta.size());
     centre = xbeta;
     centre.block(0,0,x.size(), 1) << x + shiftdir;
     std::cout << "Approximation centre is:" << shiftdir.transpose() << std::endl;
     Eigen::Matrix<T, Eigen::Dynamic, 1> out(0);
     out = taylorapprox(f, centre, order, xbeta);
     return(out);
  }

};
