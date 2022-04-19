
namespace { // begin the empty namespace

    template <class a1type, class a2type>
    a1type ll(const Eigen::Matrix<a2type, Eigen::Dynamic, 1> &a,
	       const Eigen::Matrix<a1type, Eigen::Dynamic, 1> &u)
    {   size_t n  = a.size();
        a1type y(0.);  // initialize summation
        for(size_t i = 0; i < n; i++)
        {   y   += a[i] * log(u[i]);
        }
        return y;
    }

    template <class a1type, class a2type>
    a1type llS(const Eigen::Matrix<a2type, Eigen::Dynamic, 1> &a,
	       const Eigen::Matrix<a1type, Eigen::Dynamic, 1> &z)
    {   size_t n  = a.size();
	      Eigen::Matrix<a1type, Eigen::Dynamic, 1> u(z.size());
	      u = Spos::fromS(z);
        a1type y(0.);  // initialize summation
        y += ll(a, u);
        y += Spos::logdetJ_fromS(z);
        return y;
    }

    // define a function that tapes the above function
    template <class a1type>
      CppAD::ADFun<a1type> tapellS(Eigen::Matrix<a1type, Eigen::Dynamic, 1> zbeta){
      typedef CppAD::AD<a1type> a2type;  // for second (inner) level of taping
      size_t n = 3;                  // number of dimensions

      Eigen::Matrix<a1type, Eigen::Dynamic, 1> beta(n); // vector of exponents in the outer type
      //declare dummy internal level of taping variables:
      Eigen::Matrix<a2type, Eigen::Dynamic, 1> z(n); // vector of domain space variables
      for(int i = 0; i < n; i++){
         beta[i] = zbeta[i + n];
         z[i] = zbeta[i];
      }

      //tape relationship between x and log-likelihood
      CppAD::Independent(z);

      // range space vector
      size_t m = 1;               // number of ranges space variables
      Eigen::Matrix<a2type, Eigen::Dynamic, 1> y(m); // vector of ranges space variables
      Eigen::Matrix<a2type, Eigen::Dynamic, 1> u(z.size());
      u = Spos::fromS(z);
      y.setZero();
      y[0] += ll(beta, u);
      y[0] += Spos::logdetJ_fromS(z);

      CppAD::ADFun<a1type> f;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
      f.Dependent(z, y);
      return(f);
  }

    template <class a1type>
    CppAD::ADFun<a1type> tapellsimplex(Eigen::Matrix<a1type, Eigen::Dynamic, 1> xbeta){
      typedef CppAD::AD<a1type> a2type;  // for second (inner) level of taping
      size_t n = 3;                  // number of dimensions

      Eigen::Matrix<a1type, Eigen::Dynamic, 1> a(n); // vector of exponents in the outer type
      //declare dummy internal level of taping variables:
      Eigen::Matrix<a2type, Eigen::Dynamic, 1> ax(n); // vector of domain space variables
      for(int i = 0; i < n; i++){
        a[i] = xbeta[i + n];
        ax[i] = 3.;
      }

      //tape relationship between x and log-likelihood
      CppAD::Independent(ax);

      // range space vector
      size_t m = 1;               // number of ranges space variables
      Eigen::Matrix<a2type, Eigen::Dynamic, 1> ay(m); // vector of ranges space variables
      ay.setZero();
      ay[0] += ll(a, ax);

      CppAD::ADFun<a1type> f;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
      f.Dependent(ax, ay);
      return(f);
    }


}

