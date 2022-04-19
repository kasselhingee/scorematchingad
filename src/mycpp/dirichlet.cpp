
namespace { // begin the empty namespace

    a2type ll(const veca1 &beta,
	       const veca2 &u)
    {   size_t n  = beta.size();
        a2type y(0.);  // initialize summation
        for(size_t i = 0; i < n; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }

    // define a function that tapes the above function
      CppAD::ADFun<a1type> tapellS(veca1 zbeta,
                                   a2type (*llf)(const veca1 &, const veca2 &)
                                   ){
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
      y[0] += llf(beta, u);
      y[0] += Spos::logdetJ_fromS(z);

      CppAD::ADFun<a1type> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
      tape.Dependent(z, y);
      return(tape);
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

