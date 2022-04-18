
    template <class a1type, class a2type>
    a1type llsimplex(const Eigen::Matrix<a2type, Eigen::Dynamic, 1> &a,
	       const Eigen::Matrix<a1type, Eigen::Dynamic, 1> &z)
    {   size_t n  = a.size();
	      Eigen::Matrix<a1type, Eigen::Dynamic, 1> u(z.size());
	      u = simplex::fromM(z);
        a1type y(0.);  // initialize summation
        y += ll(a, u);
        y += simplex::logdetJ_fromM(z);
        return y;
    }

    // define a function that tapes the above function
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
      ay[0] = llsimplex(a, ax);     // record operations that compute ay[0]
      CppAD::ADFun<a1type> f;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
      f.Dependent(ax, ay);
      return(f);
  }


