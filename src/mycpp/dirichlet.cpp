
namespace { // begin the empty namespace

    // define the log likelihood, with transformation to the sphere, for the Dirichlet distribution
    template <class a1type, class a2type>
    a1type llS(const Eigen::Matrix<a2type, Eigen::Dynamic, 1> &a,
	       const Eigen::Matrix<a1type, Eigen::Dynamic, 1> &x)
    {   size_t n  = a.size();
        a1type y(0.);  // initialize summation
        for(size_t i = 0; i < n; i++)
        {   y   += (1 + 2 * a[i]) * log(x[i]);  
        }
        return y;
    }
   
    // define a function that tapes the above function 
    template <class a1type>
      CppAD::ADFun<a1type> tapellS(Eigen::Matrix<a1type, Eigen::Dynamic, 1> xbeta){
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
      ay[0] = llS(a, ax);     // record operations that compute ay[0]
      CppAD::ADFun<a1type> f;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
      f.Dependent(ax, ay);
      return(f);
  }
    



  // tape grad.grad
  template <class a1type>
  CppAD::ADFun<double> tapegraddotgrad(Eigen::Matrix<a1type, Eigen::Dynamic, 1> xbeta){
    size_t n = xbeta.size() / 2;

    //P matrix instantiate
    Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> Pmat(3, 3); 

    // range space vector
    size_t m = 1;               // number of ranges space variables

    //start outer taping: relationship between exponents and log-likelihood via data vector ax
    CppAD::Independent(xbeta);

    // taping llS store operation sequence in f: X -> Y and stop recording
    CppAD::ADFun<a1type> f;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
    f = tapellS(xbeta);




    // tape derivative calculation, for each beta
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> jac(m * n); // Jacobian of f (m by n matrix)
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> jacsq(1); // Jacobian of f (m by n matrix)
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> x(n);       // domain space vector - must be a1type for plugging into Jacobian function below. Assigned from the outer container xbeta for tracking dependence
    for(int i=0; i < n; i++){
      x[i] = xbeta[i];
    }

    jac  = f.Jacobian(x);      // Jacobian for operation sequence
    Pmat = Pmat_S(x);
    jacsq[0] = (Pmat * jac).dot(jac);

    //finish taping
    CppAD::ADFun<double> g;
    g.Dependent(xbeta, jacsq);
    g.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
    g.check_for_nan(false); //no error if some of the results of the Jacobian are nan.
    return(g);
  }

  //tape laplacian
  template <class a1type>
  CppAD::ADFun<double> tapelaplacian(Eigen::Matrix<a1type, Eigen::Dynamic, 1> xbeta){
    size_t n = xbeta.size() / 2;
    CppAD::Independent(xbeta);
    CppAD::ADFun<a1type> f = tapellS(xbeta); //first get gradient
    Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n); 
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> x(n);
    x = xbeta.block(0,0,n,1);
    //std::cout << xbeta.block(0,0,n,1);
    Pmat = Pmat_S(x);
    Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> dPmat_1(n, n); 
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> lapl(1);
    lapl[0] = 0.;
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> jac(n);
    jac = f.Jacobian(x);
    for(int i=0; i < n; i++){
       lapl[0] += Pmat.row(i) * dPmat_S(x, i) * jac;
    }

    Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> hess(n * n, 1); 
    hess = f.Hessian(x, 0); //the zero here is something about selecting the range-space component of f, 0 selects the first and only component, I presume.
    hess.resize(n, n);
    std::cout << "Hessian is: " << hess << std::endl;
    lapl[0] += (Pmat*hess).trace();
    std::cout << lapl << std::endl;

    //finish taping
    CppAD::ADFun<double> laplf;
    laplf.Dependent(xbeta, lapl);
    laplf.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
    return(laplf);
  }
    
  

  ///////////////////////End Weight Function ////////////

  //tape gradhPgrad
  template <class a1type>
  CppAD::ADFun<double> tapegradhPgrad(Eigen::Matrix<a1type, Eigen::Dynamic, 1> xbeta){
    size_t n = xbeta.size()/2;
    CppAD::Independent(xbeta);
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> x(n);
    x = xbeta.block(0,0,n,1);
    CppAD::ADFun<a1type> f; 
    f = tapellS(xbeta);
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> jac(n);
    jac = f.Jacobian(x);
    Eigen::Matrix<a1type, Eigen::Dynamic, 1> ghPg;
    ghPg = gradprodsq(x).transpose() * Pmat_S(x) * jac;//jac; //gradprodsq(x).transpose().eval() *

    CppAD::ADFun<double> g;
    g.Dependent(xbeta, ghPg);
    g.optimize();
    return(g);
  }
}
