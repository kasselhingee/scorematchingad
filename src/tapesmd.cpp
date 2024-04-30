# include "tapesmd.h"

// function that tapes the score-matching objective
CppAD::ADFun<double> tapesmd(veca1 u, //a vector. The composition measurement for taping
                             veca1 theta, //a vector of parameters for taping
                             CppAD::ADFun<double> & lltape, //tape of ll on the *start* manifold
                             transform<a1type> &tran,
                             manifold<a1type> &M,
                             a1type (*h2fun)(const veca1 &, const double &), // the weight function h^2
                             const double & acut, //the acut constraint for the weight functions
                             bool verbose
                             ){
    size_t d(lltape.Domain()); //manifold is embedded in Euclidean space of dimension d (may be different to simplex)

    //Projection matrix
    mata1 Pmat(d, d);

    //check inputs and ll tape match
    if (lltape.Domain() != (unsigned)u.size()){Rcpp::stop("Dimension of x is %i, which does not match domain size of log likelihood function of %i.", u.size(), lltape.Domain());}
    if (lltape.size_dyn_ind() != theta.size()){Rcpp::stop("Size of dynamic parameter vector (=%i) does not match parameter size of log likelihood tape (=%i).", theta.size(), lltape.size_dyn_ind());}

    //convert to tape from manifold
    CppAD::ADFun<double> lltapeman;
    lltapeman = tapellman(u, //a vector. The measurement for taping in the interior of the domain of lltape
                           theta, //a vector of the dynamic parameters for taping
                           lltape,
                           tran,
                           verbose);

    // make ll tape higher order (i.e. as if it was taped using a2type instead of a1type)
    CppAD::ADFun<a1type, double> lltapehigher;
    lltapehigher = lltapeman.base2ad();

    //get h2 tape
    CppAD::ADFun<double> dh2tape;
    veca1 z(d); //size of z changed by toM result below
    z = tran.toM(u); //transform u to the end manifold
    CppAD::ADFun<a1type, double> h2tape; //The second type here 'double' is for the 'RecBase' in ad_fun.hpp. It doesn't seem to change the treatment of the object.
    h2tape = tapeh2(z, h2fun, acut).base2ad(); //convert to a function of a1type rather than double

    //START TAPING
    CppAD::Independent(theta, u); //differentiate wrt beta, dynamic parameter is u
    if (verbose){
      PrintForVec("\n\n(START) theta is: ", theta);
      PrintForVec("\n\n(START) u is: ", u);
    }

    // veca1 u(n);
    z = tran.toM(u); //transform u to the manifold

    Pmat = M.Pmatfun(z);
    if (verbose){PrintForMatrix("\nThe value of Pmat is: \n", Pmat);}
    veca1 h2(1);
    veca1 gradh2(d); //will be set by h2tape.Jacobian below anyway
    h2 = h2tape.Forward(0, z);
    gradh2 = h2tape.Jacobian(z);
    if (verbose){
      PrintForVec("\n h2 is: ", h2);
      PrintForVec("\n gradh2 is: ", gradh2);
    }

    //update parameters ('dynamic' values) of lltape
    lltapehigher.new_dynamic(theta);

    //grad(ll)
    veca1 jac(d); // Jacobian of ll
    jac  = lltapehigher.Jacobian(z);      // Jacobian for operation sequence
    if (verbose){PrintForVec("\nThe value of ll gradient is: ", jac);}

    //hgPg
    veca1 hgPg(1);
    hgPg[0] = 0.5 * h2[0] * (Pmat * jac).dot(jac);
    if (verbose){PrintForVec("\n(1)The value of hgPg is:", hgPg);}

    //hlap
    veca1 lapl(1);
    lapl[0] = 0.;
    for(size_t i=0; i < d; i++){
       lapl[0] += Pmat.row(i) * M.dPmatfun(z, i) * jac;
       if (verbose){PrintForMatrix("\nThe value of dPmat is: \n", M.dPmatfun(z, i));}
    }
    mata1 hess(d * d, 1);
    hess = lltapehigher.Hessian(z, 0); //the zero here is something about selecting the range-space component of f, 0 selects the first and only component, I presume.
    hess.resize(d, d);
    if (verbose){PrintForMatrix("\nThe value of ll Hessian is: \n", hess);}
    lapl[0] += (Pmat*hess).trace();
    lapl[0] *= h2[0]; //weight by h2
    if (verbose){PrintForVec("\n(2)The value of hlapl is:", lapl);}


    //ghPg
    veca1 ghPg(1);
    ghPg = gradh2.transpose() * Pmat * jac;//jac; //gradprodsq(x).transpose().eval() *
    if (verbose){PrintForVec("\n(3)The value of ghPg is:", ghPg);}

    //combine components
    veca1 smd(1);
    smd = lapl + hgPg + ghPg;
    if (verbose){PrintForVec("\n(END)The value of smd is:", smd);}

    //finish taping
    CppAD::ADFun<double> smdfun;
    smdfun.Dependent(theta, smd);
    //smdfun.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore. But asserts errors.
    smdfun.check_for_nan(false); //no error if some of the results of the Jacobian are nan.
    return(smdfun);
}


Rcpp::XPtr< CppAD::ADFun<double> > ptapesmd(veca1 u_ad,
                                      veca1 theta_ad,
                                      Rcpp::XPtr< CppAD::ADFun<double> > pll,
                                      transform_a1type & tran,
                                      manifold_a1type & man,
                                      std::string weightname,
                                      const double acut,
                                      bool verbose){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer

  //choose weight function
  a1type (*h2fun)(const veca1 &, const double &) = nullptr;
  if (weightname.compare("prodsq") == 0){
    h2fun = bdryweight::prodsq;
  }
  if (weightname.compare("minsq") == 0){
    h2fun = bdryweight::minsq;
  }
  if (weightname.compare("ones") == 0){
    h2fun = bdryweight::oneweights;
  }

  //check weight function
  if (h2fun == nullptr){
    throw std::invalid_argument("Matching weight function not found");
  }

  *out = tapesmd(u_ad,
                 theta_ad,
                 *pll,
                 tran,
                 man,
                 h2fun,
                 acut,
                 verbose);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

// function that combines a ll tape with a manifold transformation
CppAD::ADFun<double> tapellman(veca1 x, //a vector. The measurement for taping in the interior of the domain of lltape
                             veca1 thetavar, //a vector of the dynamic parameters for taping
                             CppAD::ADFun<double> & lltape,
                             transform<a1type> &tran,
                             bool verbose
                             ){
  // convert taping specified by x to the new manifold
  veca1 z(0); //0 here because size dictated by toM
  z = tran.toM(x);
  if (verbose){
    Rcpp::Rcout << "z=toM(x) is: " << z.transpose() << std::endl;
  }

  // make ll tape higher order (i.e. as if it was taped using a2type instead of a1type)
  CppAD::ADFun<a1type, double> lltapehigher;
  lltapehigher = lltape.base2ad();
  
  //tape relationship between z and log-likelihood
  CppAD::Independent(z, thetavar);  //for this tape, theta must be altered using new_dynamic
  if (verbose){
    Rcpp::Rcout << "x on end manifold is: " << z.transpose() << std::endl;
    PrintForVec("\n x on end manifold is: ", z);
    Rcpp::Rcout << "dynamic theta elements are: " << thetavar.transpose() << std::endl;
    PrintForVec("\n dynamic theta elements are: ", thetavar);
  }

  veca1 u(0); //0 here because size dictated by fromM
  u = tran.fromM(z);
  if (verbose){
    Rcpp::Rcout << "x on start manifold is: " << u.transpose() << std::endl;
    PrintForVec("\n x on start manifold is: ", u);
  }
  
  
  //Now evaluate
  //update parameters ('dynamic' values) of lltape
  lltapehigher.new_dynamic(thetavar);
  veca1 y(1); // vector of ranges space variables
  y = lltapehigher.Forward(0, u);

  //get log determinant of fromM
  y[0] += tran.logdetJfromM(z);

  //finish taping
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  if (verbose){
    Rcpp::Rcout << "tape has " << tape.size_dyn_ind() << " independent dynamic parameters" << std::endl;
    Rcpp::Rcout << "tape requires vectors of length " << tape.Domain() << std::endl;
    Rcpp::Rcout << "tape returns vectors of length " << tape.Range() << std::endl;
  }
  return(tape);
}

Rcpp::XPtr< CppAD::ADFun<double> > ptapellman(veca1 x, //a vector. The measurement for taping in the interior of the domain of lltape
                             veca1 thetavar, //a vector of the dynamic parameters for taping
                             Rcpp::XPtr< CppAD::ADFun<double> > plltape,
                             transform<a1type> &tran,
                             bool verbose
                             ){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapellman(x,
                 thetavar,
                 *plltape,
                 tran,
                 verbose);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);

}
