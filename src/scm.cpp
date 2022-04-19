# include "cdabyppi_types.h"
# include "mycpp/approx.cpp"
# include "mycpp/sm_possphere.cpp"
# include "mycpp/manifold_simplex.cpp"
# include "mycpp/hfuns.cpp"
# include "mycpp/dirichlet.cpp"
using namespace Rcpp;

// define a function that tapes a log likelihood
CppAD::ADFun<a1type> tapell(veca1 zbeta,
                               a2type (*llf)(const veca1 &, const veca2 &), //the log likelihood function
                               veca2 (*fromM)(const veca2 &), //transformation from manifold to simplex
                               a2type (*logdetJfromM)(const veca2 &) //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
                               ){
  size_t n = 3;                  // number of dimensions
  veca1 beta(n); // vector of exponents in the outer type
  //declare dummy internal level of taping variables:
  veca2 z(n); // vector of domain space variables
  for(int i = 0; i < n; i++){
     beta[i] = zbeta[i + n];
     z[i] = zbeta[i];
  }

  //tape relationship between x and log-likelihood
  CppAD::Independent(z);
  // range space vector
  size_t m = 1;               // number of ranges space variables
  veca2 y(m); // vector of ranges space variables
  veca2 u(z.size());
  u = fromM(z);
  y.setZero();
  y[0] += llf(beta, u);
  y[0] += logdetJfromM(z);
  CppAD::ADFun<a1type> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  return(tape);
}

// function that tapes the score-matching objective
CppAD::ADFun<double> tapesmo(svecd ubetain, //a vector. The first n elements is the measurement, the remaining elements are the parameters
                             size_t n, //the dimension of the measurements (number of components)
                             a2type (*llf)(const veca1 &, const veca2 &), //the log likelihood function
                             veca1 (*toM)(const veca1 &), //map from simplex to manifold
                             mata1 (*Pmatfun)(const veca1 &), //projection matrix for manifold
                             mata1 (*dPmatfun)(const veca1 &, const int &),//elementwise derivative of projection matrix for manifold
                             veca2 (*fromM)(const veca2 &), //transformation from manifold to simplex
                             a2type (*logdetJfromM)(const veca2 &), //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
                             a1type (*h2fun)(const veca1 &), // the weight function h^2
                             veca1 (*gradh2fun)(const veca1 &)// the gradient of the weight function h^2
                             ){
    veca1 ubeta(ubetain.size());
    for (int i=0; i < ubetain.size(); i++){
       ubeta[i] = ubetain[i];
    }

    //Projection matrix
    mata1 Pmat(n, n);

    //START TAPING
    CppAD::Independent(ubeta);

    veca1 u(n);
    u = ubeta.block(0,0,n,1);
    veca1 z(n);
    z = toM(u); //transform u to the manifold
    veca1 zbeta(ubeta.size());
    for (size_t i=0; i<n; i++){
      zbeta[i] = z[i];
      zbeta[n + i] = ubeta[n+i];
    }

    Pmat = Pmatfun(z);
    a1type h2;
    h2 = h2fun(z);

    // taping ll (log likelihood) store operation sequence
    CppAD::ADFun<a1type> lltape;
    lltape = tapell(zbeta, llf, fromM, logdetJfromM);

    //grad(ll)
    veca1 jac(n); // Jacobian of ll
    jac  = lltape.Jacobian(z);      // Jacobian for operation sequence

    //hgPg
    veca1 hgPg(1);
    hgPg[0] = 0.5 * h2 * (Pmat * jac).dot(jac);

    //hlap
    veca1 lapl(1);
    lapl[0] = 0.;
    for(size_t i=0; i < n; i++){
       lapl[0] += Pmat.row(i) * dPmatfun(z, i) * jac;
    }
    mata1 hess(n * n, 1);
    hess = lltape.Hessian(z, 0); //the zero here is something about selecting the range-space component of f, 0 selects the first and only component, I presume.
    hess.resize(n, n);
    lapl[0] += (Pmat*hess).trace();
    lapl[0] *= h2; //weight by h2


    //ghPg
    veca1 ghPg(1);
    ghPg = gradh2fun(z).transpose() * Pmat * jac;//jac; //gradprodsq(x).transpose().eval() *

    //combine components
    veca1 smo(1);
    smo = lapl + hgPg + ghPg;

    //finish taping
    CppAD::ADFun<double> smofun;
    smofun.Dependent(ubeta, smo);
    smofun.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
    smofun.check_for_nan(false); //no error if some of the results of the Jacobian are nan.
    return(smofun);
}

//in R store a pointer to the ADFun object
//' @title The score matching objective calculator.
//' @param xbetain a concatenated vector of sqrt(x) and beta
//' @param n The dimension of x.
//' @return An RCpp::XPtr object pointing to the ADFun
//' @export
// [[Rcpp::export]]
XPtr< CppAD::ADFun<double> > ptapesmo(svecd xbetain,
                                      size_t n,
                                      std::string manifoldname,
                                      std::string weightname){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer

  //choose weight function
  a1type (*h2fun)(const veca1 &) = nullptr;
  veca1 (*gradh2fun)(const veca1 &) = nullptr;// the gradient of the weight function h^2
  if (weightname.compare("prodsq") == 0){
    h2fun = prodsq;
    gradh2fun = gradprodsq;
  }
  if (weightname.compare("prod1") == 0){
    h2fun = hprod;
    gradh2fun = gradhprod;
  }
  //check weight function
  if (h2fun == nullptr){
    throw std::invalid_argument("Matching weight function not found");
  }
  if (gradh2fun == nullptr){
    throw std::invalid_argument("Matching weight function gradient not found");
  }

  if (manifoldname.compare("sphere") == 0){
  *out = tapesmo(xbetain, n, ll,
                 Spos::toS, Spos::Pmat_S, Spos::dPmat_S,
                 Spos::fromS, Spos::logdetJ_fromS,
                 h2fun, gradh2fun);
  }
  if (manifoldname.compare("simplex") == 0){
  *out = tapesmo(xbetain, n, ll,
                 simplex::toM, simplex::Pmat_M, simplex::dPmat_M,
                 simplex::fromM, simplex::logdetJ_fromM,
                 h2fun, gradh2fun);
  }

  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

//calc smo and additions
//use a pointer to an ADFun object to compute the function evaluated  at a location
//' @title The score matching objective calculator.
//' @param u A vector in the simplex.
//' @param betain
//' @return The score matching objective value
//' @export
// [[Rcpp::export]]
double psmo(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain){
  //convert input to an Eigen vectors
  vecd u_e(u.size());
  for (int i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  vecd beta_e(betain.size());
  for (int i=0; i<betain.size(); i++){
    beta_e[i] = betain[i];
  }

  vecd xbetain(u_e.size() + beta_e.size());
  xbetain << u_e, beta_e;
  vecd smo_val(1);
  smo_val = pfun->Forward(0, xbetain);  //treat the XPtr as a regular pointer
  return(smo_val[0]);
}

//calc smo and additions
//use a pointer to an ADFun object to compute the function evaluated  at a location
//' @title The score matching objective calculator.
//' @param u A vector in the simplex.
//' @param betain
//' @return The score matching objective value
//' @export
// [[Rcpp::export]]
svecd psmograd(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain){
  //convert input to an Eigen vectors
  vecd u_e(u.size());
  for (int i=0; i<u.size(); i++){
    u_e[i] = u[i];
  }
  vecd beta_e(betain.size());
  for (int i=0; i<betain.size(); i++){
    beta_e[i] = betain[i];
  }


  vecd xbetain(u_e.size() + beta_e.size());
  xbetain << u_e, beta_e;
  vecd sc_grad(xbetain.size());
  vecd out_e(beta_e.size());
  svecd out(beta_e.size());
  sc_grad = pfun->Jacobian(xbetain);  //treat the XPtr as a regular pointer
  out_e = sc_grad.block(u_e.size(),0,beta_e.size(),1);

  //convert to std::vector
  for (int i = 0; i<u_e.size(); i++){
    out[i] = out_e[i];
  }
  return(out);
}

