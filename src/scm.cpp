# include "cdabyppi_types.h"
# include "mycpp/approx.cpp"
# include "mycpp/sm_possphere.cpp"
# include "mycpp/hfuns.cpp"
# include "mycpp/dirichlet.cpp"
using namespace Rcpp;

// smo functions
CppAD::ADFun<double> tapesmo(svecd ubetain, size_t n,
                             veca1 (*toM)(const veca1 &),
                             mat1 (*Pmatfun)(const veca1 &)
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
    for (int i=0; i<n; i++){
      zbeta[i] = z[i];
      zbeta[n + i] = ubeta[n+i];
    }

    Pmat = Pmatfun(z);
    a1type h2;
    h2 = prodsq(z);

    // taping ll (log likelihood) store operation sequence
    CppAD::ADFun<a1type> lltape;
    lltape = tapell(zbeta, ll, Spos::fromS, Spos::logdetJ_fromS);

    //grad(ll)
    veca1 jac(n); // Jacobian of ll
    jac  = lltape.Jacobian(z);      // Jacobian for operation sequence

    //hgPg
    veca1 hgPg(1);
    hgPg[0] = 0.5 * h2 * (Pmat * jac).dot(jac);

    //hlap
    veca1 lapl(1);
    lapl[0] = 0.;
    for(int i=0; i < n; i++){
       lapl[0] += Pmat.row(i) * Spos::dPmat_S(z, i) * jac;
    }
    mata1 hess(n * n, 1);
    hess = lltape.Hessian(z, 0); //the zero here is something about selecting the range-space component of f, 0 selects the first and only component, I presume.
    hess.resize(n, n);
    lapl[0] += (Pmat*hess).trace();
    lapl[0] *= h2; //weight by h2


    //ghPg
    veca1 ghPg(1);
    ghPg = gradprodsq(z).transpose() * Pmat * jac;//jac; //gradprodsq(x).transpose().eval() *

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
XPtr< CppAD::ADFun<double> > ptapesmo(svecd xbetain, size_t n){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapesmo(xbetain, n, Spos::toS, Spos::Pmat_S);
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

