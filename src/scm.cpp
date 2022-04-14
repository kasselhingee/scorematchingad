# include <iostream>        // standard input/output
# include <vector>          // standard vector
# include <cppad/example/cppad_eigen.hpp>  //load eigen
# include <cppad/cppad.hpp> // the CppAD package
# include "mycpp/sm_possphere.cpp"
# include "mycpp/dirichlet.cpp"
# include <Rcpp.h>
using namespace Rcpp;

typedef std::vector<double> svecd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecd; //a vector of a1type values
typedef CppAD::AD<double> a1type;   // for first (outer) level of taping
typedef Eigen::Matrix<a1type, Eigen::Dynamic, 1> veca1; //a vector of a1type values
typedef CppAD::AD<a1type> a2type;  // for second (inner) level of taping
typedef Eigen::Matrix<a2type, Eigen::Dynamic, 1> veca2;

// smo functions
CppAD::ADFun<double> tapesmo(svecd xbetain, size_t n){
    veca1 xbeta(xbetain.size());
    for (int i=0; i < xbetain.size(); i++){
       xbeta[i] = xbetain[i];
    }

    //Projection matrix
    Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);

    //START TAPING
    CppAD::Independent(xbeta);

    veca1 x(n);
    x = xbeta.block(0,0,n,1);
    Pmat = Pmat_S(x);
    a1type h2;
    h2 = prodsq(x);

    // taping ll (log likelihood) store operation sequence
    CppAD::ADFun<a1type> ll;
    ll = tapellS(xbeta);

    //grad(ll)
    veca1 jac(n); // Jacobian of ll
    jac  = ll.Jacobian(x);      // Jacobian for operation sequence

    //hgPg
    veca1 hgPg(1);
    hgPg[0] = 0.5 * h2 * (Pmat * jac).dot(jac);

    //hlap
    veca1 lapl(1);
    lapl[0] = 0.;
    for(int i=0; i < n; i++){
       lapl[0] += Pmat.row(i) * dPmat_S(x, i) * jac;
    }
    Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> hess(n * n, 1);
    hess = ll.Hessian(x, 0); //the zero here is something about selecting the range-space component of f, 0 selects the first and only component, I presume.
    hess.resize(n, n);
    lapl[0] += (Pmat*hess).trace();
    lapl[0] *= h2; //weight by h2


    //ghPg
    veca1 ghPg(1);
    ghPg = gradprodsq(x).transpose() * Pmat_S(x) * jac;//jac; //gradprodsq(x).transpose().eval() *

    //combine components
    veca1 smo(1);
    smo = lapl + hgPg + ghPg;

    //finish taping
    CppAD::ADFun<double> smofun;
    smofun.Dependent(xbeta, smo);
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
  CppAD::ADFun<double> outobj;
  outobj = tapesmo(xbetain, n);
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapesmo(xbetain, n);
  // out = &outobj;
  XPtr< CppAD::ADFun<double> > pout(out, true);
  std::cout << "XPtr Domain: " << pout->Domain() << std::endl;
  return(pout);
}

//calc smo and additions
//use a pointer to an ADFun object to compute the function evaluated  at a location
//' @title The score matching objective calculator.
//' @param xbetain a concatenated vector of sqrt(x) and beta
//' @param n The dimension of x.
//' @return An RCpp::XPtr object pointing to the ADFun
//' @export
// [[Rcpp::export]]
double psmo_n_grad(XPtr< CppAD::ADFun<double> > pfun, svecd xin, svecd betain){
  vecd xbetain(xin.size() + betain.size());
  for (int i=0; i<xin.size(); i++){
    xbetain[i] = xin[i];
  }
  for (int i=0; i<betain.size(); i++){
    xbetain[i + xin.size()] = betain[i];
  }
  vecd smo_val(1);
  std::cout << "Attempting to copy pfun" << std::endl;
  // CppAD::ADFun<double> f;
  // f = *pfun;
  std::cout << "Attempting to use pfun" << std::endl;
  std::cout << pfun->Domain() << std::endl;
  smo_val = pfun->Forward(0, xbetain);  //treat the XPtr as a regular pointer
  return(smo_val[0]);
}



//' @title The value of the score matching objective.
//'
//' @param xin the composition after sqrt transform
//' @param betain the beta values
//' @return the score matching objective for `xin`
//' @export
// [[Rcpp::export]]
double smo_n_grad(svecd xin, svecd betain){
    // vector of exponents and values for taping
    size_t n = 3;                  // number of dimensions
    svecd xbeta(2*n); //outer container of arguments. First half the x, second half the beta
    for(int i = 0; i < n; i++){
        xbeta[i] = 3.;                 // value at which function is recorded
        xbeta[i + n] = 1.;                 // value of exponents
    }

    //tape of smo
    CppAD::ADFun<double> smofun;
    smofun = tapesmo(xbeta, n);
    std::cout << smofun.Domain() << std::endl;

    /////////////////////Calculate SMO///////////
    //compute score matching objective
    vecd xbetain(xin.size() + betain.size());
    for (int i=0; i<xin.size(); i++){
      xbetain[i] = xin[i];
    }
    for (int i=0; i<betain.size(); i++){
      xbetain[i + xin.size()] = betain[i];
    }
    vecd smo_val(1);
    smo_val = smofun.Forward(0, xbetain);

    //if xin on boundary and smo_val is nan, interpolate from nearby
    if (std::isnan(smo_val[0])){
      std::cout << "SMO is NAN" << std::endl;
      if (xbetain.block(0,0,n,1).minCoeff() < 1E-5){ //can't used xin because it is svecd, not part of Eigen
        smo_val = taylorapprox_bdry(smofun, 100, xbetain, 1E-4);
      }
    }

    vecd sc_grad(2 * n);
    sc_grad = smofun.Jacobian(xbetain);
    std::cout << "Gradient: " << std::endl;
    std::cout << sc_grad.block(n,0,n,1).transpose() << std::endl;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sc_hess(2 * n * 2 * n, 1);
    sc_hess += smofun.Hessian(xbetain, 0);
    sc_hess.resize(2 * n, 2 * n);
    std::cout << "Hessian: " << std::endl;
    std::cout << sc_hess.block(n,n,n,n) << std::endl;

    std::cout << "Score objective is:" << std::endl;
    std::cout << smo_val << std::endl;
    return(smo_val[0]);
}


// main program
int main(int argc, char** argv)
{   using CppAD::AD;   // use AD as abbreviation for CppAD::AD
    //read inputs
    size_t n = 3;
    svecd xin(n);
    svecd betain(n);
    for(int i=0; i < n; i++){
       xin[i] = 3.;
       betain[i] = 1.;
    }
    if (argc > 1){
       for(int i=0; i < n; i++){
         xin[i] = strtod(argv[i + 1], NULL);
         betain[i] = strtod(argv[i + n + 1], NULL);
       }
    }
    // std::cout << "x in is " << xin.transpose() << std::endl;
    // std::cout << "beta in is " << betain.transpose() << std::endl;
    // std::cout << "h2 is " << prodsq(xin) << std::endl;
    // std::cout << "grad(h2) is " << gradprodsq(xin).transpose() << std::endl;

    double smo_val;
    smo_val = smo_n_grad(xin, betain);
    //compute score matching objective

    return 0;
}


