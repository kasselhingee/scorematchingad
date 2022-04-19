# include "cdabyppi_types.h"
# include "mycpp/approx.cpp"
# include "mycpp/sm_possphere.cpp"
# include "mycpp/manifold_simplex.cpp"
# include "mycpp/hfuns.cpp"
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
CppAD::ADFun<double> tapesmo_simplex(svecd xbetain, size_t n){
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
    Pmat = simplex::Pmat_M(x);
    a1type h2;
    h2 = prodsq(x);

    // taping ll (log likelihood) store operation sequence
    CppAD::ADFun<a1type> ll;
    ll = tapellsimplex(xbeta);

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
       lapl[0] += Pmat.row(i) * simplex::dPmat_M(x, i) * jac;
    }
    Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> hess(n * n, 1);
    hess = ll.Hessian(x, 0); //the zero here is something about selecting the range-space component of f, 0 selects the first and only component, I presume.
    hess.resize(n, n);
    lapl[0] += (Pmat*hess).trace();
    lapl[0] *= h2; //weight by h2


    //ghPg
    veca1 ghPg(1);
    ghPg = gradprodsq(x).transpose() * simplex::Pmat_M(x) * jac;//jac;

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
XPtr< CppAD::ADFun<double> > ptapesmo_simplex(svecd xbetain, size_t n){
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapesmo_simplex(xbetain, n);
  XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

