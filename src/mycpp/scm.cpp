# include "../cdabyppi_types.h"
# include "approx.cpp"
# include "sm_possphere.cpp"
# include "manifold_simplex.cpp"
# include "hfuns.cpp"
# include "dirichlet.cpp"

// define a function that tapes a log likelihood
CppAD::ADFun<double> tapell(veca1 z, //data measurement tranformed to M manifold
                            veca1 theta, //theta parameter
                               a1type (*llf)(const veca1 &, const veca1 &), //the log likelihood function
                               veca1 (*fromM)(const veca1 &), //transformation from manifold to simplex
                               a1type (*logdetJfromM)(const veca1 &) //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
                               ){

  //tape relationship between x and log-likelihood
  CppAD::Independent(z, theta);  //for this tape, beta must be altered using new_dynamic
  // range space vector
  veca1 y(1); // vector of ranges space variables
  veca1 u(z.size());
  u = fromM(z);
  y.setZero();
  y[0] += llf(theta, u);
  y[0] += logdetJfromM(z);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  return(tape);
}

CppAD::ADFun<double> tapeh2(veca1 z,
                            a1type (*h2fun)(const veca1 &, const double &),
                            const double & acut){
  //tape relationship between z and h2
  CppAD::Independent(z);
  // range space vector
  size_t m = 1;               // number of ranges space variables
  veca1 y(m); // vector of ranges space variables
  y[0] = h2fun(z, acut);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  return(tape);
}

// function that tapes the score-matching objective
CppAD::ADFun<double> tapesmo(svecd ubetain, //a vector. The first n elements is the measurement, the remaining elements are the parameters
                             size_t n, //the dimension of the measurements (number of components)
                             a1type (*llf)(const veca1 &, const veca1 &), //the log likelihood function
                             manifold<a1type> &M,
                             a1type (*h2fun)(const veca1 &, const double &), // the weight function h^2
                             const double & acut //the acut constraint for the weight functions
                             ){
    veca1 ubeta(ubetain.size());
    for (size_t i=0; i < ubetain.size(); i++){
       ubeta[i] = ubetain[i];
    }
    size_t betasize = ubetain.size() - n;

    //Projection matrix
    mata1 Pmat(n, n);

    //get h2 tape
    CppAD::ADFun<double> dh2tape;
    veca1 u(n);
    u = ubeta.block(0,0,n,1);
    veca1 z(n);
    z = M.toM(u); //transform u to the manifold
    CppAD::ADFun<a1type, double> h2tape; //The second type here 'double' is for the 'RecBase' in ad_fun.hpp. It doesn't seem to change the treatment of the object.
    h2tape = tapeh2(z, h2fun, acut).base2ad(); //convert to a function of a1type rather than double

    // taping ll (log likelihood) store operation sequence
    veca1 beta(ubeta.size() - n);
    beta = ubeta.block(n, 0, beta.size(), 1);
    CppAD::ADFun<a1type, double> lltape;
    lltape = tapell(z, beta, llf, M.fromM, M.logdetJfromM).base2ad();

    //START TAPING
    CppAD::Independent(beta, u); //differentiate wrt beta, dynamic parameter is u

    // veca1 u(n);
    z = M.toM(u); //transform u to the manifold

    Pmat = M.Pmatfun(z);
    veca1 h2(1);
    veca1 gradh2(n);
    h2 = h2tape.Forward(0, z);
    gradh2 = h2tape.Jacobian(z);
    // CppAD::PrintFor("For h2fun, z is: ", z[0]);
    // for(size_t i=1; i<z.size(); i++){
      // CppAD::PrintFor(" ", z[i]);
    // }
    // CppAD::PrintFor(" The value of h2 is: ", h2[0]);
    // CppAD::PrintFor("The value of gradh2 is: ", gradh2[0]);
    // for(size_t i=1; i<gradh2.size(); i++){
      // CppAD::PrintFor(" ", gradh2[i]);
    // }

    //update parameters ('dynamic' values) of lltape
    lltape.new_dynamic(beta);

    //grad(ll)
    veca1 jac(n); // Jacobian of ll
    jac  = lltape.Jacobian(z);      // Jacobian for operation sequence

    //hgPg
    veca1 hgPg(1);
    hgPg[0] = 0.5 * h2[0] * (Pmat * jac).dot(jac);

    //hlap
    veca1 lapl(1);
    lapl[0] = 0.;
    for(size_t i=0; i < n; i++){
       lapl[0] += Pmat.row(i) * M.dPmatfun(z, i) * jac;
    }
    mata1 hess(n * n, 1);
    hess = lltape.Hessian(z, 0); //the zero here is something about selecting the range-space component of f, 0 selects the first and only component, I presume.
    hess.resize(n, n);
    lapl[0] += (Pmat*hess).trace();
    lapl[0] *= h2[0]; //weight by h2


    //ghPg
    veca1 ghPg(1);
    ghPg = gradh2.transpose() * Pmat * jac;//jac; //gradprodsq(x).transpose().eval() *

    //combine components
    veca1 smo(1);
    smo = lapl + hgPg + ghPg;

    //finish taping
    CppAD::ADFun<double> smofun;
    smofun.Dependent(beta, smo);
    smofun.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
    smofun.check_for_nan(false); //no error if some of the results of the Jacobian are nan.
    return(smofun);
}
