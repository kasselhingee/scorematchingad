# include "cdabyppi_types.h"
# include "mycpp/approx.cpp"
# include "mycpp/sm_possphere.cpp"
# include "mycpp/manifold_simplex.cpp"
# include "mycpp/hfuns.cpp"
# include "mycpp/dirichlet.cpp"

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
  for(size_t i = 0; i < n; i++){
     beta[i] = zbeta[i + n];
     z[i] = zbeta[i];
  }

  //tape relationship between x and log-likelihood
  CppAD::Independent(z);
  std::cout << "Taping ll using the following values: ";
  std::cout << z.transpose() << " " << beta.transpose() << std::endl;
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

CppAD::ADFun<a1type> tapeh2(veca1 zin,
                            a2type (*h2fun)(const veca2 &, const double &),
                            const double & acut){
  veca2 z(zin.size());
  for (size_t i = 0; i < zin.size(); i++){
    z[i] = zin[i];
  }
  //tape relationship between z and h2
  CppAD::Independent(z);
  std::cout << "Taping h2 using the following values: ";
  std::cout << z.transpose() << std::endl;
  // range space vector
  size_t m = 1;               // number of ranges space variables
  veca2 y(m); // vector of ranges space variables
  y[0] = h2fun(z, acut);
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
                             a2type (*h2fun)(const veca2 &, const double &), // the weight function h^2
                             veca1 (*gradh2fun)(const veca1 &, const double &),// the gradient of the weight function h^2
                             const double & acut //the acut constraint for the weight functions
                             ){
    veca1 ubeta(ubetain.size());
    for (size_t i=0; i < ubetain.size(); i++){
       ubeta[i] = ubetain[i];
    }

    //Projection matrix
    mata1 Pmat(n, n);

    //START TAPING
    CppAD::Independent(ubeta);
    std::cout << "Taping smo using the following values: " << ubeta.transpose() << std::endl;

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
    veca1 h2(1);
    veca1 gradh2(z.size());
    CppAD::ADFun<a1type> h2tape;
    h2tape = tapeh2(z, h2fun, acut);
    h2 = h2tape.Forward(0, z);
    gradh2 = h2tape.Jacobian(z);
    CppAD::PrintFor("For h2fun, z is: ", z[0]);
    for(size_t i=1; i<z.size(); i++){
      CppAD::PrintFor(" ", z[i]);
    }
    CppAD::PrintFor(" The value of h2 is: ", h2[0]);
    CppAD::PrintFor("The value of gradh2 is: ", gradh2[0]);
    for(size_t i=1; i<gradh2.size(); i++){
      CppAD::PrintFor(" ", gradh2[i]);
    }

    // taping ll (log likelihood) store operation sequence
    CppAD::ADFun<a1type> lltape;
    lltape = tapell(zbeta, llf, fromM, logdetJfromM);

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
       lapl[0] += Pmat.row(i) * dPmatfun(z, i) * jac;
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
    smofun.Dependent(ubeta, smo);
    smofun.optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.
    smofun.check_for_nan(false); //no error if some of the results of the Jacobian are nan.
    return(smofun);
}
