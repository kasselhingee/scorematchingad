# include "tapesmd.h"

// function that tapes the score-matching objective
pADFun tapesmd(pADFun & uldtape,
               transform<a1type> &tran,
               manifold<a1type> &M,
               std::string weightname, //the weight function h^2 - name of a a1type (*h2fun)(const veca1 &, const double &)
               const double & acut, //the acut constraint for the weight functions
               bool verbose
               ){
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

    size_t d(uldtape.Domain()); //manifold is embedded in Euclidean space of dimension d (may be different to simplex)

    //Projection matrix
    mata1 Pmat(d, d);

    //get boundary weight function tape
    CppAD::ADFun<double> dh2tape;
    veca1 z=uldtape.xtape; //z is on the reembed manifold, uldtape should already be this
    CppAD::ADFun<a1type, double> h2tape; //The second type here 'double' is for the 'RecBase' in ad_fun.hpp. It doesn't seem to change the treatment of the object.
    h2tape = tapeh2(z, h2fun, acut).base2ad(); //convert to a function of a1type rather than double

    //check inputs and ll tape match
    if (uldtape.Domain() != (unsigned)z.size()){Rcpp::stop("Dimension of z (the input measurement on the manifold) is %i, which does not match domain size of log density function of %i.", z.size(), uldtape.Domain());}

    // make ll tape higher order (i.e. as if it was taped using a2type instead of a1type)
    CppAD::ADFun<a1type, double> uldtapehigher;
    uldtapehigher = (uldtape.get_ptr())->base2ad();



    //START TAPING
    veca1 theta = uldtape.dyntape;
    veca1 u = tran.fromM(uldtape.xtape); //xtape is on end manifold, but want smo to use values in start manifold
    CppAD::Independent(theta, u); //because looking at empirical divergence the main arguments are parameter theta, and the dynamic 'parameter' is u
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

    //update parameters ('dynamic' values) of uldtape
    uldtapehigher.new_dynamic(theta);

    //grad(ll)
    veca1 jac(d); // Jacobian of ll
    jac  = uldtapehigher.Jacobian(z);      // Jacobian for operation sequence
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
    hess = uldtapehigher.Hessian(z, 0); //the zero here is something about selecting the range-space component of f, 0 selects the first and only component, I presume.
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

    pADFun out(smdfun, theta, u, "smo:" + uldtape.name);
    return(out);
}


