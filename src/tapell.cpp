# include "tapell.hpp"

// for calculating the determinant of the transform to a manifold
CppAD::ADFun<double> tapefromM(veca1 z,
                               manifold<a1type> *pman){
  //tape relationship between z and h2
  CppAD::Independent(z);
  // range space vector
  veca1 y(0); // vector of ranges space variables - length to be set by from(z);
  y = pman->fromM(z);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  tape.optimize(); //remove some of the extra variables that were used for recording the ADFun, but aren't needed anymore (hopefully very good when logdetJ is constant.
  return(tape);
}



// define a function that tapes a log likelihood
CppAD::ADFun<double> tapell(veca1 z, //data measurement tranformed to M manifold
                            veca1 theta, //theta parameter
                               a1type (*llf)(const veca1 &, const veca1 &), //the log likelihood function
                               manifold<a1type> *pman, //it seems points must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)
                               Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)
                               bool verbose
                               ){
  if (theta.size() != fixedtheta.size()){
    Rcpp::stop("theta and fixedtheta must have the same length");
  }

  //separate fixed and variable theta
  veca1 thetavar(theta.size() - fixedtheta.sum());
  veca1 thetafxd(fixedtheta.sum());
  size_t idx_var(0);
  size_t idx_fxd(0);
  for (long int i=0;i<theta.size();i++){
    if (fixedtheta[i]){
      thetafxd[idx_fxd] = theta[i];
      idx_fxd += 1;
    } else {
      thetavar[idx_var] = theta[i];
      idx_var += 1;
    }
  }

  if (verbose){
    Rcpp::Rcout << "Fixing according to pattern: " << std::endl;
    for (long int i=0;i<fixedtheta.size();i++){
      Rcpp::Rcout << " " << fixedtheta[i];
    }
    Rcpp::Rcout << std::endl;

    Rcpp::Rcout << "Fixed theta is:";
    if (thetafxd.size() == 0){
      Rcpp::Rcout << " none" << std::endl;
    } else {
      for (long int i=0;i<thetafxd.size();i++){
        Rcpp::Rcout << " " << thetafxd[i];
      }
      Rcpp::Rcout << std::endl;
    }
  }

  //prepare fromM tape for logdetJfromM later
  CppAD::ADFun<a1type, double> fromMtape; //The second type here 'double' is for the 'RecBase' in ad_fun.hpp. It doesn't seem to change the treatment of the object.
  fromMtape = tapefromM(z, pman).base2ad(); //convert to a function of a1type rather than double
  if(fromMtape.Domain() != z.size()){Rcpp::stop("fromMtape domain size does not match input z size: %i != %i", fromMtape.Domain(), z.size());}



  //tape relationship between x and log-likelihood
  CppAD::Independent(z, thetavar);  //for this tape, theta must be altered using new_dynamic
  if (verbose){
    Rcpp::Rcout << "thetavar is: " << thetavar.transpose() << std::endl;
    PrintForVec("\n thetavar is: ", thetavar);
  }

  //combine fixed and variable theta
  veca1 thetarecom(theta.size());
  idx_var = 0;
  for (long int i=0;i<theta.size();i++){
    if (fixedtheta[i]){
      thetarecom[i] = theta[i];
    } else {
      thetarecom[i] = thetavar[idx_var];
      idx_var += 1;
    }
  }
  if (verbose){
    Rcpp::Rcout << "thetarecom is: " << thetarecom.transpose() << std::endl;
    PrintForVec("\n thetarecom is: ", thetarecom);
  }

  // range space vector
  veca1 y(1); // vector of ranges space variables
  veca1 u(0); //0 here because size dictated by fromM
  u = fromMtape.Forward(0, z);
  y.setZero();
  y[0] += llf(u, thetarecom);

  //get log determinant of fromM
  mata1 jacmat(z.size() * u.size(), 1);
  jacmat = fromMtape.Jacobian(z);
  jacmat.resize(z.size(), u.size()); //first row is: du1/dz1, du2/dz1, du3/dz1. Second row is du1/dz2, du2/dz2, du3/dz2
  veca1 logdet(1);
  logdet[0] = CppAD::log(jacmat.determinant());
  y[0] += logdet[0];
  if (verbose){
    PrintForVec("\n z is: ", z);
    PrintForVec("\n fromM(z) is: ", u.transpose());
    PrintForMatrix("\n jacmat is: ", jacmat);
    PrintForVec("\n jacmat logdeterminant is: ", logdet);
  }
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  tape.optimize(); //remove some of the extra variables that were used for recording the ADFun, but aren't needed anymore.
  if (verbose){
    Rcpp::Rcout << "tape has " << tape.size_dyn_ind() << " independent dynamic parameters" << std::endl;
    Rcpp::Rcout << "tape requires vectors of length " << tape.Domain() << std::endl;
    Rcpp::Rcout << "tape returns vectors of length " << tape.Range() << std::endl;
  }
  return(tape);
}


Rcpp::XPtr< CppAD::ADFun<double> > ptapell(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     std::string llname,
                                     Rcpp::XPtr< manifold<a1type> > pman,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     ){

  //choose ll function
  a1type (*ll)(const veca1 &, const veca1 &) = nullptr;
  if (llname.compare("dirichlet") == 0){
    ll = ll::ll_dirichlet;
  }
  if (llname.compare("ppi") == 0){
    ll = ll::ll_ppi;
  }
  if (llname.compare("vMF") == 0){
    ll = ll::ll_vMF;
  }
  if (llname.compare("Bingham") == 0){
    ll = ll::ll_Bingham;
  }
  if (llname.compare("FB") == 0){
    ll = ll::ll_FB;
  }
  if (llname.compare("Rivest") == 0){
    ll = ll::ll_Rivest;
  }
  //check ll function
  if (ll == nullptr){
    throw std::invalid_argument("Matching ll function not found");
  }

  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  *out = tapell(z_ad,
                theta_ad,
                ll,
                pman.checked_get(),
                fixedtheta,
                verbose);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}




