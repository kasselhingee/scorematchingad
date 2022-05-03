# include "scm.cpp"

int main(int argc, char** argv)
 {   using CppAD::AD;   // use AD as abbreviation for CppAD::AD
       //read inputs
       size_t n = 3;
       veca1 z_ad(n);
       for(int i=0; i < n; i++){
           z_ad[i] = 0.1;
       }
       veca1 theta_ad(8);
       for(int i=0; i < theta_ad.size(); i++){
           theta_ad[i] = 1. + i/2.;
       }

       Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta_e(theta_ad.size());
       for(int i=0; i < fixedtheta_e.size() - 2; i++){
           fixedtheta_e[i] = 0;
       }
       fixedtheta_e[fixedtheta_e.size() - 1] = 1;
       fixedtheta_e[fixedtheta_e.size() - 2] = 1;

       manifold<a1type> sphere = {
           Spos::toS, Spos::Pmat_S, Spos::dPmat_S,
           Spos::fromS, Spos::logdetJ_fromS,
       };

   std::cout << "Preparing to tape" << std::endl;
   CppAD::ADFun<double> out; //returning a pointer
   out = tapell(z_ad,
                 theta_ad,
                 ll_ppi,
                 sphere.fromM, //transformation from manifold to simplex
                 sphere.logdetJfromM, //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
                 fixedtheta_e,
                 true);

    return 0;
 }
