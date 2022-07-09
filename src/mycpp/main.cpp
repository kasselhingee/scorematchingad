# include "scm.cpp"

int main(int argc, char** argv)
 {   using CppAD::AD;   // use AD as abbreviation for CppAD::AD
       std::cout << "Reading inputs" << std::endl;
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

       std::cout << "Creating manifold object" << std::endl;
       manifold<a1type> * man;
       man = new Spos<a1type>();

   std::cout << "Preparing to tape" << std::endl;
   CppAD::ADFun<double> out; //returning a pointer
   out = tapell(z_ad,
                 theta_ad,
                 ll_ppi,
                 man,
                 fixedtheta_e,
                 true);

   delete man;

   return 0;
 }
