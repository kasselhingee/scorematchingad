# include "scm.cpp"

int main(int argc, char** argv)
 {   using CppAD::AD;   // use AD as abbreviation for CppAD::AD
       //read inputs
       size_t n = 3;
       vecd xin(n);
       vecd betain(n);
       for(int i=0; i < n; i++){
            xin[i] = 3.;
            betain[i] = 1.;
         }

     CppAD::ADFun<double> smotape;
     svecd xbetain(2 * n);
     for(int i=0; i < n; i++){
             xbetain[i] = xin[i];
             xbetain[i + n] = betain[i];
     }
     smotape = tapesmo(xbetain, n, ll,
                       Spos::toS, Spos::Pmat_S, Spos::dPmat_S,
                       Spos::fromS, Spos::logdetJ_fromS,
                       minsq, gradminsq, 0.1);

    //accept command line arguments for evaluation
       if (argc > 1){
            for(int i=0; i < n; i++){
                xin[i] = strtod(argv[i + 1], NULL);
                betain[i] = strtod(argv[i + n + 1], NULL);
              }
         }

     //update xbetain
       for(int i=0; i < n; i++){
               xbetain[i] = xin[i];
               xbetain[i + n] = betain[i];
       }
       std::cout << "x in is " << xin << std::endl;
       std::cout << "beta in is " << betain << std::endl;
       // std::cout << "h2 is " << minsq(xin, 1.) << std::endl;
       // std::cout << "grad(h2) is " << gradminsq(xin, 1.) << std::endl;

        svecd smo_val(1);
        smo_val = smotape.Forward(0, xbetain);
        std::cout << "smo is: " << smo_val[0] << std::endl;

         return 0;
 }
