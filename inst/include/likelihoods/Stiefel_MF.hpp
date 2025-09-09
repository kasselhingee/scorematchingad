#ifndef likelihoods_Stiefel_MF
#define likelihoods_Stiefel_MF
#include <RcppEigen.h>

namespace ll { // namespace for log-likelihood functions

    template <class T>
    // The density given by Eq 1.1 of Hoff 2009 that is proportional to exp(tr(M^T X))
    // because of an elegant property of the vec operator
    // tr(M^T X) = vec(M)^T vec(X)
    T ll_Stiefel_MF(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u, //u is vectorised matrix
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta){//theta is a vectorised matrix too
        Eigen::Matrix<T, 1, 1> ymat = theta.transpose() * u; //ymat should be 1 x 1
        return ymat(0,0);
  }
    

} // namespace ll

#endif
