
namespace { // begin the empty namespace

    template <class T>
    T ll_dirichlet(const Eigen::Matrix<T, Eigen::Dynamic, 1> &beta,
	       const Eigen::Matrix<T, Eigen::Dynamic, 1> &u)
    {   size_t d  = u.size();
        T y(0.);  // initialize summation
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }

    template <class T>
    T ll_ppi(const Eigen::Matrix<T, Eigen::Dynamic, 1> &beta,
             const Eigen::Matrix<T, Eigen::Dynamic, 1> & u){
        size_t d  = u.size();
        //assume the parameter vector beta is encoded as:
        //c(diag(ALs), ALs[upper.tri(ALs)], bL, betafordirichlet)
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Amat(d, d);
        Amat.setZero();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ALmat(d-1, d-1);
        ALmat.diagonal() = beta.block(0,0,d-1,1);
        //the upper triangle has (d-1)-1 rows, the rows have 1 - ()d-2) elements. Arithmetic series:
        //(d-2) * (2 + d - 3) / 2 = (d-2) * (d-1)/2
        size_t upptrisize = (d-2) * (d-1)/2;
        ALmat.triangularView<Eigen::Upper>() = beta.block(d-1, 0, upptrisize, 1);
        ALmat.triangularView<Eigen::Lower>() = beta.block(d-1, 0, upptrisize, 1);
        Amat.block(0,0,d-1,d-1) = ALmat;
        Eigen::Matrix<T, Eigen::Dynamic, 1> bvec(d);
        bvec.setZero();
        bvec.block(0,0,d-1,1) = beta.block(d-1 + upptrisize, 0, d-1, 1);

        T y(0.);  // initialize summation
        y += u.transpose()*Amat*u;
        y += bvec.transpose()*u;
        //dirichlet component
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }


}

