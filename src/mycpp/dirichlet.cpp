// # include <cppad/example/atomic_three/mat_mul.hpp> // for matrix multiplication tapesmo's matrix multiplication don't seem to need this though :S

namespace { // begin the empty namespace

    template <class T>
    T ll_dirichlet(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
	       const Eigen::Matrix<T, Eigen::Dynamic, 1> &beta)
    {   size_t d  = u.size();
        T y(0.);  // initialize summation
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }

    template <class T>
    T ll_ppi(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta){
        size_t d  = u.size();
        //assume the parameter vector theta is encoded as:
        //c(diag(ALs), ALs[upper.tri(ALs)], bL, beta)
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Amat(d, d);
        Amat.setZero();
        // Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ALmat(d-1, d-1);
        //populate the upper and lower triangles
        //the upper triangle has (d-1)-1 rows, the rows have 1 - ()d-2) elements. Arithmetic series:
        //(d-2) * (2 + d - 3) / 2 = (d-2) * (d-1)/2
        Eigen::Matrix<T, Eigen::Dynamic, 1> upptriblock((d-2) * (d-1)/2);
        upptriblock = theta.block(d-1, 0, upptriblock.size(), 1);
        size_t vecidx = 0;
        for (size_t col=1; col < d-1; col++){
          for (size_t row=0; row < col; row++){
            Amat(row, col) = upptriblock(vecidx, 1);
            Amat(col, row) = upptriblock(vecidx, 1);
            vecidx +=1;
          }
        }
        //populate the diagonal
        vecidx = 0;
        for (size_t row=0; row < d-1; row++){
            Amat(row,row) = theta(vecidx, 1);
            vecidx +=1;
        }
        // ALmat.diagonal() << theta.block(0,0,d-1,1);
        // Amat.block(0,0,d-1,d-1) = ALmat;
        std::cout << "A matrix:" << std::endl;
        std::cout << Amat << std::endl;

        Eigen::Matrix<T, Eigen::Dynamic, 1> bvec(d);
        bvec.setZero();
        bvec.block(0,0,d-1,1) = theta.block(d-1 + upptriblock.size(), 0, d-1, 1);
        std::cout << "b vector:" << std::endl;
        std::cout << bvec.transpose() << std::endl;

        Eigen::Matrix<T, Eigen::Dynamic, 1> beta(d);
        beta = theta.block(theta.size() - d, 0, d,1);
        std::cout << "beta vector:" << std::endl;
        std::cout << beta.transpose() << std::endl;

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

