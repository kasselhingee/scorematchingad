namespace { // begin the empty namespace

  template <class T>
  T ll_vMF(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
                 const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta)
    //lklhood is log(k*mu.x), where mu is a unit vector
    //for this software theta is the vector k*mu, which is unrestricted (except mu shouldn't be zero)
  {
    T y(0.);  // initialize summation
    for(size_t i = 0; i < u.size(); i++)
    {   y   += theta[i] * u[i];
    }
    return y;
  }

  template <class T>
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
  BinghamMatrix(const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta){
    // A is symmetric, and the sum of the diagonals is zero
    //assume the parameter vector theta is encoded as:
    //c(diag(A)[1:(p-1)], A[upper.tri(ALs)]
    size_t d = (-1 + std::sqrt(1 + 4 * 2 * (1+theta.size()))) / 2 + 0.5;//the +0.5 makes sure truncation gets to the correct integer
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Amat(d, d);
    Amat.setZero();
    //populate the diagonal
    size_t vecidx = 0;
    for (size_t row=0; row < d-1; row++){
      Amat(row,row) = theta[vecidx];
      vecidx +=1;
    }
    Amat(d-1, d-1) = -theta.block(0, 0, d - 1, 1).sum(); //so trace of A is zero

    //populate the upper and lower triangles
    //the upper triangle has d-1 rows, the rows have 1 to (d-1) elements. Arithmetic series:
    //(d-1)/2 [2 + (d-1−1)] = (d-1) (2 + d - 2)/2 = (d - 1) d/2
    Eigen::Matrix<T, Eigen::Dynamic, 1> upptriblock((d - 1) * d/2);
    upptriblock = theta.block(d-1, 0, upptriblock.size(), 1);
    vecidx = 0;
    for (size_t col=1; col < d; col++){
      for (size_t row=0; row < col; row++){
        Amat(row, col) = upptriblock[vecidx]; //bug fix - column index is 0!!
        Amat(col, row) = upptriblock[vecidx];
        vecidx +=1;
      }
    }
    return(Amat);
  }

  template <class T>
  T ll_Bingham(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
           const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta)
    //lklhood is log(u*A*u) - from Mardia et al 2016
    // A is symmetric, and the sum of the diagonals is zero
  {
    size_t d  = u.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Amat(d, d);
    Amat = BinghamMatrix(theta);

    Eigen::Matrix<T, 1, 1> out_e;
    out_e = u.transpose() * Amat * u;
    T out(out_e[0]);
    return out;
  }

  template <class T>
  T ll_FB(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
               const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta)
    //lklhood is log(u*A*u + km*u) - from Mardia et al 2016
    // A is symmetric, and the sum of the diagonals is zero
    //m*k is any vector
  {
    size_t d  = u.size();
    size_t Binghamthetasize = d - 1 + (d-1)*d/2;
    Eigen::Matrix<T, Eigen::Dynamic, 1> Btheta;
    Btheta = theta.block(0,0, Binghamthetasize, 1);
    Eigen::Matrix<T, Eigen::Dynamic, 1> Ftheta;
    Ftheta = theta.block(Binghamthetasize,0, d, 1);
    T out;
    out = ll_Bingham(u, Btheta);
    out += ll_vMF(u, Ftheta);
    return(out);
  }

  template <class T>
  T ll_Rivest(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
          const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta){
    //assume the first part of theta defines the matrix A
    // then the next element is concentration parameter, and
    // the final element is the eigenvector/value to use
    std::cout << "Theta is: " << theta.transpose() << std::endl;
    size_t d  = u.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Amat(d, d);
    Eigen::Matrix<T, Eigen::Dynamic, 1> Btheta;
    Btheta = theta.block(0,0, d - 1 + (d-1)*d/2, 1);
    std::cout << "Btheta is: " << Btheta.transpose() << std::endl;
    Amat = BinghamMatrix(Btheta);
    Eigen::Matrix<T, 1, 1> out_e;
    out_e = u.transpose() * Amat * u;
    T out(out_e[0]);
    std::cout << "Out: " << out_e << std::endl;


    std::cout << "Amat is:" << std::endl << Amat << std::endl;
    std::cout << "k is:" << theta.block(Btheta.size(), 0, 1, 1) << std::endl;
    std::cout << "idx is:" << theta.block(Btheta.size() + 1, 0, 1, 1) << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> eigensolver(Amat);
    if (eigensolver.info() != Eigen::Success) return(out);
    std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
              << "corresponding to these eigenvalues:\n"
              << eigensolver.eigenvectors() << std::endl;

    Eigen::Matrix<T, Eigen::Dynamic, 1> evals;
    evals = eigensolver.eigenvalues().cwiseAbs();
    std::cout << "The eigenvalues sizes of Amat are:\n" << evals.transpose() << std::endl;
    std::vector<T> sevals(evals.size());
    for (size_t i=0; i<evals.size(); i++){
      sevals[i] = evals[i];
    }
    std::vector<size_t> ind(evals.size());
    CppAD::index_sort(sevals, ind);
    std::cout << "eval order: " << ind[0] << std::endl;
    std::cout << "eval order: " << ind[1] << std::endl;
    std::cout << "eval order: " << ind[2] << std::endl;

    Eigen::Matrix<T, Eigen::Dynamic, 1> m;
    // CppAD::AD<int> eidx;
    // eidx = theta[Btheta.size() + 1] + 0.5;
    m = eigensolver.eigenvectors().block(0, ind[0], d, 1);
    std::cout << "m is: " << m << std::endl;
    //extra: put a condition that the first element of m is negative (or positive) so that don't get directional uncertainty
    //but I suspect this is already fixed by eigen

    out_e += theta[Btheta.size()] * m.transpose() * u;
    std::cout << "Out: " << out_e << std::endl;

    out = out_e[0];
    return(out);
  }

}
