# include "PrintFor.cpp"

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

  //given a vector x, return a vector of indices that give
  //x in increasing order
  //with the whole thing taped
  template <class T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> incorder(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x){
    size_t xsize(x.size());
    //A matrix for storing the results of the lt comparison
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ltmat(xsize, xsize);
    ltmat.setZero();
    T zero(0), one(1);
    for (size_t col=0; col<xsize; col++){
      for (size_t row=col+1; row<xsize; row++){
        ltmat(row, col) = CondExpLt(x[col], x[row], one, zero); //if the x[col] item is strictly less than x[row] then ltmat(row,col) is 1
        ltmat(col, row) = 1 - ltmat(row, col); //because the opposite is false
      }
    }
    //index order stored in below object
    Eigen::Matrix<T, Eigen::Dynamic, 1> order(xsize);
    order = ltmat.colwise().sum();
    //a high value in this vector means the corresponding x value is LOWER than many of the other x values
    //order[i] == xsize-1 means x[i] is the LOWEST
    //order[i] == 0 means x[i] is the HIGHEST

    Eigen::Matrix<T, Eigen::Dynamic, 1> order2(xsize);
    order2 = order.array() * (-1) + xsize - 1;
    //xsize-1-order[i] == xsize-1 means x[i] is the highest
    //xsize-1-order[i] == 0 means x[i] is the lowest
    return(order2);
  }

  template <class T>
  T ll_Rivest(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
          const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta){
    //assume the first part of theta defines the matrix A
    // then the next element is concentration parameter, and
    // the final element is the eigenvector/value to use
    size_t d  = u.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Amat(d, d);
    Eigen::Matrix<T, Eigen::Dynamic, 1> Btheta;
    Btheta = theta.block(0,0, d - 1 + (d-1)*d/2, 1);
    Amat = BinghamMatrix(Btheta);
    Eigen::Matrix<T, 1, 1> out_e;
    out_e = u.transpose() * Amat * u;
    T out(out_e[0]);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> eigensolver(Amat);
    if (eigensolver.info() != Eigen::Success) return(20 * abs(out));
    //According to https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html,
    //as the size isn't known at compile time:
    //This SelfAdjointEigenSolver uses a symmetric QR algorithm. The matrix is first reduced
    //to tridiagonal form using the Tridiagonalization class. The tridiagonal matrix is then
    //brought to diagonal form with implicit symmetric QR steps with Wilkinson shift. Details
    //can be found in Section 8.3 of Golub & Van Loan, Matrix Computations.
    // PrintForMatrix("\n Amat is: ", Amat);
    // PrintForMatrix("\n The eigenvectors of Amat are: ", eigensolver.eigenvectors());
    // PrintForVec("\n The eigenvalues of Amat are: ", eigensolver.eigenvalues());

    Eigen::Matrix<T, Eigen::Dynamic, 1> evals;
    evals = eigensolver.eigenvalues().cwiseAbs();//eigenvalues() presents the results in increasing order (negative -> 0 -> positive)
    // PrintForVec("\n The eigenvalues sizes of Amat are: ", evals);

    //ordering
    Eigen::Matrix<T, Eigen::Dynamic, 1> evalorder;
    evalorder = incorder(evals);

    /////////////////choose the largest vector////////////
    //A hack using ones and zeros, could maybe use VecAD in the future
    T zero(0), one(1);
    Eigen::Matrix<T, Eigen::Dynamic, 1> eselector(d); //selector vector in eigen type
    T sizethwanted;
    sizethwanted = theta[Btheta.size() + 1] - 1;
    for (size_t i=0; i<d; i++){
      eselector[i] = CondExpLt(evalorder[i], sizethwanted + 0.1, one, zero); //if less than 0.1 + idx wanted then one
      eselector[i] = eselector[i] * CondExpLt(evalorder[i] + 1, sizethwanted + 0.1, zero, one); //if 1+ is also less than index wanted then multiple by zero
    }
    // PrintForVec("\n eselector is: ", eselector);

    Eigen::Matrix<T, Eigen::Dynamic, 1> m;
    m = eigensolver.eigenvectors() * eselector;
    // PrintForVec("\nThe original m is: ", m);
    //force first element of m to be positive, not sure how well this works for differentiation later
    T multiplier(1.);
    multiplier = CondExpLt(m[0], zero, -one, one);
    m = multiplier * m;

    // PrintForVec("\n m is: ", m);
    ////////////finished getting the eigenvector////////////


    out_e += theta[Btheta.size()] * m.transpose() * u;

    out = out_e[0];
    return(out);
  }

}
