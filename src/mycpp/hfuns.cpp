  //////////////////////////////////////////
  // weight function and grad(h^2) functions

  //prodsq
  template <class Type>
  Type prodsq(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const double & acut){
    Type prd;
    prd = x.array().square().prod();
    //constraint
    Type acutb(acut * acut);
    Type out = CppAD::CondExpLe(prd, acutb, prd, acutb);
    return(out);
  }

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> gradprodsq(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const double & acut){
    size_t n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(n);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> avoidone(n-1);
    Type prodx;
    prodx = 2 * x.array().prod();
    for (size_t i=0; i < n; i++){
    	avoidone << x.head(i), x.tail(n-i-1);
    	out[i] = prodx * avoidone.prod();
    }
    //apply constraint, need to do prodsq again
    Type acutb(acut * acut);
    Type prd;
    prd = x.array().square().prod();
    Type one(1.0);
    Type mult = CppAD::CondExpLe(prd, acutb, one, one * 0.);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> outc = out * mult;
    return(outc);
  }

  //minsq
  template <class Type>
  Type minsq(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const double & acut){
    Type minval;
    minval = x.array().square().minCoeff(); //this version may need retaping!
    //constraint
    Type acutb(acut * acut);
    Type out = CppAD::CondExpLe(minval, acutb, minval, acutb);
    return(out);
  }

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> gradminsq(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const double & acut){
    size_t n = x.size();
    typename Eigen::Matrix<Type, Eigen::Dynamic, 1>::Index min_index;
    Type gradsize;
    gradsize = 2.0 * x.array().minCoeff(&min_index);

    //apply constraint to size
    Type acutb(acut * acut);
    Type one(1.0);
    Type mult = CppAD::CondExpLe(gradsize/2, acutb, one, one * 0.);
    gradsize *= mult;

    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(n);
    out.setZero();
    out(min_index) = gradsize;
    return(out);
  }

  //hprod
  template <class Type>
  Type hprod(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const double & acut){
    Type prd;
    prd = x.array().prod();
    //constraint
    Type acutb(acut);
    Type out = CppAD::CondExpLe(prd, acutb, prd, acutb);
    return(out);
  }

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> gradhprod(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const double & acut){
    size_t n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(n);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> avoidone(n-1);
    for (size_t i=0; i < n; i++){
      avoidone << x.head(i), x.tail(n-i-1);
      out[i] = avoidone.prod();
    }
    //apply constraint, need to do prod again
    Type acutb(acut);
    Type prd;
    prd = prd = x.array().prod();
    Type one(1.0);
    Type mult = CppAD::CondExpLe(prd, acutb, one, one * 0.);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> outc = out * mult;
    return(outc);
  }
