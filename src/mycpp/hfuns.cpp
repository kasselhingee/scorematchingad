  //////////////////////////////////////////
  // weight function and grad(h^2) functions
  // squared without constraint
  template <class Type>
  Type prodsq(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
    Type out;
    out = x.array().square().prod();
    return(out);
  }

  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> gradprodsq(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x){
    size_t n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(n);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> avoidone(n-1);
    Type prodx;
    prodx = 2 * x.array().prod();
    for (int i=0; i < n; i++){
	avoidone << x.head(i), x.tail(n-i-1);
	out[i] = prodx * avoidone.prod();
    }
    return(out);
  }
