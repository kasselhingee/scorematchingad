a1type dirich(const veca1 &u, const veca1 &beta) {
  size_t d  = u.size();
  a1type y(0.);  // initialize summation at 0
  for(size_t i = 0; i < d; i++)
  {   y   += beta[i] * CppAD::log(u[i]);
  }
  return y;
}
