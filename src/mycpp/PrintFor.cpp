template <class T> //T must be something CppAD::PrintFor can print
void PrintForVec(const char* before, const T & printvec){
  CppAD::PrintFor(before, printvec[0]);
  for(size_t i=1; i<printvec.size(); i++){
    CppAD::PrintFor(" ", printvec[i]);
  }
}
