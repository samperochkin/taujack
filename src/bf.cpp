#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

//* 
//* Functions to "naively" compute concordance
//* 

void conquer(const std::vector<std::vector<double>>& X,
             std::vector<int>& C) {
  int n = X.size();
  int d = X[0].size();
  bool conc;
  
  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      conc = true;
      if(X[i][0] < X[j][0]){
        for (size_t k = 1; k < d; k++){
          if(X[i][k] > X[j][k]){
            conc = false;
            break;
          }
        }  
      }else{
        for (size_t k = 1; k < d; k++){
          if(X[i][k] < X[j][k]){
            conc = false;
            break;
          }
        }  
      }
      if (conc) {
        C[i] += 1;
        C[j] += 1;
      }
    }
  }
}

 

// [[Rcpp::export]]
IntegerVector bruteForce(NumericMatrix X) {
  int n = X.nrow();
  int d = X.ncol();
  
  // initialize C
  std::vector<int> C(n, 0);
  
  // initialize X_vec
  std::vector<std::vector<double>> X_vec(n, std::vector<double>(d));
  
  // fill X_vec
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      X_vec[i][j] = X(i, j);
    }
  }  
  
  conquer(X_vec, C);

  IntegerVector C_out(n);
  for (int i = 0; i < n; i++) {
    C_out(i) = C[i];
  }
  
  return C_out;
}