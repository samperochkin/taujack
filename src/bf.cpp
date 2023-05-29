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
  bool all_less, all_greater;
  
  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      all_less = true;
      all_greater = true;
      for (size_t k = 0; k < d; k++) {
        if (X[i][k] < X[j][k]) {
          all_greater = false;
        } else if (X[i][k] > X[j][k]) {
          all_less = false;
        }
        if (!all_less && !all_greater) {
          break;
        }
      }
      if (all_less || all_greater) {
        C[i] += 1;
        C[j] += 1;
      }
    }
  }
}

 

// [[Rcpp::export]]
NumericVector bruteForce(NumericMatrix X) {
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

  // convert C and return it
  return Rcpp::wrap(C);
}