#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

//* 
//* Functions to "naively" compute concordance
//* 

void conquer(const std::vector<std::vector<double>>& X,
             std::vector<std::vector<int>>& C) {
  int n = X.size();
  int d = X[0].size();

  int ind;
  int temp;
  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      ind = 0;
      temp = 1;
      if(X[i][0] < X[j][0]){
        for (size_t k = 1; k < d; k++) {
          ind += (X[i][k] < X[j][k]) * temp;
          temp *= 2;
        }
      }else{
        for (size_t k = 1; k < d; k++) {
          ind += (X[i][k] > X[j][k]) * temp;
          temp *= 2;
        }
      }
      C[i][ind] += 1;
      C[j][ind] += 1;
    }
  }
}

// [[Rcpp::export]]
IntegerMatrix bruteForceFull(NumericMatrix X) {
  int n = X.nrow();
  int d = X.ncol();
  
  // initialize C
  int m = static_cast<int>(std::pow(2, d-1));
  std::vector<std::vector<int>> C(n, std::vector<int>(m));
  
  // initialize X_vec
  std::vector<std::vector<double>> X_vec(n, std::vector<double>(d));
  
  // fill X_vec
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      X_vec[i][j] = X(i, j);
    }
  }  
  
  conquer(X_vec, C);

  // convert C and return it as IntegerMatrix
  IntegerMatrix result(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      result(i, j) = C[i][j];
    }
  }
  
  return result;
}
