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
        C[i][0] += 1;
        C[j][0] += 1;
      }
    }
  }
}

void conquer_seq(const std::vector<std::vector<double>>& X,
                          std::vector<std::vector<int>>& D) {
  int n = X.size();
  int d = X[0].size();
  int ind;

  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      ind = -1;
      if(X[i][0] < X[j][0]){
        for (size_t k = 1; k < d; k++){
          if(X[i][k] > X[j][k]){
            ind = k;
            break;
          }
        }  
      }else{
        for (int k = 1; k < d; k++){
          if(X[i][k] < X[j][k]){
            ind = k;
            break;
          }
        }  
      }
      if (ind != -1) {
        D[i][ind-1] += 1;
        D[j][ind-1] += 1;
      }
    }
  }
}


void conquer_all(const std::vector<std::vector<double>>& X,
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
IntegerVector bruteForce(NumericMatrix X, bool seq = false, bool all = false) {
  int n = X.nrow();
  int d = X.ncol();
  
  // initialize X_vec
  std::vector<std::vector<double>> X_vec(n, std::vector<double>(d));
  
  // fill X_vec
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      X_vec[i][j] = X(i, j);
    }
  }  
  
  int K;
  if(all){
    K = static_cast<int>(std::pow(2, d-1));
  } else if(seq){
    K = d-1;
  } else{
    K = 1;
  }
  
  // initialize C or D
  std::vector<std::vector<int>> CoD(n, std::vector<int>(K,0));
  
  if(all){
    conquer_all(X_vec, CoD);
  }else if(seq){
    conquer_seq(X_vec, CoD);
  }else{
    conquer(X_vec, CoD);
  }
    
  IntegerMatrix C_out(n,K);
  if(seq & !all){
    for (int i = 0; i < n; i++) {
      C_out(i, 0) = CoD[i][0];
      for (int j = 1; j < K; j++) {
        C_out(i, j) = C_out(i, j-1) + CoD[i][j];
      }
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < K; j++) {
        C_out(i, j) = n - 1 - C_out(i, j);
      }
    }
  } else{
    for (int i = 0; i < n; i++) {
      for(int j=0; j < K; j++){
        C_out(i,0) = CoD[i][0];
      }
    }
  }
  return C_out;
}