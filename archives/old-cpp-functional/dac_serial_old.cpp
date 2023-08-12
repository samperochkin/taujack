#include <Rcpp.h>
#include <vector>

using namespace Rcpp;


//* 
//* Functions to "naively" compute concordance
//* 

void conquer(const std::vector<std::vector<double>>& X,
             std::vector<std::vector<int>>& D,
             const std::vector<int>& ids) {
  int n = ids.size();
  int d = X[0].size();

  // Make sure X[,1] is sorted
  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      for (size_t k = 1; k < d; k++){
        if(X[ids[i]][k] > X[ids[j]][k]){
          D[ids[i]][k-1] += 1;
          D[ids[j]][k-1] += 1;
          break;
        }
      }  
    }
  }
}


void conquer1(const std::vector<std::vector<double>>& X,
              std::vector<std::vector<int>>& D,
              const std::vector<int>& ids0,
              const std::vector<int>& ids1,
              const int k0) {
  int n0 = ids0.size();
  int n1 = ids1.size();
  int d = X[0].size();

  // Make sure X[][0] is sorted
  for (int i : ids0) {
    for (int j : ids1) {
      for (int k = k0; k < d; k++) {
        if (X[i][k] > X[j][k]){
          D[i][k-1] += 1;
          D[j][k-1] += 1;
          break;
        }
      }
    }
  }
}


//* 
//* Main functions
//* 

void divideAndConquer1(const std::vector<std::vector<double>>& X,
                           std::vector<std::vector<int>>& D,
                           const int thresh,
                           std::vector<int>& ids0,
                           std::vector<int>& ids1,
                           const int k0) {

  int n0 = ids0.size();
  int n1 = ids1.size();
  int d = X[0].size();

  // k0 = d means that the given pairs are concordant
  if(k0 == d) return;

  double sqrt_nn = std::sqrt(n0)*std::sqrt(n1);

  if(sqrt_nn <= thresh || n0 < 3 || n1 < 3) {
    // If n is small enough, perform brute force computation
    conquer1(X, D, ids0, ids1, k0);
    return;
  }

  // Otherwise, split Ids into 2 sets relative to the midpoint
  double mid;
  double mi;
  double ma;
  mi = std::numeric_limits<double>::infinity() ;
  ma = std::numeric_limits<double>::infinity()*(-1);
  for (int i : ids0) {
    mi = std::min(mi, X[i][k0]);
    ma = std::max(ma, X[i][k0]);
  }
  for (int i : ids1) {
    mi = std::min(mi, X[i][k0]);
    ma = std::max(ma, X[i][k0]);
  }
  mid = (ma + mi)/2;

  std::vector<int> ids00;
  std::vector<int> ids01;
  std::vector<int> ids10;
  std::vector<int> ids11;
  for (int i : ids0) {
    if(X[i][k0] <= mid){
      ids00.push_back(i);
    }else{
      ids01.push_back(i);
    }
  }
  for (int i : ids1) {
    if(X[i][k0] <= mid){
      ids10.push_back(i);
    }else{
      ids11.push_back(i);
    }
  }

  int n01 = ids01.size();
  int n10 = ids10.size();

  // necessarily discordant
  for(int i : ids01) D[i][k0-1] += n10;
  for(int j : ids10) D[j][k0-1] += n01;

  // need further investigation
  divideAndConquer1(X, D, thresh, ids00, ids10, k0);
  divideAndConquer1(X, D, thresh, ids01, ids11, k0);

  // can move on to next dimension
  divideAndConquer1(X, D, thresh, ids00, ids11, k0+1);
}


// no need to compute midpoint, but need to know zeros and ones (zo)
void divideAndConquer2(const std::vector<std::vector<double>>& X,
                       std::vector<std::vector<int>>& D,
                       const int thresh,
                       std::vector<int>& ids,
                       std::vector<int>& zo) {

  int n = ids.size();
  int d = X[0].size();

  if(n <= thresh) {
    conquer(X, D, ids);
    return;
  }

  std::vector<int> ids0;
  std::vector<int> ids1;
  std::vector<int> zo0;
  std::vector<int> zo1;
  std::vector<int> ids00;
  std::vector<int> ids01;
  std::vector<int> ids10;
  std::vector<int> ids11;
  for(int i=0; i < n/2; i++){
    ids0.push_back(ids[i]);
    zo0.push_back(zo[i]);
    if(zo[i] == 0){
      ids00.push_back(ids[i]);
    }else{
      ids01.push_back(ids[i]);
    }
  }
  for(int i = n/2; i < n; i++){
    ids1.push_back(ids[i]);
    zo1.push_back(zo[i]);
    if(zo[i] == 0){
      ids10.push_back(ids[i]);
    }else{
      ids11.push_back(ids[i]);
    }
  }

  int n01 = ids01.size();
  int n10 = ids10.size();

  // necessarily discordant
  for(int i : ids01) D[i][0] += n10;
  for(int j : ids10) D[j][0] += n01;

  // need further investigation (within blocks)
  divideAndConquer2(X, D, thresh, ids0, zo0);
  divideAndConquer2(X, D, thresh, ids1, zo1);

  // need further investigation (between blocks)
  divideAndConquer1(X, D, thresh, ids00, ids10, 1);
  divideAndConquer1(X, D, thresh, ids01, ids11, 1);
  
  // can move on to next dimension
  divideAndConquer1(X, D, thresh, ids00, ids11, 2);
}


void divideAndConquer(const std::vector<std::vector<double>>& X,
                      std::vector<std::vector<int>>& D,
                      const int thresh,
                      const std::vector<int> ids) {

  int n = ids.size();
  int d = X[0].size();

  if(n <= thresh) {
    conquer(X, D, ids);
    return;
  }

  // Split Ids into sets relative to the midpoint
  double mid;
  double mi;
  double ma;
  mi = std::numeric_limits<double>::infinity() ;
  ma = std::numeric_limits<double>::infinity()*(-1);
  for (int i : ids) {
    mi = std::min(mi, X[i][1]);
    ma = std::max(ma, X[i][1]);
  }
  mid = (mi + ma)/2.0;

  std::vector<int> ids0;
  std::vector<int> ids1;
  std::vector<int> zo0;
  std::vector<int> zo1;
  std::vector<int> ids00;
  std::vector<int> ids01;
  std::vector<int> ids10;
  std::vector<int> ids11;
  for(int i=0; i < n/2; i++){
    ids0.push_back(ids[i]);
    if(X[ids[i]][1] <= mid){
      zo0.push_back(0);
      ids00.push_back(ids[i]);
    }else{
      zo0.push_back(1);
      ids01.push_back(ids[i]);
    }
  }
  for(int i = n/2; i < n; i++){
    ids1.push_back(ids[i]);
    if(X[ids[i]][1] <= mid){
      zo1.push_back(0);
      ids10.push_back(ids[i]);
    }else{
      zo1.push_back(1);
      ids11.push_back(ids[i]);
    }
  }

  int n01 = ids01.size();
  int n10 = ids10.size();

  // necessarily discordant
  for(int i : ids01) D[i][0] += n10;
  for(int j : ids10) D[j][0] += n01;

  // need further investigation (within blocks)
  // divideAndConquer(X, D, thresh, ids0);
  // divideAndConquer(X, D, thresh, ids1);
  divideAndConquer2(X, D, thresh, ids0, zo0);
  divideAndConquer2(X, D, thresh, ids1, zo1);
  
  // need further investigation (between blocks)
  divideAndConquer1(X, D, thresh, ids00, ids10, 1);
  divideAndConquer1(X, D, thresh, ids01, ids11, 1);

  // can move on to next dimension
  divideAndConquer1(X, D, thresh, ids00, ids11, 2);
}


// [[Rcpp::export]]
IntegerMatrix dac_serial(NumericMatrix X, int thresh = 100, bool brute_force = false) {
  int n = X.nrow();
  int d = X.ncol();

  // initialize D
  std::vector<std::vector<int>> D(n, std::vector<int>(d - 1, 0));
  
  // initialize X_vec
  std::vector<std::vector<double>> X_vec(n, std::vector<double>(d));

  // fill X_vec
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      X_vec[i][j] = X(i, j);
    }
  }

  // Reorder X wrt first column
  std::vector<std::pair<double, int>> val_ind;
  for (int i = 0; i < n; i++) {
    val_ind.push_back({X_vec[i][0], i});
  }
  std::sort(val_ind.begin(), val_ind.end());

  std::vector<int> ids(n);
  for (int i = 0; i < n; i++) {
    ids[i] = val_ind[i].second;
  }

  if(brute_force){
    conquer(X_vec, D, ids);
  }else{
    divideAndConquer(X_vec, D, thresh, ids);
  }

  // Convert D to Rcpp::IntegerMatrix and return it
  Rcpp::IntegerMatrix D_matrix(n, d - 1);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d - 1; j++) {
      D_matrix(i, j) = D[i][j];
    }
  }
  
  return D_matrix;
}

