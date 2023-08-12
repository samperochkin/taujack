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
                       const int k0,
                       const bool need_sort) {
  
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
  

  // Otherwise, split each ids sets into 2
  int count0 = 0;
  int count1 = 0;
  int split0;
  int split1;
  int split_quality = n0*n1;
  int new_split_quality;

  if(need_sort) {
    // Reorder ids0 and ids1 wrt variable k0
    std::vector<std::tuple<double, int, int>> val_ind_01;
    for (int i : ids0) val_ind_01.push_back({X[i][k0], i, 0});
    for (int i : ids1) val_ind_01.push_back({X[i][k0], i, 1});

    std::sort(val_ind_01.begin(), val_ind_01.end());
    
    for (const auto& triplet : val_ind_01) {
      
      if(std::get<2>(triplet) == 0){
        ids0[count0] = std::get<1>(triplet);
        count0 += 1;
      }else{
        ids1[count1] = std::get<1>(triplet);
        count1 += 1;
      }
      
      new_split_quality = count0*count1 + (n0 - count0)*(n1 - count1);
      if(new_split_quality <= split_quality){
        split_quality = new_split_quality;
        split0 = count0;
        split1 = count1;
      }
    }
  }else{
    
    while(count0<n0 && count1<n1){
      if(X[ids0[count0]][k0] < X[ids1[count1]][k0]){
        count0 += 1;
      }else{
        count1 += 1;
      }
      new_split_quality = count0*count1 + (n0 - count0)*(n1 - count1);
      if(new_split_quality <= split_quality){
        split_quality = new_split_quality;
        split0 = count0;
        split1 = count1;
      }
    }
  }

  std::vector<int> ids00(ids0.begin(), ids0.begin() + split0);
  std::vector<int> ids01(ids0.begin() + split0, ids0.end());
  std::vector<int> ids10(ids1.begin(), ids1.begin() + split1);
  std::vector<int> ids11(ids1.begin() + split1, ids1.end());
  
  // necessarily discordant
  for(int i : ids01) D[i][k0-1] += split1;
  for(int j : ids10) D[j][k0-1] += n0 - split0;
  
  // need further investigation
  // should already be in the right ordering
  divideAndConquer1(X, D, thresh, ids00, ids10, k0, false);
  divideAndConquer1(X, D, thresh, ids01, ids11, k0, false);
  
  // can move on to next dimension
  divideAndConquer1(X, D, thresh, ids00, ids11, k0+1, true);
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
  
  std::vector<int> ids0(ids.begin(), ids.begin() + n/2);
  std::vector<int> ids1(ids.begin() + n/2, ids.end());
  int n0 = n/2;
  int n1 = n - n/2;
  
  std::vector<int> zo0;
  std::vector<int> zo1;
  std::vector<int> ids00;
  std::vector<int> ids01;
  std::vector<int> ids10;
  std::vector<int> ids11;
  for(int i=0; i < n0; i++){
    zo0.push_back(zo[i]);
    if(zo[i] == 0){
      ids00.push_back(ids0[i]);
    }else{
      ids01.push_back(ids0[i]);
    }
  }
  for(int i = 0; i < n1; i++){
    zo1.push_back(zo[n0 + i]);
    if(zo[n0 + i] == 0){
      ids10.push_back(ids1[i]);
    }else{
      ids11.push_back(ids1[i]);
    }
  }
  
  int n01 = ids01.size();
  int n10 = ids10.size();
  
  // necessarily discordant
  for(int i : ids01) D[i][0] += n10;
  for(int j : ids10) D[j][0] += n01;
  
  divideAndConquer2(X, D, thresh, ids0, zo0);
  divideAndConquer2(X, D, thresh, ids1, zo1);
  
  // need further investigation (between blocks)
  divideAndConquer1(X, D, thresh, ids00, ids10, 1, true);
  divideAndConquer1(X, D, thresh, ids01, ids11, 1, false);
  
  // can move on to next dimension
  divideAndConquer1(X, D, thresh, ids00, ids11, 2, true);
}


void divideAndConquer(const std::vector<std::vector<double>>& X,
                      std::vector<std::vector<int>>& D,
                      const int thresh,
                      const std::vector<int>& ids) {
  
  int n = ids.size();
  int d = X[0].size();
  
  if(n <= thresh) {
    conquer(X, D, ids);
    return;
  }
  
  // Split Ids into sets
  std::vector<int> ids0(ids.begin(), ids.begin() + n/2);
  std::vector<int> ids1(ids.begin() + n/2, ids.end());
  
  std::vector<int> zo0;
  std::vector<int> zo1;
  std::vector<int> ids00;
  std::vector<int> ids01;
  std::vector<int> ids10;
  std::vector<int> ids11;
  for(int i : ids0){
    if(X[i][1] <= n/2){
      zo0.push_back(0);
      ids00.push_back(i);
    }else{
      zo0.push_back(1);
      ids01.push_back(i);
    }
  }
  for(int i : ids1){
    if(X[i][1] <= n/2){
      zo1.push_back(0);
      ids10.push_back(i);
    }else{
      zo1.push_back(1);
      ids11.push_back(i);
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
  divideAndConquer1(X, D, thresh, ids00, ids10, 1, true);
  divideAndConquer1(X, D, thresh, ids01, ids11, 1, false);
  
  // can move on to next dimension
  divideAndConquer1(X, D, thresh, ids00, ids11, 2, true);
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
  Rcpp::IntegerMatrix C_matrix(n, d - 1);
  for (int i = 0; i < n; i++) {
    C_matrix(i, 0) = D[i][0];
    for (int j = 1; j < d - 1; j++) {
      C_matrix(i, j) = C_matrix(i, j-1) + D[i][j];
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d - 1; j++) {
      C_matrix(i, j) = n - 1 - C_matrix(i, j);
    }
  }
  
  return C_matrix;
}

