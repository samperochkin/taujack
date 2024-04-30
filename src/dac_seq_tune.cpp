#include <Rcpp.h>
#include <vector>
#include <algorithm>  // For sorting
// #include <limits>
// #include <thread>  // for std::this_thread::sleep_for

using namespace Rcpp;


//* 
//* Functions to "naively" compute concordance
//* 

void conquer(const std::vector<std::vector<int>>& X,
             std::vector<std::vector<int>>& D,
             const std::vector<int>& ids,
             const int l, const int r,
             int k0) {
  
  int d = X[0].size();
  for (int i=l; i<r-1; i++) {
    for (int j=i+1; j<r; j++) {
      for (int k = k0; k < d; k++){
        if(X[ids[i]][k] > X[ids[j]][k]){
          D[ids[i]][k-1] += 1;
          D[ids[j]][k-1] += 1;
          break;
        }
      }  
    }
  }
}

void conquer2(const std::vector<std::vector<int>>& X,
              std::vector<std::vector<int>>& D,
              const std::vector<int>& ids,
              const int l0, const int r0,
              const int l1, const int r1,
              int k0) {
  int d = X[0].size();
  
  // Make sure X[][0] is sorted
  for (int i=l0; i<r0; i++) {
    for (int j=l1; j<r1; j++) {
      for (int k = k0; k < d; k++) {
        if (X[ids[i]][k] > X[ids[j]][k]){
          D[ids[i]][k-1] += 1;
          D[ids[j]][k-1] += 1;
          break;
        }
      }
    }
  }
}

void merge(const std::vector<std::vector<int>>& X,
           std::vector<std::vector<int>>& D,
           std::vector<int>& ids,
           const int l, const int mid, const int r,
           const int k0) {
  
  int n0 = mid-l;
  int n1 = r-mid;
  int count0 = 0;
  int count1 = 0;
  std::vector<int> temp_ids(n0+n1);
  
  while (count0 < n0 && count1 < n1) {
    if (X[ids[l+count0]][k0] < X[ids[mid+count1]][k0]) {
      D[ids[l+count0]][k0-1] += count1;
      temp_ids[count0+count1] = ids[l+count0];
      count0++;
    } else {
      D[ids[mid+count1]][k0-1] += n0-count0;
      temp_ids[count0+count1] = ids[mid+count1];
      count1++;
    }
  }
  for(int i=count0; i<n0; i++){
    D[ids[l+i]][k0-1] += n1; 
    temp_ids[n1+i] = ids[l+i];
  }
  for(int i=count1; i<n1; i++) temp_ids[n0+i] = ids[mid+i];
  for(int i=0; i<n0+n1; i++) ids[l+i] = temp_ids[i];
}


void merge2(const std::vector<std::vector<int>>& X,
           std::vector<std::vector<int>>& D,
           const std::vector<int>& ids,
           const int l0, const int r0,
           const int l1, const int r1,
           const int k0) {
  
  int n0 = r0-l0;
  int n1 = r1-l1;
  int count0 = 0;
  int count1 = 0;
  
  while (count0 < n0 && count1 < n1) {
    if (X[ids[l0+count0]][k0] < X[ids[l1+count1]][k0]) {
      D[ids[l0+count0]][k0-1] += count1;
      count0++;
    } else {
      D[ids[l1+count1]][k0-1] += n0-count0;
      count1++;
    }
  }
  for(int i=count0; i<n0; i++) D[ids[l0+i]][k0-1] += n1;
}


std::pair<int, int> findSplit(const std::vector<std::vector<int>>& X,
                              const std::vector<int>& ids,
                              const int l0, const int r0,
                              const int l1, const int r1,
                              const int k0) {
  
  int n0 = r0-l0;
  int n1 = r1-l1;
  int count0 = 0;
  int count1 = 0;
  int split0;
  int split1;
  long long num_side_to_side = n0;
  num_side_to_side = num_side_to_side * n1;
  long long temp_nss;
  
  
  while (count0 < n0 && count1 < n1) {
    if (X[ids[l0+count0]][k0] < X[ids[l1+count1]][k0]) {
      count0 ++;
    } else {
      count1 ++;
    }
    temp_nss = count0 * count1 + (n0 - count0) * (n1 - count1);
    if (temp_nss <= num_side_to_side) {
      num_side_to_side = temp_nss;
      split0 = count0;
      split1 = count1;
    }
  }
  
  return std::make_pair(split0, split1);
}


std::pair<int, int> findSplitSort(const std::vector<std::vector<int>>& X,
                                  const std::vector<int>& ids,
                                  std::vector<int>& idss,
                                  const int l0, const int r0,
                                  const int l1, const int r1,
                                  const int k0) {
  
  int n0 = r0-l0;
  int n1 = r1-l1;
  int count0 = 0;
  int count1 = 0;
  int split0;
  int split1;
  long long num_side_to_side = n0;
  num_side_to_side = num_side_to_side * n1;
  long long temp_nss;
  
  while (count0 < n0 && count1 < n1) {
    if (X[ids[l0+count0]][k0] < X[ids[l1+count1]][k0]) {
      idss[l0+count0+count1] = ids[l0+count0];
      count0 ++;
    } else {
      idss[l0+count0+count1] = ids[l1+count1];
      count1 ++;
    }
    
    temp_nss = count0 * count1 + (n0 - count0) * (n1 - count1);
    if (temp_nss <= num_side_to_side) {
      num_side_to_side = temp_nss;
      split0 = count0;
      split1 = count1;
    }
  }
  
  for(int i=count0; i<n0; i++) idss[l0+n1+i] = ids[l0+i];
  for(int i=count1; i<n1; i++) idss[l0+n0+i] = ids[l1+i];

  return std::make_pair(split0, split1);
}



//* 
//* Main functions
//* 

void divideAndConquer2(const std::vector<std::vector<int>>& X,
                       std::vector<std::vector<int>>& D,
                       const int thresh,
                       std::vector<int>& ids,
                       const int l0, const int r0,
                       const int l1, const int r1,
                       const int k0,
                       const bool need_sort,
                       const int tune) {
  
  int n0 = r0-l0;
  int n1 = r1-l1;
  if(n0 == 0 | n1 == 0) return;
  int d = X[0].size();
  if(k0 == d) return; // end of recursion - concordant sets
  
  double sqrt_nn = std::sqrt(n0)*std::sqrt(n1);
  if(sqrt_nn <= thresh || n0 < tune || n1 < tune) {
    conquer2(X, D, ids, l0, r0, l1, r1, k0);
    return;
  }
  
  // When skipping to a new dimension, need to ensure proper ordering
  if(need_sort){
    
    bool okay0 = true;
    bool okay1 = true;

    for(int i=l0; i<r0-1; i++){
      if(X[ids[i]][k0] > X[ids[i+1]][k0]){
        okay0 = false;
        break;
      }
    }
    for(int i=l1; i<r1-1; i++){
      if(X[ids[i]][k0] > X[ids[i+1]][k0]){
        okay1 = false;
        break;
      }
    }

    if(!okay0){
      std::sort(ids.begin() + l0, ids.begin() + r0, [&](int i, int j) {
        return X[i][k0] < X[j][k0];
      });
    }
    if(!okay1){
      std::sort(ids.begin() + l1, ids.begin() + r1, [&](int i, int j) {
        return X[i][k0] < X[j][k0];
      });
    }
  }  
  
  // if appropriate perform merge phase of merge sort
  if(k0+1 == d){
    merge2(X, D, ids, l0, r0, l1, r1, k0);
    return;
  }

  // find split value
  std::pair<int, int> splitIndices = findSplit(X, ids, l0, r0, l1, r1, k0);
  int split0 = splitIndices.first;
  int split1 = splitIndices.second;
  
  // rec necessarily discordant
  for(int i = l0 + split0; i < r0; i++) D[ids[i]][k0-1] += split1;
  for(int j = l1; j < l1 + split1; j++) D[ids[j]][k0-1] += (n0-split0);
  
  // need further investigation on k0+1 (re-divide)
  divideAndConquer2(X, D, thresh, ids, l0, l0+split0, l1, l1+split1, k0, false, tune);
  divideAndConquer2(X, D, thresh, ids, l0+split0, r0, l1+split1, r1, k0, false, tune);
  
  // can move on to next dimension (need_sort = true will break ordering)
  if(k0+1 < d){
    divideAndConquer2(X, D, thresh, ids, l0, l0+split0, l1+split1, r1, k0+1, true, tune);
  }
}

void divideAndConquer(const std::vector<std::vector<int>>& X,
                      std::vector<std::vector<int>>& D,
                      const int thresh,
                      std::vector<int>& ids,
                      std::vector<int>& idss,
                      const int l, const int r,
                      const int k0,
                      const int tune) {
  
  int n = r - l;
  int d = X[0].size();
  
  if(n <= thresh) {
    conquer(X, D, ids, l, r, k0);
    // end of recursion, order indices for later merging
    std::sort(ids.begin() + l, ids.begin() + r, [&](int i, int j) {
      return X[i][k0] < X[j][k0];
    });
    return;
  }
  
  // divide -----------------------------------------------
  int mid = l + n/2;
  divideAndConquer(X, D, thresh, ids, idss, l, mid, k0, tune);
  divideAndConquer(X, D, thresh, ids, idss, mid, r, k0, tune);
  
  // merge/conquer ----------------------------------------
  if(k0+1 == d){
    merge(X, D, ids, l, mid, r, k0);
    return;
  }

  // std::vector<int> ids_sorted(ids.begin()+l, ids.begin()+r);
  std::pair<int, int> splitIndices = findSplitSort(X, ids, idss, l, mid, mid, r, k0);
  int split0 = splitIndices.first;
  int split1 = splitIndices.second;
  
  for(int i = split0; i < mid-l; i++) D[ids[l+i]][k0-1] += split1;
  for(int j = 0; j < split1; j++) D[ids[mid+j]][k0-1] += mid - l - split0;
  
  // need further investigation on k0 (re-divide)
  divideAndConquer2(X, D, thresh, ids, l, l+split0, mid, mid+split1, k0, false, tune);
  divideAndConquer2(X, D, thresh, ids, l+split0, mid, mid+split1, r, k0, false, tune);
  
  // can move on to next dimension
  if(k0+1 < d){
    divideAndConquer2(X, D, thresh, ids, l, l+split0, mid+split1, r, k0+1, true, tune);
  }
  
  // enforce ordering on k0
  for(int i=0; i<r-l; i++) ids[l+i] = idss[l+i];
}


std::vector<int> computeRanks(const std::vector<double>& values) {
  std::vector<int> ranks(values.size());
  std::iota(ranks.begin(), ranks.end(), 0);
  
  // Sort the indices based on the values
  std::sort(ranks.begin(), ranks.end(), [&](int i, int j) {
    return values[i] < values[j];
  });
  
  // Create the rank vector
  std::vector<int> rankVector(values.size());
  for (int i = 0; i < values.size(); ++i) {
    rankVector[ranks[i]] = i;  // 0-based
  }
  
  return rankVector;
}

// Function to compute column ranks of a matrix
std::vector<std::vector<int>> columnRanks(const std::vector<std::vector<double>>& X) {
  int n = X.size();
  int d = X[0].size();
  std::vector<std::vector<int>> rankMatrix(n, std::vector<int>(d));
  
  // Compute ranks for each column
  for (int j = 0; j < d; ++j) {
    std::vector<double> columnValues(n);
    for (int i = 0; i < n; ++i) columnValues[i] = X[i][j];
    std::vector<int> columnRanks = computeRanks(columnValues);
    for (int i = 0; i < n; ++i) rankMatrix[i][j] = columnRanks[i];
  }
  return rankMatrix;
}

// Rcpp wrapper function
// This function takes an R numeric matrix X, computes the discordance matrix, and returns an integer matrix.
// [[Rcpp::export]]
IntegerMatrix dac_seq_tune(NumericMatrix X, int thresh = 25, bool brute_force = false, int tune = 3) {
  int n = X.nrow();
  int d = X.ncol();
  
  // initialize D
  std::vector<std::vector<int>> D(n, std::vector<int>(d - 1, 0));
  
  // initialize X_vec
  std::vector<std::vector<double>> X_vec(n, std::vector<double>(d));
  // fill X_vec
  for (int i = 0; i < n; i++)
    for (int j = 0; j < d; j++) X_vec[i][j] = X(i, j);
  std::vector<std::vector<int>> U = columnRanks(X_vec);
  
  std::vector<int> ids(n);
  for(int i=0; i<n; i++) ids[U[i][0]] = i;
  std::vector<int> idss = ids;
  
  if(brute_force){
    conquer(U, D, ids, 0, n, 1);
  }else{
    divideAndConquer(U, D, thresh, ids, idss, 0, n, 1, tune);
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
