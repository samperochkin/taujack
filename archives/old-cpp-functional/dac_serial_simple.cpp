#include <Rcpp.h>
#include <vector>

using namespace Rcpp;


//* 
//* Functions to "naively" compute concordance
//* 

void conquer(const std::vector<std::vector<double>>& X,
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

void conquer2(const std::vector<std::vector<double>>& X,
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

void merge(const std::vector<std::vector<double>>& X,
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


std::pair<int, int> findSplit(const std::vector<std::vector<double>>& X,
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
  int num_side_to_side = n0*n1; // to minimize
  int temp_nss;
  
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


std::pair<int, int> findSplitSort(const std::vector<std::vector<double>>& X,
                              std::vector<int>& ids,
                              const std::vector<int>& ids_copy,
                              const int l0, const int r0,
                              const int l1, const int r1,
                              const int k0) {

  int n0 = r0-l0;
  int n1 = r1-l1;
  int count0 = 0;
  int count1 = 0;
  int split0;
  int split1;
  int num_side_to_side = n0*n1; // to minimize
  int temp_nss;

  std::vector<int> ids_sorted(n0+n1);
  
  while (count0 < n0 && count1 < n1) {
    if (X[ids_copy[count0]][k0] < X[ids_copy[n0+count1]][k0]) {
      ids_sorted[count0+count1] = ids_copy[count0];
      count0 ++;
    } else {
      ids_sorted[count0+count1] = ids_copy[n0+count1];
      count1 ++;
    }
    temp_nss = count0 * count1 + (n0 - count0) * (n1 - count1);
    if (temp_nss <= num_side_to_side) {
      num_side_to_side = temp_nss;
      split0 = count0;
      split1 = count1;
    }
  }

  for(int i=count0; i<n0; i++) ids_sorted[n1+i] = ids[l0+i];
  for(int i=count1; i<n1; i++) ids_sorted[n0+i] = ids[l1+i];
  for(int i=0; i<n0; i++) ids[l0+i] = ids_sorted[i];
  for(int i=0; i<n1; i++) ids[l1+i] = ids_sorted[n0+i];
  
  return std::make_pair(split0, split1);
}

//* 
//* Main functions
//* 

void divideAndConquer2(const std::vector<std::vector<double>>& X,
                       std::vector<std::vector<int>>& D,
                       const int thresh,
                       std::vector<int>& ids,
                       const int l0, const int r0,
                       const int l1, const int r1,
                       const int k0,
                       const bool need_sort) {
  
  int n0 = r0-l0;
  int n1 = r1-l1;
  int d = X[0].size();
  
  // k0 = d means that the given pairs are concordant
  if(k0 == d) return;
  
  double sqrt_nn = std::sqrt(n0)*std::sqrt(n1);
  if(sqrt_nn <= thresh || n0 < 3 || n1 < 3) {
    conquer2(X, D, ids, l0, r0, l1, r1, k0);
    return;
  }

  // if(k0+1 == d){
  //   merge(X, D, ids, l0, r0, l1, r1, k0);
  //   return;
  // }
  
  // When skipping to a new dimension, need to ensure proper ordering
  if(need_sort){
    std::sort(ids.begin() + l0, ids.begin() + r0, [&](int i, int j) {
      return X[i][k0] < X[j][k0];
    });
    std::sort(ids.begin() + l1, ids.begin() + r1, [&](int i, int j) {
      return X[i][k0] < X[j][k0];
    });
  }  
  
  // find split value
  std::pair<int, int> splitIndices = findSplit(X, ids, l0, r0, l1, r1, k0);
  int split0 = splitIndices.first;
  int split1 = splitIndices.second;

  // rec necessarily discordant
  for(int i = l0 + split0; i < r0; i++) D[ids[i]][k0-1] += split1;
  for(int j = l1; j < l1 + split1; j++) D[ids[j]][k0-1] += (n0-split0);
  
  // need further investigation on k0+1 (re-divide)
  divideAndConquer2(X, D, thresh, ids, l0, l0+split0, l1, l1+split1, k0, false);
  divideAndConquer2(X, D, thresh, ids, l0+split0, r0, l1+split1, r1, k0, false);
  
  // can move on to next dimension (need_sort = true will break ordering)
  if(k0+1 < d){
    std::vector<int> ids_copy(ids.begin() + l0, ids.begin() + r1);
    divideAndConquer2(X, D, thresh, ids_copy, 0, split0, l1-l0+split1, r1-l0, k0+1, true);
  }
}

void divideAndConquer(const std::vector<std::vector<double>>& X,
                      std::vector<std::vector<int>>& D,
                      const int thresh,
                      std::vector<int>& ids,
                      const int l, const int r,
                      const int k0) {
  
  int n = r - l;
  int d = X[0].size();
  
  if(n <= thresh) {
    conquer(X, D, ids, l, r, k0);
    if(k0 < d){
      // end of recursion, order indices for merging
      std::sort(ids.begin() + l, ids.begin() + r, [&](int i, int j) {
        return X[i][k0] < X[j][k0];
      });
    }
    return;
  }
  
  
  // divide -----------------------------------------------
  int mid = l + n/2;
  divideAndConquer(X, D, thresh, ids, l, mid, k0);
  divideAndConquer(X, D, thresh, ids, mid, r, k0);

  
  // merge/conquer ----------------------------------------
  if(k0+1 == d){
    merge(X, D, ids, l, mid, mid, r, k0);
    return;
  }
  
  std::vector<int> ids_copy(ids.begin()+l, ids.begin()+r);
  std::pair<int, int> splitIndices = findSplitSort(X, ids, ids_copy, l, mid, mid, r, k0);
  int split0 = splitIndices.first;
  int split1 = splitIndices.second;

  // new limits
  int n0 = mid-l;
  int n1 = r-mid;
  
  for(int i = split0; i < n0; i++) D[ids_copy[i]][k0-1] += split1;
  for(int j = 0; j < split1; j++) D[ids_copy[n0+j]][k0-1] += n0 - split0;
  

  // need further investigation on k0+1 (re-divide)
  divideAndConquer2(X, D, thresh, ids_copy, 0, split0, n0, n0+split1, k0, false);
  divideAndConquer2(X, D, thresh, ids_copy, split0, n0, n0+split1, n0+n1, k0, false);

  // can move on to next dimension
  if(k0+1 < d){
    divideAndConquer2(X, D, thresh, ids_copy, 0, split0, n0+split1, n0+n1, k0+1, true);
  }
}


// [[Rcpp::export]]
IntegerMatrix dac_sserial(NumericMatrix X, int thresh = 100, bool brute_force = false) {
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
    conquer(X_vec, D, ids, 0, n, 1);
  }else{
    divideAndConquer(X_vec, D, thresh, ids, 0, n, 1);
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

