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
  
  if(k0 == 0) k0++;
    
  // Make sure X[,1] is sorted
  for (int i=l; i<r-1; i++) {
    for (int j=i+1; j<r; j++) {
      for (size_t k = k0; k < d; k++){
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
  
  if(k0 == 0) k0++;
  
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

void merge(const std::vector<std::vector<double>>& X,
           std::vector<int>& ids,
           const int l, const int mid, const int r,
           const int k0) {
  
  int n = r - l;
  std::vector<int> mergedIds(n);
  
  int i = l;        // Index for the first subarray
  int j = mid;      // Index for the second subarray
  int mergedIndex = 0;
  
  while (i < mid && j < r) {
    // Compare values in the subarrays and merge them in order
    if (X[ids[i]][k0] < X[ids[j]][k0]) {
      mergedIds[mergedIndex++] = ids[i++];
    } else {
      mergedIds[mergedIndex++] = ids[j++];
    }
  }
  
  // Copy remaining elements from the first subarray, if any
  while (i < mid) {
    mergedIds[mergedIndex++] = ids[i++];
  }
  
  // Copy remaining elements from the second subarray, if any
  while (j < r) {
    mergedIds[mergedIndex++] = ids[j++];
  }
  
  // Copy the merged ids back to the original ids vector
  for (int k = 0; k < n; k++) {
    ids[l + k] = mergedIds[k];
  }
}



//* 
//* Main functions
//* 

bool areIndicesValid(int l, int r, int n) {
  return l >= 0 && l <= n && r >= 0 && r <= n && l <= r;
}
void testIndicesValidity(int l, int r, int n) {
  if (l < 0 || r > n || l > r) {
    throw std::invalid_argument("Invalid indices: l and r must be within [0, n] and l <= r.");
  }
}
bool isSplitValid(int split0, int split1, int n0, int n1) {
  return split0 >= 0 && split0 < n0 && split1 >= 0 && split1 < n1;
}
void testSplitValidity(int split0, int split1, int n0, int n1) {
  if (!isSplitValid(split0, split1, n0, n1)) {
    throw std::invalid_argument("Invalid split values: split0 and split1 must be within [0, n0) and [0, n1), respectively.");
  }
}

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
    
    Rcpp::Rcout << " -------------------------------- " << std::endl;
    Rcpp::Rcout << " base case 2 (w/ k0=" << k0 << ") on " << std::endl;
    for(int i=l0; i<r0; i++)
      Rcpp::Rcout << X[ids[i]][k0] << " ";
    Rcpp::Rcout << std::endl;;
    for(int i=l1; i<r1; i++)
      Rcpp::Rcout << X[ids[i]][k0] << " ";
    Rcpp::Rcout << std::endl;;
    Rcpp::Rcout << " -------------- " << std::endl;;
    
    Rcpp::Rcout << "original ids: ";
    for(int j=0; j < ids.size(); j++) Rcpp::Rcout << ids[j] << " ";
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << std::endl;
    
    Rcpp::Rcout << "D in:  ";
    for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][0] << " ";
    Rcpp::Rcout << " -- ";
    if(d > 2) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][1] << " ";
    Rcpp::Rcout << " -- ";
    if(d > 3) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][2] << " ";
    Rcpp::Rcout << std::endl;
    
    // If n is small enough, perform brute force computation
    conquer2(X, D, ids, l0, r0, l1, r1, k0);
    
    Rcpp::Rcout << "D out: ";
    for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][0] << " ";
    Rcpp::Rcout << " -- ";
    if(d > 2) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][1] << " ";
    Rcpp::Rcout << " -- ";
    if(d > 3) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][2] << " ";
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << std::endl;
    
    return;
  }
  
  if(need_sort){
    Rcpp::Rcout << "sorts ids (after dive in dimension " << k0 << "): " << std::endl;
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
  
  int l00 = l0;
  int r00 = l0 + split0;
  int l01 = l0 + split0;
  int r01 = r0;
  int l10 = l1;
  int r10 = l1 + split1;
  int l11 = l1 + split1;
  int r11 = r1;
  
  // rec necessarily discordant
  for(int i = l01; i < r01; i++) D[ids[i]][k0-1] += split1;
  for(int j = l10; j < r10; j++) D[ids[j]][k0-1] += (n0-split0);
  
  // need further investigation on k0+1 (re-divide)
  divideAndConquer2(X, D, thresh, ids, l00, r00, l10, r10, k0, false);
  divideAndConquer2(X, D, thresh, ids, l01, r01, l11, r11, k0, false);

  // can move on to next dimension (need_sort = true will break ordering)
  if(k0+1 < d){
    std::vector<int> ids_copy(ids.begin() + l00, ids.begin() + r11);
    divideAndConquer2(X, D, thresh, ids_copy, 0, r00-l00, l11-l00, r11-l00, k0+1, true);
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

  // base case --------------------------------------------
  if(n <= thresh) {
    
    Rcpp::Rcout << " -------------------------------- " << std::endl;
    Rcpp::Rcout << " bases case (w/ k0=" << k0 << ") on " << std::endl;
    for(int i=l; i<r; i++)
      Rcpp::Rcout << X[ids[i]][k0] << " ";
    Rcpp::Rcout << std::endl;;
    Rcpp::Rcout << " -------------- " << std::endl;;
    
    
    Rcpp::Rcout << "original ids: ";
    for(int j=0; j < ids.size(); j++) Rcpp::Rcout << ids[j] << " ";
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << std::endl;
    
    Rcpp::Rcout << "D in:  ";
    for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][0] << " ";
    Rcpp::Rcout << " -- ";
    if(d > 2) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][1] << " ";
    Rcpp::Rcout << " -- ";
    if(d > 3) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][2] << " ";
    Rcpp::Rcout << std::endl;
    
    conquer(X, D, ids, l, r, k0);
    if(k0 < d){
      // end of recursion, order indices for next dimension
      std::sort(ids.begin() + l, ids.begin() + r, [&](int i, int j) {
        return X[i][k0] < X[j][k0];
      });
    }
    
    Rcpp::Rcout << "D out: ";
    for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][0] << " ";
    Rcpp::Rcout << " -- ";
    if(d > 2) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][1] << " ";
    Rcpp::Rcout << " -- ";
    if(d > 3) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][2] << " ";
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << std::endl;
    
    Rcpp::Rcout << "new ids: ";
    for(int j=0; j < ids.size(); j++) Rcpp::Rcout << ids[j] << " ";
    Rcpp::Rcout << std::endl;
    
    return;
  }
  
  
  // divide -----------------------------------------------
  int mid = l + n/2;

  divideAndConquer(X, D, thresh, ids, l, mid, k0);
  divideAndConquer(X, D, thresh, ids, mid, r, k0);
  
  // merge/conquer ----------------------------------------
  // find split value
  std::pair<int, int> splitIndices = findSplit(X, ids, l, l + n/2, l + n/2, r, k0);
  int split0 = splitIndices.first;
  int split1 = splitIndices.second;
  
  // new limits
  int mid0 = l + split0;
  int mid1 = mid + split1;

  // rec necessarily discordant (index k0 is for dims (k0, k0+1))
  Rcpp::Rcout << " -------------------------------- " << std::endl;
  Rcpp::Rcout << " merge update (w/ k0=" << k0 << ") on " << std::endl;
  for(int i=l; i<r; i++)
    Rcpp::Rcout << X[ids[i]][k0] << " ";
  Rcpp::Rcout << std::endl;;
  Rcpp::Rcout << " -------------- " << std::endl;;
  
  Rcpp::Rcout << "original ids: ";
  for(int j=0; j < ids.size(); j++) Rcpp::Rcout << ids[j] << " ";
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << std::endl;

  Rcpp::Rcout << "split: ";
  Rcpp::Rcout << std::endl;
  for(int j=l; j < mid0; j++) Rcpp::Rcout << ids[j] << " ";
  Rcpp::Rcout << "| ";
  for(int j=mid0; j < mid; j++) Rcpp::Rcout << ids[j] << " ";
  Rcpp::Rcout << std::endl;
  for(int j=mid; j < mid1; j++) Rcpp::Rcout << ids[j] << " ";
  Rcpp::Rcout << "| ";
  for(int j=mid1; j < r; j++) Rcpp::Rcout << ids[j] << " ";
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << std::endl;

  Rcpp::Rcout << "D in:  ";
  for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][0] << " ";
  Rcpp::Rcout << " -- ";
  if(d > 2) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][1] << " ";
  Rcpp::Rcout << " -- ";
  if(d > 3) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][2] << " ";
  Rcpp::Rcout << std::endl;
  
  for(int i = mid0; i < mid; i++) D[ids[i]][k0-1] += split1;
  for(int j = mid; j < mid1; j++) D[ids[j]][k0-1] += (mid-l) - split0;

  Rcpp::Rcout << "D out: ";
  for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][0] << " ";
  Rcpp::Rcout << " -- ";
  if(d > 2) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][1] << " ";
  Rcpp::Rcout << " -- ";
  if(d > 3) for(int j=0; j < D.size(); j++) Rcpp::Rcout << D[j][2] << " ";
  Rcpp::Rcout << std::endl;
  Rcpp::Rcout << std::endl;
  
  Rcpp::Rcout << "new ids: ";
  for(int j=0; j < ids.size(); j++) Rcpp::Rcout << ids[j] << " ";
  Rcpp::Rcout << std::endl;
  
  // need further investigation on k0+1 (re-divide)
  divideAndConquer2(X, D, thresh, ids, l, mid0, mid, mid1, k0, false);
  divideAndConquer2(X, D, thresh, ids, mid0, mid, mid1, r, k0, false);
  
  // can move on to next dimension
  if(k0+1 < d){
    std::vector<int> ids_copy(ids.begin() + l, ids.begin() + r);
    divideAndConquer2(X, D, thresh, ids_copy, 0, mid0-l, mid1-l, r-l, k0+1, true);
  }
  
  Rcpp::Rcout << "merge operation (sort the two sorted vectors)" << std::endl;
  merge(X, ids, l, mid, r, k0);
  Rcpp::Rcout << "new ids: ";
  for(int j=0; j < ids.size(); j++) Rcpp::Rcout << ids[j] << " ";
  Rcpp::Rcout << std::endl;
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

