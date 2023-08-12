#include <Rcpp.h>
#include <vector>

using namespace Rcpp;


//* 
//* Functions to "naively" compute concordance
//* 

void sortI(const std::vector<std::vector<double>>& X,
           std::vector<std::vector<int>>& I,
           const std::vector<int>& ids,
           const int k0) {
  for (size_t j = k0; j < I[0].size(); ++j) {
    
    std::vector<std::pair<int, int>> col_data;
    for (int i : ids) {
      col_data.emplace_back(I[i][j], i);
    }
    
    // Sort col_data based on values
    std::sort(col_data.begin(), col_data.end(), [&](int i, int j) {
      return X[i][k0+1] < X[j][k0+1];
    });
    
    // Update the matrix with sorted values
    for (size_t i = 0; i < ids.size(); ++i) {
      I[ids[i]][j] = col_data[i].first;
    }
  }
}

void conquer(const std::vector<std::vector<double>>& X,
             std::vector<std::vector<int>>& D,
             const std::vector<int>& ids,
             int k0) {
  int n = ids.size();
  int d = X[0].size();
  
  if(k0 == 0) k0++;
  
  // Make sure X[,1] is sorted
  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
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

// HERE FIND A WAY TO UPDATE ids ids (if necessary)
std::pair<int, int> findSplit(const std::vector<std::vector<double>>& X,
                              std::vector<std::vector<int>>& I,
                              const std::vector<int>& ids0,
                              const std::vector<int>& ids1,
                              const int k0,
                              const bool merge) {
  int n0 = ids0.size();
  int n1 = ids1.size();
  int count0 = 0;
  int count1 = 0;
  int split0;
  int split1;
  int num_side_to_side = 0; // to minimize
  int nss_temp;

  if(merge){
    int temp_ids(n0+n1);
    while (count0 < n0 && count1 < n1) {
      if (X[ids0[count0]][k0] < X[ids1[count1]][k0]) {
        temp_ids[count0 + count1] = ids0[count0];
        count0++;
      } else {
        temp_ids[count0 + count1] = ids1[count1];
        count1++;
      }
      nss_temp = count0 * count1 + (n0 - count0) * (n1 - count1);
      if (nss_temp <= num_side_to_side) {
        num_side_to_side = nss_temp;
        split0 = count0;
        split1 = count1;
      }
    }
    
    for(int j=0; j<n0+n1; j++) ids[j] = temp_ids[j];
    
  }else{
    while (count0 < n0 && count1 < n1) {
      if (X[ids0[count0]][k0] < X[ids1[count1]][k0]) {
        count0 += 1;
      } else {
        count1 += 1;
      }
      nss_temp = count0 * count1 + (n0 - count0) * (n1 - count1);
      if (nss_temp <= num_side_to_side) {
        num_side_to_side = nss_temp;
        split0 = count0;
        split1 = count1;
      }
    }
  }
  
  return std::make_pair(split0, split1);
}


//* 
//* Main functions
//* 

void divideAndConquer(const std::vector<std::vector<double>>& X,
                      std::vector<std::vector<int>>& I,
                      std::vector<std::vector<int>>& D,
                      const int thresh,
                      std::vector<int>& ids,
                      const int k0) {
  
  int n = ids.size();
  int d = X[0].size();
  
  
  // base case --------------------------------------------
  if(n <= thresh) {
    conquer(X, D, ids, k0);
    sortI(X, I, ids, k0+1);
    return;
  }
  
  
  // divide -----------------------------------------------
  std::vector<int> ids0(ids.begin(), ids.begin() + n/2);
  std::vector<int> ids1(ids.begin() + n/2, ids.end());
  divideAndConquer(X, I, D, thresh, ids0, k0);
  divideAndConquer(X, I, D, thresh, ids1, k0);
  // At this point ids0 and ids1 are ordered wrt to column k0+1
      
  // merge/conquer ----------------------------------------
  // find split value
  std::pair<int, int> splitIndices = findSplit(X, I, ids0, ids1, k0+1, false);
  int split0 = splitIndices.first;
  int split1 = splitIndices.second;
  
  // new vectors (to keep ids0 and ids1 intact for later)
  std::vector<int> ids00(ids0.begin(), ids0.begin() + split0);
  std::vector<int> ids01(ids0.begin() + split0, ids0.end());
  std::vector<int> ids10(ids1.begin(), ids1.begin() + split1);
  std::vector<int> ids11(ids1.begin() + split1, ids1.end());
  // std::vector<int>& ids00 = ids0;
  // std::vector<int>& ids01 = ids0;
  // std::vector<int>& ids10 = ids1;
  // std::vector<int>& ids11 = ids1;
  // ids01.erase(ids01.begin(), ids01.begin() + split0);
  // ids10.erase(ids10.begin() + split1, ids10.end());
  // ids11.erase(ids11.begin(), ids11.begin() + split1);
  
  // rec necessarily discordant (index k0 is for dims (k0, k0+1))
  for(int i : ids01) D[i][k0] += split1;
  for(int j : ids10) D[j][k0] += n/2 - split0;
  
  // need further investigation on k0+1 (re-divide)
  divideAndConquer2(X, D, thresh, ids00, ids10, k0+1, false);
  divideAndConquer2(X, D, thresh, ids01, ids11, k0+1, false);
  
  // can move on to next dimension
  divideAndConquer1(X, D, thresh, ids00, ids11, k0+2, true);
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

