#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

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
  
  // // Sort columns of the rankMatrix based on the first column
  // std::vector<int> sortOrder(n);
  // std::iota(sortOrder.begin(), sortOrder.end(), 0);
  // std::sort(sortOrder.begin(), sortOrder.end(), [&](int i, int j) {
  //   return rankMatrix[i][0] < rankMatrix[j][0];
  // });
  // 
  // // Permute the rankMatrix based on the sortOrder
  // std::vector<std::vector<int>> permutedRankMatrix(n, std::vector<int>(d));
  // for (int j = 0; j < d; ++j) {
  //   for (int i = 0; i < n; ++i) {
  //     permutedRankMatrix[i][j] = rankMatrix[sortOrder[i]][j];
  //   }
  // }  
  // return permutedRankMatrix;
}


std::vector<int> cheatSort(const std::vector<std::vector<int>>& P,
                           const std::vector<int>& ids,
                           int l, int r,
                           int k_to) {
  
  int n = ids.size();
  std::vector<int> new_ids(r-l);
  std::vector<int> dummy(n,-1);
  for(int i=l; i<r; i++) dummy[P[ids[i]][k_to]] = ids[i];
  
  int count = 0;
  int hit = 0;
  while(count < n && hit < r-l){
    if(dummy[count] > -1) new_ids[hit++] = dummy[count];
    count++;
  }
  
  return new_ids;
}



// Function to reorder matrix X based on the ranks of its columns
// [[Rcpp::export]]
IntegerVector testCS(NumericMatrix X, int l, int r, int k0) {
  int n = X.nrow();
  int d = X.ncol();
  
  // Initialize X_vec
  std::vector<std::vector<double>> X_vec(n, std::vector<double>(d));
  
  // Fill X_vec
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
  
  std::vector<std::vector<int>> P = columnRanks(X_vec);
  // IntegerMatrix P_mat(n, d);
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < d; j++) {
  //     P_mat(i,j) = P[i][j];
  //   }
  // }
  // return P_mat;
  
  // Compute the reordered matrix
  std::vector<int> new_ids = cheatSort(P, ids, l, r, k0); // Sort based on column i
  IntegerVector new_ids_vec(r-l, d);
  for (int i = 0; i < r-l; i++) new_ids_vec(i) = new_ids[i];

  return new_ids_vec;
}
