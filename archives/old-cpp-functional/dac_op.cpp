#include <Rcpp.h>
using namespace Rcpp;

// Function to compute concordance using Rcpp sugar functions
void conquer(const IntegerMatrix& X, IntegerVector& C, const IntegerVector& ids) {
  int n = ids.size();
  
  for (int i = 0; i < n - 1; i++) {
    for (int j = i + 1; j < n; j++) {
      LogicalVector conc = X(ids[i], _) <= X(ids[j], _);
      if (is_true(all(conc))) {
        C[ids[i]] += 1;
        C[ids[j]] += 1;
      }
    }
  }
}


// Function to compute concordance for two index sets using specific columns
void conquer2(const IntegerMatrix& X, IntegerVector& C, const IntegerVector& ids1, const IntegerVector& ids2, const IntegerVector& ks) {
  int n1 = ids1.size();
  int n2 = ids2.size();
  
  for (int i : ids1) {
    for (int j : ids2) {
      bool conc = true;
      for (int k = 0; k < ks.size(); ++k) {
        if (X(i, ks[k]) > X(j, ks[k])) {
          conc = false;
          break;
        }
      }
      if (conc) {
        C[i] += 1;
        C[j] += 1;
      }
    }
  }
}


// Function performing the merge phase of mergeSort for cases when there remains one dimension to look at
void merge(const IntegerMatrix& X, IntegerVector& C, IntegerVector& ids1, IntegerVector& ids2, const int k) {
  int n1 = ids1.size();
  int n2 = ids2.size();
  
  std::sort(ids1.begin(), ids1.end(), [&X, k](int i, int j){return X(i,k) < X(i,k);});
  std::sort(ids2.begin(), ids2.end(), [&X, k](int i, int j){return X(i,k) < X(i,k);});
  
  int i = 0, j = 0;
  while (i < n1 && j < n2) {
    if (X(ids2[j], k) < X(ids1[i], k)) {
      C[ids2[j]] += i;
      j++;
    } else {
      C[ids1[i]] += n2 - j;
      i++;
    }
  }
  while (j < n2) {
    C[ids2[j]] += n1;
    j++;
  }
}

// Divide and conquer function for pairs of index sets
void divideAndConquerLower(const IntegerMatrix& X, IntegerVector& C, const int thresh, IntegerVector& ids1, IntegerVector& ids2, IntegerVector& ks) {
  int n1 = ids1.size();
  int n2 = ids2.size();
  int d = ks.size();
  
  double sqrt_nn = sqrt(n1) * sqrt(n2);
  
  if (sqrt_nn <= thresh || n1 < 3 || n2 < 3) {
    conquer2(X, C, ids1, ids2, ks);
    return;
  }
  
  NumericVector mid(d);
  
  for (int k = 0; k < d; k++) {
    IntegerVector X_col = X(_, ks[k]);
    mid[k] = (max(X_col) + min(X_col)) / 2.0;
  }
  
  IntegerVector pos_vec(d, 0);
  std::map<IntegerVector, IntegerVector> pos1_map;
  std::map<IntegerVector, IntegerVector> pos2_map;
  
  for (int i : ids1) {
    IntegerVector pos(d);
    for (int k = 0; k < d; k++) {
      if (X(i, ks[k]) > mid[k]) {
        pos[k] = 1;
      }
    }
    pos1_map[pos].push_back(i);
  }
  
  for (int j : ids2) {
    IntegerVector pos(d);
    for (int k = 0; k < d; k++) {
      if (X(j, ks[k]) > mid[k]) {
        pos[k] = 1;
      }
    }
    pos2_map[pos].push_back(j);
  }
  
  for (const auto& entry : pos1_map) {
    IntegerVector ids1_subset = entry.second;
    IntegerVector ids2_subset = pos2_map[entry.first];
    
    if (ids1_subset.size() > 0 && ids2_subset.size() > 0) {
      divideAndConquerLower(X, C, thresh, ids1_subset, ids2_subset, ks);
    }
  }
  
  for (const auto& entry : pos2_map) {
    IntegerVector ids2_subset = entry.second;
    IntegerVector ids1_subset = pos1_map[entry.first];
    
    if (ids1_subset.size() > 0 && ids2_subset.size() > 0) {
      divideAndConquerLower(X, C, thresh, ids1_subset, ids2_subset, ks);
    }
  }
  
  merge(X, C, ids1, ids2, ks[d - 1]);
}


void divideAndConquer(const IntegerMatrix& X, IntegerVector& C, const int thresh, const IntegerVector& ids) {
  int n = ids.size();
  int d = X.ncol();
  
  if (n <= thresh) {
    conquer(X, C, ids);
    return;
  }
  
  NumericVector mid(d);
  
  for (int k = 0; k < d; k++) {
    IntegerVector X_col = X(_, k);
    mid[k] = (max(X_col) + min(X_col)) / 2.0;
  }
  
  IntegerVector pos_vec(d, 0);
  std::map<IntegerVector, IntegerVector> pos_map;
  
  for (int i : ids) {
    IntegerVector pos(d);
    for (int k = 0; k < d; k++) {
      if (X(i, k) > mid[k]) {
        pos[k] = 1;
      }
    }
    pos_map[pos].push_back(i);
  }
  
  for (auto it1 = pos_map.begin(); it1 != pos_map.end(); ++it1) {
    auto pos1 = it1->first;
    auto ids1 = it1->second;
    int n1 = ids1.size();
    
    divideAndConquer(X, C, thresh, ids);
    
    for (auto it2 = std::next(it1); it2 != pos_map.end(); ++it2) {
      auto pos2 = it2->first;
      auto ids2 = it2->second;
      int n2 = ids2.size();
      
      int v_e = 0;
      int v_g = 0;
      for (int i = 0; i < d; i++) {
        v_e += (pos1[i] == pos2[i]);
        v_g += (pos1[i] > pos2[i]);
      }
      
      if ((v_g > 0) && (v_e + v_g < d)) { // blocks that are necessarily discordant
        continue;
      } else if (v_e == 0) { // blocks necessarily concordant
        for (int i : ids1) C[i] += n2;
        for (int j : ids2) C[j] += n1;
      } else if (v_e == 1) { // one-dimensional problem
        int k;
        for (int i = 0; i < pos1.size(); i++) {
          if (pos1[i] == pos2[i]) {
            k = i;
            break;
          }
        }
        merge(X, C, ids1, ids2, k);
      } else { // k-dimensional problem (k > 1)
        IntegerVector ks;
        for (int i = 0; i < d; i++) {
          if (pos1[i] == pos2[i]) {
            ks.push_back(i);
          }
        }
        divideAndConquerLower(X, C, thresh, ids1, ids2, ks);
      }
    }
  }
}


// Main function to compute concordance using divide-and-conquer approach with Rcpp sugar functions
// [[Rcpp::export]]
IntegerVector dacOP(IntegerMatrix X, int thresh = 100, bool brute_force = false) {
  int n = X.nrow();
  int d = X.ncol();
  
  // initialize C
  IntegerVector C(n, 0);
  IntegerVector ids(n);
  std::iota(ids.begin(), ids.end(), 0);

  if(brute_force){
    conquer(X, C, ids);
  }else{
    divideAndConquer(X, C, thresh, ids);
  }
  
  // convert C and return it
  return C;
}