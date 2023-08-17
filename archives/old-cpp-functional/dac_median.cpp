#include <Rcpp.h>
#include <vector>

using namespace Rcpp;


//* 
//* Functions to "naively" compute concordance
//* 

void conquer(const std::vector<std::vector<int>>& X,
             std::vector<int>& C,
             const std::vector<int>& ids) {
  int n = ids.size();
  int d = X[0].size();
  bool all_less, all_greater;
  
  for (int i=0; i<n-1; i++) {
    for (int j=i+1; j<n; j++) {
      all_less = true;
      all_greater = true;
      for (size_t k = 0; k < d; k++) {
        if (X[ids[i]][k] < X[ids[j]][k]) {
          all_greater = false;
        } else if (X[ids[i]][k] > X[ids[j]][k]) {
          all_less = false;
        }
        if (!all_less && !all_greater) {
          break;
        }
      }
      if (all_less || all_greater) {
        C[ids[i]] += 1;
        C[ids[j]] += 1;
      }
    }
  }
}

void conquer2(const std::vector<std::vector<int>>& X,
              std::vector<int>& C,
              const std::vector<int>& ids1,
              const std::vector<int>& ids2,
              std::vector<int>& ks) {
  int n1 = ids1.size();
  int n2 = ids2.size();
  int d = ks.size();
  bool all_less, all_greater;
  
  for (int i : ids1) {
    for (int j : ids2) {
      all_less = true;
      for (int k = 0; k < d; k++) {
        if (X[i][ks[k]] > X[j][ks[k]]) all_less = false;
      }
      if (all_less) {
        C[i] += 1;
        C[j] += 1;
      }
    }
  }
}



//* 
//* Function performing the merge phase of mergeSort
//* for cases when there remains one dimension to look at.
//* 

void merge(const std::vector<std::vector<int>>& X,
           std::vector<int>& C,
           std::vector<int>& ids1,
           std::vector<int>& ids2,
           const int k) {
  
  int n1 = ids1.size();
  int n2 = ids2.size();
  
  // Sort ids1 and ids2 by the values in X[:,k]
  std::sort(ids1.begin(), ids1.end(), [&X, k](int i, int j){return X[i][k] < X[j][k];});
  std::sort(ids2.begin(), ids2.end(), [&X, k](int i, int j){return X[i][k] < X[j][k];});
  
  // merge operation in merge sort
  int i = 0, j = 0;
  while (i < n1 && j < n2) {
    if (X[ids2[j]][k] < X[ids1[i]][k]) {
      C[ids2[j]] += i;
      j++;
    } else {
      C[ids1[i]] += n2-j;
      i++;
    }
  }
  while (j < n2) {
    C[ids2[j]] += n1;
    j++;
  }
}




//* 
//* Divide and conquer function for pairs of index sets.
//* 

void divideAndConquerLower(const std::vector<std::vector<int>>& X,
                           std::vector<int>& C,
                           const int thresh,
                           std::vector<int>& ids1,
                           std::vector<int>& ids2,
                           std::vector<int>& ks) {
                           // const bool median) {
  
  int n1 = ids1.size();
  int n2 = ids2.size();
  int d = ks.size();
  
  double sqrt_nn = std::sqrt(n1)*std::sqrt(n2);
  
  if(sqrt_nn <= thresh || n1 < 3 || n2 < 3) {
    // If n is small enough, perform brute force computation
    conquer2(X, C, ids1, ids2, ks);
    return;
  }
  
  // Otherwise, split Ids into 2^d sets relative to the midpoint
  std::vector<double> mid(d);
  // if(median){
  //   mid = midpoint2(X, ids1, ids2, ks);
  // }else{
  int mi;
  int ma;
  for (int k = 0; k < d; k++) {
    mi = X.size();
    ma = 1;
    for (int i : ids1) {
      mi = std::min(mi, X[i][ks[k]]);
      ma = std::max(ma, X[i][ks[k]]);
    }
    for (int i : ids2) {
      mi = std::min(mi, X[i][ks[k]]);
      ma = std::max(ma, X[i][ks[k]]);
    }
    mid[k] = (ma + mi)/2;
  }
  // }

  std::vector<int> pos_vec(d, 0);
  std::map<std::vector<int>, std::vector<int>> pos1_map;
  std::map<std::vector<int>, std::vector<int>> pos2_map;
  
  for (int i : ids1) {
    for (int k = 0; k < d; k++) pos_vec[k] = (X[i][ks[k]] >= mid[k]) ? 1 : 0;
    pos1_map[pos_vec].push_back(i); // add index i to the corresponding pos vector in the map
  }
  for (int i : ids2) {
    for (int k = 0; k < d; k++) pos_vec[k] = (X[i][ks[k]] >= mid[k]) ? 1 : 0;
    pos2_map[pos_vec].push_back(i); // add index i to the corresponding pos vector in the map
  }
  
  
  // Loop over all pairs of groups
  for (auto it1 = pos1_map.begin(); it1 != pos1_map.end(); ++it1) {
    
    // Access the key and value of the first group
    auto pos1i = it1->first;
    auto ids1i = it1->second;
    int n1i = ids1i.size();
    
    for(auto it2 = pos2_map.begin(); it2 != pos2_map.end(); ++it2) {
      
      // Access the key and value of the second group
      auto pos2j = it2->first;
      auto ids2j = it2->second;
      int n2j = ids2j.size();
      
      // record relation between pos1i and pos2j
      int v_e = 0;
      int v_g = 0;
      for(int i=0; i<d; i++){
        v_e += (pos1i[i] == pos2j[i]);
        v_g += (pos1i[i] > pos2j[i]);
      }
      
      if(v_g > 0){ // blocks necessarily discordant 
        continue;
      }
      else if(v_e == 0) { // blocks necessarily concordant
        for(int i : ids1i){
          C[i] += n2j;
        }
        for(int j : ids2j){
          C[j] += n1i;
        }
      }
      else if(v_e == 1) { // one-dimensional problem
        int k;
        for (int i = 0; i < d; i++) {
          if (pos1i[i] == pos2j[i]) {
            k = ks[i];
            break;
          }
        }
        
        merge(X, C, ids1i, ids2j, k);
      }
      else if(v_e == d) { // still populates the same box.
        // divideAndConquerLower(X, C, thresh, ids1i, ids2j, ks, median);
        divideAndConquerLower(X, C, thresh, ids1i, ids2j, ks);
      }
      else{ // lower-dimensional problem
        std::vector<int> new_ks;
        for (int i = 0; i < d; i++) {
          if (pos1i[i] == pos2j[i]) new_ks.push_back(ks[i]);
        }
        // divideAndConquerLower(X, C, thresh, ids1i, ids2j, new_ks, median);
        divideAndConquerLower(X, C, thresh, ids1i, ids2j, new_ks);
        
      }
    } 
  } 
}


void divideAndConquer(const std::vector<std::vector<int>>& X,
                      std::vector<int>& C,
                      const int thresh,
                      const std::vector<int> ids) {
                      // const bool median) {
  
  int n = ids.size();
  int d = X[0].size();
  
  if(n <= thresh) {
    conquer(X, C, ids);
    return;
  }
  
  // Split Ids into (at most) 2^d sets relative to the midpoint
  
  
  // compute midpoint
  std::vector<double> mid(d);
  // if(median){
  //   mid = midpoint(X, ids);
  // }else{
  int mi;
  int ma;
  for (int k = 0; k < d; k++) {
    mi = X.size();
    ma = 1;
    for (int i : ids) {
      mi = std::min(mi, X[i][k]);
      ma = std::max(ma, X[i][k]);
    }
    mid[k] = (mi + ma)/2.0;
  }
  // }
  
  // create a partition based on position relative to midpoint 
  std::vector<int> pos_vec(d, 0);
  std::map<std::vector<int>, std::vector<int>> pos_map;
  for (int i : ids) {
    for (int j = 0; j < d; j++) {
      pos_vec[j] = (X[i][j] >= mid[j]) ? 1 : 0; // compute pos vector
    }
    pos_map[pos_vec].push_back(i); // add index i to the corresponding pos vector in the map
  }
  
  
  // loop over all pairs of groups
  for (auto it1 = pos_map.begin(); it1 != pos_map.end(); ++it1) {
    // Access the key and value of the first group
    auto pos1 = it1->first;
    auto ids1 = it1->second;
    int n1 = ids1.size();
    
    // divideAndConquer(X, C, thresh, ids1, median);
    divideAndConquer(X, C, thresh, ids1);
    
    for (auto it2 = std::next(it1); it2 != pos_map.end(); ++it2) {
      
      // Access the key and value of the second group
      auto pos2 = it2->first;
      auto ids2 = it2->second;
      int n2 = ids2.size();
      
      // record relation between pos1 and pos2
      int v_l = 0;
      int v_e = 0;
      int v_g = 0;
      for(int i=0; i<pos1.size(); i++){
        v_l += (pos1[i] < pos2[i]);
        v_e += (pos1[i] == pos2[i]);
        v_g += (pos1[i] > pos2[i]);
      }
      
      if((v_l > 0) && (v_g > 0)){ // blocks that are necessarily discordant
        continue;
      }
      else if(v_e == 0) { // blocks necessarily concordant
        for(int i=0; i < n1; i++){
          C[ids1[i]] += n2;
        }
        for(int j=0; j < n2; j++){
          C[ids2[j]] += n1;
        }
      }
      else if(v_e == 1) { // one-dimensional problem
        int k;
        for (int i = 0; i < pos1.size(); i++) {
          if (pos1[i] == pos2[i]) {
            k = i;
            break;
          }
        }
        
        // Note that std::map behaves so that if pos1[k] < pos2[k] here 
        merge(X, C, ids1, ids2, k);
      }
      else{ // k-dimensional problem (k > 1)
        std::vector<int> ks;
        for (int i = 0; i < pos1.size(); i++) {
          if (pos1[i] == pos2[i]) {
            ks.push_back(i);
          }
        }
        
        // Note that std::map behaves so that if pos1[k] <= pos2[k] here 
        // divideAndConquerLower(X, C, thresh, ids1, ids2, ks, median);
        divideAndConquerLower(X, C, thresh, ids1, ids2, ks);
      } 
    } 
  } 
} 


// [[Rcpp::export]]
// IntegerVector dac(IntegerMatrix X, int thresh = 100, bool brute_force = false, bool median = true) {
IntegerVector dac(IntegerMatrix X, int thresh = 100, bool brute_force = false) {
  int n = X.nrow();
  int d = X.ncol();
  
  // initialize C
  std::vector<int> C(n, 0);
  
  // initialize X_vec
  std::vector<std::vector<int>> X_vec(n, std::vector<int>(d));
  
  // fill X_vec
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      X_vec[i][j] = X(i, j);
    }
  }  
  
  std::vector<int> ids(n);
  std::iota(ids.begin(), ids.end(), 0);  // initialize ids to 0, 1, ..., n-1
  
  if(brute_force){
    conquer(X_vec, C, ids);
  }else{
    // divideAndConquer(X_vec, C, thresh, ids, median);
    divideAndConquer(X_vec, C, thresh, ids);
  }
  
  // convert C and return it
  return Rcpp::wrap(C);
}







//* 
//* Functions to compute median as midpoints (not worth it.)
//* 

// std::vector<double> midpoint(const std::vector<std::vector<int>>& X,
//                              const std::vector<int>& ids) {
//   
//   int n = ids.size();
//   int d = X[0].size();  
//   std::vector<int> x(n);
//   std::vector<double> medians(d);
//   
//   for (int k = 0; k < d; k++) {
//     for (int i=0; i<n; i++) {
//       x[i] = X[ids[i]][k];
//     }
//     
//     // Find the median using std::nth_element
//     int midIndex = n / 2;
//     std::nth_element(x.begin(), x.begin() + midIndex, x.end());
//     
//     // Store the median in the medians array
//     medians[k] = x[midIndex];
//   }
//   
//   return(medians);
// }
// 
// 
// std::vector<double> midpoint2(const std::vector<std::vector<int>>& X,
//                               const std::vector<int>& ids1,
//                               const std::vector<int>& ids2,
//                               const std::vector<int>& ks) {
//   
//   int d = ks.size();  
//   int n1 = ids1.size();
//   int n2 = ids2.size();
//   std::vector<int> x1(n1);
//   std::vector<int> x2(n2);
//   std::vector<double> medians(d);
//   
//   for (int k = 0; k < d; k++) {
//     for (int i=0; i<n1; i++) x1[i] = X[ids1[i]][ks[k]];
//     for (int i=0; i<n2; i++) x2[i] = X[ids2[i]][ks[k]];
//     
//     // Find the median using std::nth_element
//     int i1 = n1 / 2;
//     int i2 = n2 / 2;
//     std::nth_element(x1.begin(), x1.begin() + i1, x1.end());
//     std::nth_element(x2.begin(), x2.begin() + i2, x2.end());
//     
//     // Store the median in the medians array
//     medians[k] = (x1[i1] + x2[i2])/2;
//   }
//   
//   return(medians);
// }

