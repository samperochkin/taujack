#include <Rcpp.h>
using namespace Rcpp;

#include <string.h>
#include <vector>
using namespace std;
// #include <stdio.h>

// This code is loosely based on functions from the pcaPP package.
// See https://github.com/cran/pcaPP/blob/master/src/cov.kendall.cpp

void insertionSort(double* y, int* ids, int* swapCounts, size_t len) {
  size_t maxJ, i;
  
  if(len < 2) {
    return;
  }
  
  maxJ = len - 1;
  for(i = len - 2; i < len; --i) {
    size_t j = i;
    
    double tempVal = y[i];
    int tempId = ids[i];
    for(; j < maxJ && y[j + 1] < tempVal; ++j) {
      y[j] = y[j + 1];
      ids[j] = ids[j + 1];
      swapCounts[ids[j + 1]] += 1;
    }
    
    y[j] = tempVal;
    ids[j] = tempId;
    swapCounts[tempId] += (j - i);
  }
  
}


static void merge(double* y, double* x,
                  int* ids, int* swapCounts,
                  size_t middle, size_t len) {
  
  vector<int> copyIds;
  double* left;
  double* right;
  size_t bufIndex, leftLen, rightLen;
  
  for(size_t k=0; k<len; ++k) copyIds.push_back(ids[k]);
  left = y;
  right = y + middle;
  leftLen = middle;
  rightLen = len - middle;
  bufIndex = 0;
  
  while(leftLen && rightLen) {
    
    if(right[0] < left[0]) {
      x[bufIndex] = right[0];
      ids[bufIndex] = copyIds[len - rightLen];
      swapCounts[copyIds[len - rightLen]] += leftLen;
      
      rightLen--;
      right++;
    } else {
      x[bufIndex] = left[0];
      ids[bufIndex] = copyIds[middle - leftLen];
      swapCounts[copyIds[middle - leftLen]] += len - middle - rightLen;
      
      leftLen--;
      left++;
    }
    bufIndex++;
  }
  
  if(leftLen) {
    for(size_t k=0; k<leftLen; ++k){
      ids[bufIndex + k] = copyIds[middle - leftLen + k];
      swapCounts[copyIds[middle - leftLen + k]] += len - middle;
    }
    memcpy(x + bufIndex, left, leftLen * sizeof(double));
  } else if(rightLen) {
    memcpy(x + bufIndex, right, rightLen * sizeof(double));
  }
  
  return;
}


void mergeSort(double* y, double* x,
               int* ids, int* swapCounts,
               size_t len) {
  
  if(len < 2) {
    return;
  }
  if(len < 10) {
    insertionSort(y, ids, swapCounts, len);
    return;
  }
  
  size_t half;
  half = len / 2;
  mergeSort(y, x, ids, swapCounts, half);
  mergeSort(y + half, x + half, ids + half, swapCounts, len - half);
  merge(y, x, ids, swapCounts, half, len);
  
  memcpy(y, x, len * sizeof(double));
  return;
}



// [[Rcpp::export]]
NumericVector countSwaps(NumericVector y) {
  
  int n = y.size();
  
  // Allocate memory for arrays
  double* y_ptr = new double[n];
  double* x_ptr = new double[n];
  int* ids_ptr = new int[n];
  int* swapCounts_ptr = new int[n];
  
  // Initialize arrays
  for (int i = 0; i < n; i++) {
    y_ptr[i] = y[i];
    x_ptr[i] = 0;
    ids_ptr[i] = i;
    swapCounts_ptr[i] = 0;
  }
  
  // Call mergeSort
  mergeSort(y_ptr, x_ptr, ids_ptr, swapCounts_ptr, n);
  
  // Copy sorted array to output vector
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    out[i] = swapCounts_ptr[i];
  }
  
  // Free memory
  delete[] y_ptr;
  delete[] x_ptr;
  delete[] ids_ptr;
  delete[] swapCounts_ptr;
  
  return out;
}
