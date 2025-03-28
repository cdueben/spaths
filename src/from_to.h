#ifndef FROMTO_H
#define FROMTO_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include "structs.h"

// shift "to" cell numbers (return -1 for NA cells)
inline int shift_to(Rcpp::IntegerVector& cell_numbers, int to, const int n_not_na_cells) {
  int left = 0;
  int right = n_not_na_cells - 1;
  
  while(left <= right) {
    int center = left + (right - left) / 2;
    
    if(cell_numbers[center] == to) {
      return center;
    }
    
    if(cell_numbers[center] < to) {
      left = center + 1;
    } else {
      right = center - 1;
    }
  }
  
  return -1;
}

#endif
