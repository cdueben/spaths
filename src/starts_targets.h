#ifndef STARTSTARGETS_H
#define STARTSTARGETS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>

inline std::vector<int> get_starts_i(Rcpp::List& starts_targets) {
  Rcpp::IntegerVector starts = starts_targets["starts"];
  starts_targets["starts"] = R_NilValue;
  return Rcpp::as<std::vector<int> >(starts);
}

inline std::vector<unsigned short int> get_starts_u(Rcpp::List& starts_targets) {
  Rcpp::IntegerVector starts = starts_targets["starts"];
  starts_targets["starts"] = R_NilValue;
  return Rcpp::as<std::vector<unsigned short int> >(starts);
}

inline std::vector<int> get_targets_i(Rcpp::List& starts_targets) {
  if(starts_targets.containsElementNamed("targets")) {
    Rcpp::IntegerVector targets_r = starts_targets["targets"];
    const std::vector<int> targets_c = Rcpp::as<std::vector<int> >(targets_r);
    starts_targets["targets"] = R_NilValue;
    return targets_c;
  } else {
    return std::vector<int>();
  }
}

inline std::vector<unsigned short int> get_targets_u(Rcpp::List& starts_targets) {
  if(starts_targets.containsElementNamed("targets")) {
    Rcpp::IntegerVector targets_r = starts_targets["targets"];
    const std::vector<unsigned short int> targets_c = Rcpp::as<std::vector<unsigned short int> >(targets_r);
    starts_targets["targets"] = R_NilValue;
    return targets_c;
  } else {
    return std::vector<unsigned short int>();
  }
}

std::vector<int> get_starting_indices_i(Rcpp::List& starts_targets, int starts_size, const bool no_targets_not_directed, const bool pairwise);

inline int compute_n_paths(Rcpp::List& starts_targets, const bool directed, const bool pairwise) {
  const int n_starts = starts_targets["n_starts"];
  const int n_targets = starts_targets["n_targets"];
  int n_paths;
  
  // from starts to specific targets
  if(pairwise) {
    n_paths = n_targets;
  // from all starts to all other starts
  } else if(n_targets == 0) {
    n_paths = n_starts * (n_starts - 1);
    if(!directed) {
      n_paths /= 2;
    }
  // from all starts to all targets
  } else {
    n_paths = n_starts * n_targets;
  }
  
  return n_paths;
}

#endif
