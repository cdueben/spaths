// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include "starts_targets.h"

// extract R vectors and convert them to C++ vectors
// inline std::vector<int> get_starts_i
// inline std::vector<unsigned short int> get_starts_u
// inline std::vector<int> get_targets_i
// inline std::vector<unsigned short int> get_targets_u
// inline std::vector<int> get_starting_indices_i
// inline std::size_t compute_n_paths
// inline std::size_t compute_n_paths

// format of R inputs
// pairwise: unique starts, targets sorted by start, n_paths_per_start denote how many targets each start links to
// no targets, directed: starts, no targets, no n_paths_per_start
// no targets, not directed: starts, no targets, no n_paths_per_start
// not pairwise: starts, targets, no n_paths_per_start

std::vector<int> get_starting_indices_i(Rcpp::List& starts_targets, int starts_size, const bool no_targets_not_directed, const bool pairwise) {
  if(no_targets_not_directed) {
    --starts_size;
    std::vector<int> starting_indices (starts_size);
    int n_targets = starts_size;
    int starting_index {0};
    for(int i = 1; i < starts_size; ++i) {
      starting_index += n_targets;
      --n_targets;
      starting_indices[i] = starting_index;
    }
    return starting_indices;
  } else if(pairwise) {
    Rcpp::IntegerVector n_paths_per_start = starts_targets["n_paths_per_start"];
    std::vector<int> starting_indices (starts_size);
    --starts_size;
    int starting_index {0};
    for(int i = 0; i < starts_size; ++i) {
      starting_index += n_paths_per_start[i];
      starting_indices[i + 1] = starting_index;
    }
    return starting_indices;
  } else {
    return std::vector<int>();
  }
}
