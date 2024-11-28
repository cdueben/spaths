#ifndef TARGETSSET_H
#define TARGETSSET_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_set>

// targets sets checked in the early_stopping case

template <typename T> // T: targets type
inline std::unordered_set<T> create_targets_set(const std::vector<T>& targets, const bool affected_paths_empty, const int starting_index,
  const int n_targets, const int begin_target, const int exclude_index) {
  // static paths
  if(affected_paths_empty) {
    // pairwise
    if(n_targets != -1) {
      return std::unordered_set<T>(targets.begin() + starting_index, targets.begin() + starting_index + n_targets);
    // no targets and directed
    } else if(exclude_index != -1) {
      std::unordered_set<T> targets_set;
      const int targets_size_1 = targets.size() - 1;
      targets_set.reserve(targets_size_1);
      if(exclude_index != 0) {
        targets_set.insert(targets.begin(), targets.begin() + (exclude_index - 1));
      }
      if(exclude_index != targets_size_1) {
        targets_set.insert(targets.begin() + (exclude_index + 1), targets.end());
      }
      return targets_set;
    // no targets and not directed
    } else if(begin_target != -1) {
      return std::unordered_set<T>(targets.begin() + begin_target, targets.end());
    // not pairwise
    } else {
      return std::unordered_set<T>(targets.begin(), targets.end());
    }
  // updated paths
  } else {
    return std::unordered_set<T>(targets.begin(), targets.end());
  }
}

#endif
