#ifndef TARGETDISTANCES_H
#define TARGETDISTANCES_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include "stat_target_distances.h"
#include "upd_target_distances.h"

// distances of the target cells
// functions are overloaded with double, float, int, and unsigned short int distances, int and unsigned short targets,
// and int and unsigned short int affected_paths

template <typename D, typename A, typename T> // D: vertex_distance type, A: affected_paths, T: targets
inline void target_distances(const std::vector<D>& vertex_distance, const std::vector<A>& affected_paths, const std::vector<T>& targets,
  const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<D>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

#endif
