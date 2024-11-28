#ifndef UPDSTARTSTARGETSMAP_H
#define UPDSTARTSTARGETSMAP_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstddef>
#include "upd_starts_targets.h"

// starts linked to the respective targets and paths
// functions are overloaded with int and unsigned short int starts and targets
// functions are overloaded with int and unsigned short int path numbers
template <typename A, typename S> // A: affected_paths, S: starts type
void upd_st_map(const std::vector<A>& affected_paths, const std::vector<S>& starts, const std::vector<S>& targets, const bool pairwise, const bool directed,
  const std::vector<int>& starting_indices, std::unordered_map<S, std::vector<S> >& u_starts_targets,
  std::unordered_map<S, std::vector<A> >& u_affected_paths) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  
  // starts and targets of affected paths
  std::vector<S> upd_starts (n_affected_paths);
  std::vector<S> upd_targets (n_affected_paths);
  if(pairwise) {
    upd_starts_targets_pairwise(starting_indices, affected_paths, starts, targets, upd_starts, upd_targets);
  } else if(targets.empty()) {
    if(directed) {
      upd_starts_targets_no_targets_directed(affected_paths, starts, upd_starts, upd_targets);
    } else {
      upd_starts_targets_no_targets_not_directed(starting_indices, affected_paths, starts, upd_starts, upd_targets);
    }
  } else {
    upd_starts_targets_not_pairwise(affected_paths, starts, targets, upd_starts, upd_targets);
  }
  
  // map with the starts as keys and the respective targets as values
  if(directed || std::unordered_set<S>(upd_starts.begin(), upd_starts.end()).size() < std::unordered_set<S>(upd_targets.begin(),
    upd_targets.end()).size()) {
    for(std::size_t p = 0; p < n_affected_paths; ++p) {
      const S start = upd_starts[p];
      u_starts_targets[start].push_back(upd_targets[p]);
      u_affected_paths[start].push_back(affected_paths[p]);
    }
  // map with the targets as keys and the respective starts as values
  } else {
    for(std::size_t p = 0; p < n_affected_paths; ++p) {
      const S target = upd_targets[p];
      u_starts_targets[target].push_back(upd_starts[p]);
      u_affected_paths[target].push_back(affected_paths[p]);
    }
  }
}

#endif
