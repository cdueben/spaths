#ifndef TARGETPATHS_H
#define TARGETPATHS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include "stat_target_paths.h"
#include "upd_target_paths.h"

// header-only functions, because they are inline functions

inline void target_paths(const std::vector<int>& predecessor, const int start, const std::vector<int>& targets, const std::unordered_set<int>& graph_to_0,
  const std::vector<int>& affected_paths, const bool all_visited, const int ncores, const int starting_index, const int n_targets, const int begin_target,
  const int exclude_index, const bool show_progress, std::vector<std::vector<int> >& paths) {
  
  if(all_visited) {
    if(affected_paths.empty()) {
      stat_target_paths(predecessor, start, targets, ncores, starting_index, n_targets, begin_target, exclude_index, show_progress, paths);
    } else {
      upd_target_paths(predecessor, start, targets, ncores, affected_paths, paths);
    }
  } else {
    if(affected_paths.empty()) {
      stat_target_paths(graph_to_0, predecessor, start, targets, ncores, starting_index, n_targets, begin_target, exclude_index, show_progress, paths);
    } else {
      upd_target_paths(graph_to_0, predecessor, start, targets, ncores, affected_paths, paths);
    }
  }
}

inline void target_paths(const std::vector<unsigned short int>& predecessor, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const std::unordered_set<unsigned short int>& graph_to_0, const std::vector<int>& affected_paths, const bool all_visited, const int ncores,
  const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<std::vector<unsigned short int> >& paths) {
  
  if(all_visited) {
    if(affected_paths.empty()) {
      stat_target_paths(predecessor, start, targets, ncores, starting_index, n_targets, begin_target, exclude_index, show_progress, paths);
    } else {
      upd_target_paths(predecessor, start, targets, ncores, affected_paths, paths);
    }
  } else {
    if(affected_paths.empty()) {
      stat_target_paths(graph_to_0, predecessor, start, targets, ncores, starting_index, n_targets, begin_target, exclude_index, show_progress, paths);
    } else {
      upd_target_paths(graph_to_0, predecessor, start, targets, ncores, affected_paths, paths);
    }
  }
}

inline void target_paths(const std::vector<int>& predecessor, const int start, const std::vector<int>& targets, const std::unordered_set<int>& graph_to_0,
  const std::vector<unsigned short int>& affected_paths, const bool all_visited, const int ncores, const int starting_index, const int n_targets,
  const int begin_target, const int exclude_index, const bool show_progress, std::vector<std::vector<int> >& paths) {
  
  if(all_visited) {
    if(affected_paths.empty()) {
      stat_target_paths(predecessor, start, targets, ncores, starting_index, n_targets, begin_target, exclude_index, show_progress, paths);
    } else {
      upd_target_paths(predecessor, start, targets, ncores, affected_paths, paths);
    }
  } else {
    if(affected_paths.empty()) {
      stat_target_paths(graph_to_0, predecessor, start, targets, ncores, starting_index, n_targets, begin_target, exclude_index, show_progress, paths);
    } else {
      upd_target_paths(graph_to_0, predecessor, start, targets, ncores, affected_paths, paths);
    }
  }
}

inline void target_paths(const std::vector<unsigned short int>& predecessor, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const std::unordered_set<unsigned short int>& graph_to_0, const std::vector<unsigned short int>& affected_paths, const bool all_visited, const int ncores,
  const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<std::vector<unsigned short int> >& paths) {
  
  if(all_visited) {
    if(affected_paths.empty()) {
      stat_target_paths(predecessor, start, targets, ncores, starting_index, n_targets, begin_target, exclude_index, show_progress, paths);
    } else {
      upd_target_paths(predecessor, start, targets, ncores, affected_paths, paths);
    }
  } else {
    if(affected_paths.empty()) {
      stat_target_paths(graph_to_0, predecessor, start, targets, ncores, starting_index, n_targets, begin_target, exclude_index, show_progress, paths);
    } else {
      upd_target_paths(graph_to_0, predecessor, start, targets, ncores, affected_paths, paths);
    }
  }
}

#endif
