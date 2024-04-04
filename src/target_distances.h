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
// the functions are defined in the header file because of their inline nature

inline void target_distances(const std::vector<double>& vertex_distance, const std::vector<int>& affected_paths, const std::vector<int>& targets,
  const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<double>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<float>& vertex_distance, const std::vector<int>& affected_paths, const std::vector<int>& targets,
  const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<float>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<int>& vertex_distance, const std::vector<int>& affected_paths, const std::vector<int>& targets,
  const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<int>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<int>& affected_paths,
  const std::vector<int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<unsigned short int>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<double>& vertex_distance, const std::vector<int>& affected_paths,
  const std::vector<unsigned short int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index,
  const bool show_progress, std::vector<double>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<float>& vertex_distance, const std::vector<int>& affected_paths,
  const std::vector<unsigned short int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index,
  const bool show_progress, std::vector<float>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<int>& vertex_distance, const std::vector<int>& affected_paths,
  const std::vector<unsigned short int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index,
  const bool show_progress, std::vector<int>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<int>& affected_paths,
  const std::vector<unsigned short int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index,
  const bool show_progress, std::vector<unsigned short int>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<double>& vertex_distance, const std::vector<unsigned short int>& affected_paths,
  const std::vector<int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<double>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<float>& vertex_distance, const std::vector<unsigned short int>& affected_paths,
  const std::vector<int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<float>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<int>& vertex_distance, const std::vector<unsigned short int>& affected_paths,
  const std::vector<int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<int>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<unsigned short int>& affected_paths,
  const std::vector<int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<unsigned short int>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<double>& vertex_distance, const std::vector<unsigned short int>& affected_paths,
  const std::vector<unsigned short int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index,
  const bool show_progress, std::vector<double>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<float>& vertex_distance, const std::vector<unsigned short int>& affected_paths,
  const std::vector<unsigned short int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index,
  const bool show_progress, std::vector<float>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<int>& vertex_distance, const std::vector<unsigned short int>& affected_paths,
  const std::vector<unsigned short int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index,
  const bool show_progress, std::vector<int>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

inline void target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<unsigned short int>& affected_paths,
  const std::vector<unsigned short int>& targets, const int starting_index, const int n_targets, const int begin_target, const int exclude_index,
  const bool show_progress, std::vector<unsigned short int>& distances) {
  if(affected_paths.empty()) {
    stat_target_distances(vertex_distance, targets, starting_index, n_targets, begin_target, exclude_index, show_progress, distances);
  } else {
    upd_target_distances(vertex_distance, targets, starting_index, affected_paths, distances);
  }
}

#endif
