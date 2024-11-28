#ifndef UPDSTARTSTARGETS_H
#define UPDSTARTSTARGETS_H

// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads()  1
#define omp_get_thread_num()   0
#define omp_get_max_threads()  1
#define omp_get_thread_limit() 1
#define omp_get_num_procs()    1
#endif
#include <vector>
#include <cstddef>

// pairwise and no targets and not directed cases
// starts indices search uses binary search
template <typename A, typename S> // A: affected_paths type, S: starts type
void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<A>& affected_paths, const std::vector<S>& starts,
  const std::vector<S>& targets, std::vector<S>& upd_starts, std::vector<S>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    bool found {false};
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        upd_starts[p] = starts[i];
        found = true;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    if(!found) {
      upd_starts[p] = starts[right];
    }
    upd_targets[p] = targets[p];
  }
}

// no targets and directed
template <typename A, typename S> // A: affected_paths, S: starts type
void upd_starts_targets_no_targets_directed(const std::vector<A>& affected_paths, const std::vector<S>& starts, std::vector<S>& upd_starts,
  std::vector<S>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const A starts_size_1 = starts.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    const A start = affected_paths[p] / starts_size_1;
    A target = affected_paths[p] % starts_size_1;
    if(start <= target) {
      ++target;
    }
    upd_starts[p] = starts[start];
    upd_targets[p] = starts[target];
  }
}

// no targets and not directed
template <typename A, typename S> // A: affected_paths, S: starts type
void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<A>& affected_paths,
  const std::vector<S>& starts, std::vector<S>& upd_starts, std::vector<S>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        right = i;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    upd_starts[p] = starts[right];
    upd_targets[p] = starts[affected_path - starting_indices[right] + right + 1];
  }
}

// not pairwise
template <typename A, typename S> // A: affected_paths, S: starts type
void upd_starts_targets_not_pairwise(const std::vector<A>& affected_paths, const std::vector<S>& starts, const std::vector<S>& targets,
  std::vector<S>& upd_starts, std::vector<S>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const A targets_size = targets.size();
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    upd_starts[p] = starts[affected_paths[p] / targets_size];
    upd_targets[p] = targets[affected_paths[p] % targets_size];
  }
}

#endif
