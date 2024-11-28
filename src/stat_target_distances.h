#ifndef STATTARGETDISTANCES_H
#define STATTARGETDISTANCES_H

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

// static target distances
// functions are overloaded with double, float, int, and unsigned short int distances and int and unsigned short int targets
template <typename T, typename D> // T: targets type, D: distances type
void stat_target_distances(const std::vector<D>& vertex_distance, const std::vector<T>& targets, const int starting_index, const int n_targets,
  const int begin_target, const int exclude_index, const bool show_progress, std::vector<D>& distances) {
  // cases
  // 1: pairwise (check starting_indices for the respective targets)
  // 2: no targets and directed (all starts to all starts except self)
  // 3: no targets and not directed (check the starting indices)
  // 4: not pairwise (n_targets in each iteration)
  
  // pairwise
  if(n_targets != -1) {
    #pragma omp simd
    for(int i = 0; i < n_targets; ++i) {
      const int starting_index_i = starting_index + i;
      distances[starting_index_i] = vertex_distance[targets[starting_index_i]];
    }
    if(show_progress) {
      #pragma omp critical(stdrcout)
      Rcpp::Rcout << std::string(n_targets, '=');
    }
  // no targets and directed
  } else if(exclude_index != -1) {
    int n_i = targets.size();
    #pragma omp simd
    for(int i = 0; i < n_i; ++i) {
      if(i != exclude_index) {
        const int j = (i < exclude_index) ? i : (i - 1);
        distances[starting_index + j] = vertex_distance[targets[i]];
      }
    }
    --n_i;
    if(show_progress) {
      #pragma omp critical(stdrcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  // no targets and not directed
  } else if(begin_target != -1) {
    const int n_i = targets.size() - begin_target;
    #pragma omp simd
    for(int i = 0; i < n_i; ++i) {
      distances[starting_index + i] = vertex_distance[targets[begin_target + i]];
    }
    if(show_progress) {
      #pragma omp critical(stdrcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  // not pairwise
  } else {
    const int n_i = targets.size();
    #pragma omp simd
    for(int i = 0; i < n_i; ++i) {
      distances[starting_index + i] = vertex_distance[targets[i]];
    }
    if(show_progress) {
      #pragma omp critical(stdrcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  }
}

#endif
