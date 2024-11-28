#ifndef UPDTARGETDISTANCES_H
#define UPDTARGETDISTANCES_H

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

// target distances in updated grids
// functions are overloaded with double, float, int, and unsigned short int weights, int and unsigned short int targets,
// and int and unsigned short int affected_paths
template <typename D, typename T, typename A> // D: vertex_distance type, T: targets type, A: affected_paths type
void upd_target_distances(const std::vector<D>& vertex_distance, const std::vector<T>& targets, const int starting_index,
  const std::vector<A>& affected_paths, std::vector<D>& distances) {
  const A starting_index_c = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index_c + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

#endif
