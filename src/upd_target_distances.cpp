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
#include "upd_target_distances.h"

// target distances in updated grids
// functions are overloaded with double, float, int, and unsigned short in weights, int and unsigned short int targets,
// and int and unsigned short int affected_paths
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances
// void upd_target_distances

void upd_target_distances(const std::vector<double>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<double>& distances) {
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<float>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<float>& distances) {
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<int>& distances) {
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<unsigned short int>& distances) {
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<double>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<double>& distances) {
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<float>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<float>& distances) {
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<int>& distances) {
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<unsigned short int>& distances) {
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<double>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<double>& distances) {
  const unsigned short int u_starting_index = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[u_starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<float>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<float>& distances) {
  const unsigned short int u_starting_index = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[u_starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<int>& distances) {
  const unsigned short int u_starting_index = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[u_starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<unsigned short int>& distances) {
  const unsigned short int u_starting_index = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[u_starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<double>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<double>& distances) {
  const unsigned short int u_starting_index = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[u_starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<float>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<float>& distances) {
  const unsigned short int u_starting_index = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[u_starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<int>& distances) {
  const unsigned short int u_starting_index = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[u_starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}

void upd_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<unsigned short int>& distances) {
  const unsigned short int u_starting_index = starting_index;
  
  const std::size_t n_i = targets.size();
  #pragma omp simd
  for(std::size_t i = 0; i < n_i; ++i) {
    distances[u_starting_index + affected_paths[i]] = vertex_distance[targets[i]];
  }
}
