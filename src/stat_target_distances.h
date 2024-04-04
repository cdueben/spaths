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

void stat_target_distances(const std::vector<double>& vertex_distance, const std::vector<int>& targets, const int starting_index, const int n_targets,
  const int begin_target, const int exclude_index, const bool show_progress, std::vector<double>& distances);
void stat_target_distances(const std::vector<float>& vertex_distance, const std::vector<int>& targets, const int starting_index, const int n_targets,
  const int begin_target, const int exclude_index, const bool show_progress, std::vector<float>& distances);
void stat_target_distances(const std::vector<int>& vertex_distance, const std::vector<int>& targets, const int starting_index, const int n_targets,
  const int begin_target, const int exclude_index, const bool show_progress, std::vector<int>& distances);
void stat_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<unsigned short int>& distances);
void stat_target_distances(const std::vector<double>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<double>& distances);
void stat_target_distances(const std::vector<float>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<float>& distances);
void stat_target_distances(const std::vector<int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<int>& distances);
void stat_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<unsigned short int>& distances);

#endif
