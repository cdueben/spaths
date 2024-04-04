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

void upd_target_distances(const std::vector<double>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<double>& distances);
void upd_target_distances(const std::vector<float>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<float>& distances);
void upd_target_distances(const std::vector<int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<int>& distances);
void upd_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<unsigned short int>& distances);
void upd_target_distances(const std::vector<double>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<double>& distances);
void upd_target_distances(const std::vector<float>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<float>& distances);
void upd_target_distances(const std::vector<int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<int>& distances);
void upd_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<int>& affected_paths, std::vector<unsigned short int>& distances);
void upd_target_distances(const std::vector<double>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<double>& distances);
void upd_target_distances(const std::vector<float>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<float>& distances);
void upd_target_distances(const std::vector<int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<int>& distances);
void upd_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<unsigned short int>& distances);
void upd_target_distances(const std::vector<double>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<double>& distances);
void upd_target_distances(const std::vector<float>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<float>& distances);
void upd_target_distances(const std::vector<int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<int>& distances);
void upd_target_distances(const std::vector<unsigned short int>& vertex_distance, const std::vector<unsigned short int>& targets, const int starting_index,
  const std::vector<unsigned short int>& affected_paths, std::vector<unsigned short int>& distances);

#endif
