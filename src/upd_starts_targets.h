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

void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<int>& affected_paths, const std::vector<int>& starts,
  const std::vector<int>& targets, std::vector<int>& upd_starts, std::vector<int>& upd_targets);
void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<unsigned short int>& affected_paths,
  const std::vector<int>& starts, const std::vector<int>& targets, std::vector<int>& upd_starts, std::vector<int>& upd_targets);
void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<int>& affected_paths,
  const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, std::vector<unsigned short int>& upd_starts,
  std::vector<unsigned short int>& upd_targets);
void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<unsigned short int>& affected_paths,
  const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, std::vector<unsigned short int>& upd_starts,
  std::vector<unsigned short int>& upd_targets);
void upd_starts_targets_no_targets_directed(const std::vector<int>& affected_paths, const std::vector<int>& starts, std::vector<int>& upd_starts,
  std::vector<int>& upd_targets);
void upd_starts_targets_no_targets_directed(const std::vector<unsigned short int>& affected_paths, const std::vector<int>& starts,
  std::vector<int>& upd_starts, std::vector<int>& upd_targets);
void upd_starts_targets_no_targets_directed(const std::vector<int>& affected_paths, const std::vector<unsigned short int>& starts,
  std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets);
void upd_starts_targets_no_targets_directed(const std::vector<unsigned short int>& affected_paths, const std::vector<unsigned short int>& starts,
  std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets);
void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<int>& affected_paths,
  const std::vector<int>& starts, std::vector<int>& upd_starts, std::vector<int>& upd_targets);
void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<unsigned short int>& affected_paths,
  const std::vector<int>& starts, std::vector<int>& upd_starts, std::vector<int>& upd_targets);
void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<int>& affected_paths,
  const std::vector<unsigned short int>& starts, std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets);
void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<unsigned short int>& affected_paths,
  const std::vector<unsigned short int>& starts, std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets);
void upd_starts_targets_not_pairwise(const std::vector<int>& affected_paths, const std::vector<int>& starts, const std::vector<int>& targets,
  std::vector<int>& upd_starts, std::vector<int>& upd_targets);
void upd_starts_targets_not_pairwise(const std::vector<unsigned short int>& affected_paths, const std::vector<int>& starts, const std::vector<int>& targets,
  std::vector<int>& upd_starts, std::vector<int>& upd_targets);
void upd_starts_targets_not_pairwise(const std::vector<int>& affected_paths, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets);
void upd_starts_targets_not_pairwise(const std::vector<unsigned short int>& affected_paths, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets);

#endif
