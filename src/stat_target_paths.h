#ifndef STATTARGETPATHS_H
#define STATTARGETPATHS_H

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
#include <unordered_set>
#include <cstddef>

void stat_target_paths(const std::vector<int>& predecessor, const int start, const std::vector<int>& targets, const int ncores, const int starting_index,
  const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<std::vector<int> >& paths);
void stat_target_paths(const std::vector<unsigned short int>& predecessor, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const int ncores, const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<std::vector<unsigned short int> >& paths);
void stat_target_paths(const std::unordered_set<int>& graph_to_0, const std::vector<int>& predecessor, const int start, const std::vector<int>& targets,
  const int ncores, const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<std::vector<int> >& paths);
void stat_target_paths(const std::unordered_set<unsigned short int>& graph_to_0, const std::vector<unsigned short int>& predecessor,
  const unsigned short int start, const std::vector<unsigned short int>& targets, const int ncores, const int starting_index, const int n_targets,
  const int begin_target, const int exclude_index, const bool show_progress, std::vector<std::vector<unsigned short int> >& paths);

#endif
