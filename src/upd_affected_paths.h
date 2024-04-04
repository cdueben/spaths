#ifndef UPDAFFECTEDPATHS_H
#define UPDAFFECTEDPATHS_H

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

std::vector<int> upd_affected_paths_i(const std::unordered_set<int>& upd_rst, const std::vector<std::vector<int> >& paths, const int ncores);
std::vector<unsigned short int> upd_affected_paths_u(const std::unordered_set<int>& upd_rst, const std::vector<std::vector<int> >& paths,
  const int ncores);
std::vector<int> upd_affected_paths_i(const std::unordered_set<unsigned short int>& upd_rst, const std::vector<std::vector<unsigned short int> >& paths,
  const int ncores);
std::vector<unsigned short int> upd_affected_paths_u(const std::unordered_set<unsigned short int>& upd_rst,
  const std::vector<std::vector<unsigned short int> >& paths, const int ncores);

#endif
