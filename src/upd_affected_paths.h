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

// paths affected by upd_rst
// functions are overloaded with int and unsigned short int affected paths and int and unsigned short int upd_rst and paths
// upd_rst and paths contain cell indices (with a possible range from 0 to n_cells - 1)
// affected_paths contains path indices (with a possible range from 0 to n_paths - 1)
template <typename A, typename U> // A: affected_paths type, U: upd_rst type
std::vector<A> upd_affected_paths(const std::unordered_set<U>& upd_rst, const std::vector<std::vector<U> >& paths, const int ncores, const A t) {
  const A n_paths = paths.size();
  std::vector<A> affected_paths;
  
  if(ncores == 1) {
    for(A p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const U & c : paths[p]) {
        if(upd_rst.contains(c)) {
          found = true;
          break;
        }
      }
      if(found) {
        affected_paths.push_back(p);
      }
    }
  } else {
    #pragma omp parallel for num_threads(ncores) schedule(dynamic)
    for(A p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const U & c : paths[p]) {
        if(upd_rst.contains(c)) {
          found = true;
          break;
        }
      }
      if(found) {
        #pragma omp critical(apupdate)
        {
          affected_paths.push_back(p);
        }
      }
    }
  }
  return affected_paths;
}

template <typename U> // U: upd_rst type
std::vector<int> upd_affected_paths_i(const std::unordered_set<U>& upd_rst, const std::vector<std::vector<U> >& paths, const int ncores) {
  constexpr int t {};
  return upd_affected_paths(upd_rst, paths, ncores, t);
}

template <typename U> // U: upd_rst type
std::vector<unsigned short int> upd_affected_paths_u(const std::unordered_set<U>& upd_rst, const std::vector<std::vector<U> >& paths, const int ncores) {
  constexpr unsigned short int t {};
  return upd_affected_paths(upd_rst, paths, ncores, t);
}

#endif
