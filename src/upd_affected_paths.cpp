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
#include "upd_affected_paths.h"

// paths affected by upd_rst
// functions are overloaded with int and unsigned short int affected paths and int and unsigned short int upd_rst and paths and 
// upd_rst and paths contain cell indices (with a possible range from 0 to n_cells - 1)
// affected_paths contains path indices (with a possible range from 0 to n_paths - 1)
// std::vector<int> upd_affected_paths_i
// std::vector<unsigned short int> upd_affected_paths_u
// std::vector<int> upd_affected_paths_i
// std::vector<unsigned short int> upd_affected_paths_u

std::vector<int> upd_affected_paths_i(const std::unordered_set<int>& upd_rst, const std::vector<std::vector<int> >& paths, const int ncores) {
  const int n_paths = paths.size();
  std::vector<int> affected_paths;
  
  if(ncores == 1) {
    for(int p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const int & c : paths[p]) {
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
    for(int p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const int & c : paths[p]) {
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

std::vector<unsigned short int> upd_affected_paths_u(const std::unordered_set<int>& upd_rst, const std::vector<std::vector<int> >& paths,
  const int ncores) {
  const unsigned short int n_paths = paths.size();
  std::vector<unsigned short int> affected_paths;
  
  if(ncores == 1) {
    for(unsigned short int p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const int & c : paths[p]) {
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
    for(unsigned short int p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const int & c : paths[p]) {
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

std::vector<int> upd_affected_paths_i(const std::unordered_set<unsigned short int>& upd_rst, const std::vector<std::vector<unsigned short int> >& paths,
  const int ncores) {
  const int n_paths = paths.size();
  std::vector<int> affected_paths;
  
  if(ncores == 1) {
    for(int p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const unsigned short int & c : paths[p]) {
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
    for(int p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const unsigned short int & c : paths[p]) {
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

std::vector<unsigned short int> upd_affected_paths_u(const std::unordered_set<unsigned short int>& upd_rst,
  const std::vector<std::vector<unsigned short int> >& paths, const int ncores) {
  const unsigned short int n_paths = paths.size();
  std::vector<unsigned short int> affected_paths;
  
  if(ncores == 1) {
    for(unsigned short int p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const unsigned short int & c : paths[p]) {
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
    for(unsigned short int p = 0; p < n_paths; ++p) {
      bool found {false};
      for(const unsigned short int & c : paths[p]) {
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
