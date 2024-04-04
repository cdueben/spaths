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
#include <limits>
#include <unordered_set>
#include <cstddef>
#include "upd_target_paths.h"

// target paths in updated grids
// functions are overloaded with int and unsigned int predecessor, start, targets, and paths, int and unsigned short int affected_paths, and no, int and
// unsigned short int graph_to_0

// with early stopping and all targets confirmed visited
void upd_target_paths(const std::vector<int>& predecessor, const int start, const std::vector<int>& targets, const int ncores,
  const std::vector<int>& affected_paths, std::vector<std::vector<int> >& paths) {
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const int starting_index_i = affected_paths[i];
    int pre = targets[i];
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}

void upd_target_paths(const std::vector<unsigned short int>& predecessor, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const int ncores, const std::vector<int>& affected_paths, std::vector<std::vector<unsigned short int> >& paths) {
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const int starting_index_i = affected_paths[i];
    unsigned short int pre = targets[i];
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}

void upd_target_paths(const std::vector<int>& predecessor, const int start, const std::vector<int>& targets, const int ncores,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths) {
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const unsigned short int starting_index_i = affected_paths[i];
    int pre = targets[i];
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}

void upd_target_paths(const std::vector<unsigned short int>& predecessor, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const int ncores, const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<unsigned short int> >& paths) {
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const unsigned short int starting_index_i = affected_paths[i];
    unsigned short int pre = targets[i];
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}

// without targets confirmed visited
void upd_target_paths(const std::unordered_set<int>& graph_to_0, const std::vector<int>& predecessor, const int start, const std::vector<int>& targets,
  const int ncores, const std::vector<int>& affected_paths, std::vector<std::vector<int> >& paths) {
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const int starting_index_i = affected_paths[i];
    int pre = targets[i];
    if(predecessor[pre] == 0 && !graph_to_0.contains(pre)) {
      paths[starting_index_i].push_back(-1);
      continue;
    }
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}

void upd_target_paths(const std::unordered_set<unsigned short int>& graph_to_0, const std::vector<unsigned short int>& predecessor,
  const unsigned short int start, const std::vector<unsigned short int>& targets, const int ncores, const std::vector<int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths) {
  
  const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const int starting_index_i = affected_paths[i];
    unsigned short int pre = targets[i];
    if(predecessor[pre] == 0 && !graph_to_0.contains(pre)) {
      paths[starting_index_i].push_back(inf);
      continue;
    }
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}

void upd_target_paths(const std::unordered_set<int>& graph_to_0, const std::vector<int>& predecessor, const int start, const std::vector<int>& targets,
  const int ncores, const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths) {
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const unsigned short int starting_index_i = affected_paths[i];
    int pre = targets[i];
    if(predecessor[pre] == 0 && !graph_to_0.contains(pre)) {
      paths[starting_index_i].push_back(-1);
      continue;
    }
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}

void upd_target_paths(const std::unordered_set<unsigned short int>& graph_to_0, const std::vector<unsigned short int>& predecessor,
  const unsigned short int start, const std::vector<unsigned short int>& targets, const int ncores, const std::vector<unsigned short int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths) {
  
  const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const unsigned short int starting_index_i = affected_paths[i];
    unsigned short int pre = targets[i];
    if(predecessor[pre] == 0 && !graph_to_0.contains(pre)) {
      paths[starting_index_i].push_back(inf);
      continue;
    }
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}
