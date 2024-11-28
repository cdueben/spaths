#ifndef UPDTARGETPATHS_H
#define UPDTARGETPATHS_H

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
#include <type_traits>

// target paths in updated grids
// functions are overloaded with int and unsigned int predecessor, start, targets, and paths, int and unsigned short int affected_paths, and no, int and
// unsigned short int graph_to_0

// with early stopping and all targets confirmed visited
template <typename P, typename A> // P: predecessor type, A: affected_paths type
void upd_target_paths(const std::vector<P>& predecessor, const P start, const std::vector<P>& targets, const int ncores,
  const std::vector<A>& affected_paths, std::vector<std::vector<P> >& paths) {
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const A starting_index_i = affected_paths[i];
    P pre = targets[i];
    while(pre != start) {
      paths[starting_index_i].push_back(pre);
      pre = predecessor[pre];
    }
    paths[starting_index_i].push_back(start);
  }
}

// without targets confirmed visited
template <typename G, typename A> // P: graph_to_0 type, A: affected_paths type
void upd_target_paths(const std::unordered_set<G>& graph_to_0, const std::vector<G>& predecessor, const G start, const std::vector<G>& targets,
  const int ncores, const std::vector<A>& affected_paths, std::vector<std::vector<G> >& paths) {
  
  constexpr G inf = (std::is_same_v<G, int>) ? -1 : std::numeric_limits<unsigned short int>::max();
  
  const std::size_t n_i = targets.size();
  #pragma omp parallel for num_threads(ncores) if(ncores != 1)
  for(std::size_t i = 0; i < n_i; ++i) {
    const A starting_index_i = affected_paths[i];
    G pre = targets[i];
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

#endif
