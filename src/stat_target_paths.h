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

// static paths between start and targets
// with early stopping and all targets confirmed visited
template <typename P> // P: predecessors type
void stat_target_paths(const std::vector<P>& predecessor, const P start, const std::vector<P>& targets, const int ncores, const int starting_index,
  const int n_targets, const int begin_target, const int exclude_index, const bool show_progress, std::vector<std::vector<P> >& paths) {
  
  // pairwise
  if(n_targets != -1) {
    #pragma omp parallel for num_threads(ncores) if(ncores != 1)
    for(int i = 0; i < n_targets; ++i) {
      const int starting_index_i = starting_index + i;
      P pre = targets[starting_index_i];
      while(pre != start) {
        paths[starting_index_i].push_back(pre);
        pre = predecessor[pre];
      }
      paths[starting_index_i].push_back(start);
    }
    if(show_progress) {
      #pragma omp critical(stprcout)
      Rcpp::Rcout << std::string(n_targets, '=');
    }
  // no targets and directed
  } else if(exclude_index != -1) {
    int n_i = targets.size();
    #pragma omp parallel for num_threads(ncores) if(ncores != 1)
    for(int i = 0; i < n_i; ++i) {
      if(i != exclude_index) {
        int starting_index_i = (i < exclude_index) ? i : (i - 1);
        starting_index_i = starting_index + starting_index_i;
        P pre = targets[i];
        while(pre != start) {
          paths[starting_index_i].push_back(pre);
          pre = predecessor[pre];
        }
        paths[starting_index_i].push_back(start);
      }
    }
    --n_i;
    if(show_progress) {
      #pragma omp critical(stprcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  // no targets and not directed
  } else if(begin_target != -1) {
    const int n_i = targets.size() - begin_target;
    #pragma omp parallel for num_threads(ncores) if(ncores != 1)
    for(int i = 0; i < n_i; ++i) {
      const int starting_index_i = starting_index + i;
      P pre = targets[begin_target + i];
      while(pre != start) {
        paths[starting_index_i].push_back(pre);
        pre = predecessor[pre];
      }
      paths[starting_index_i].push_back(start);
    }
    if(show_progress) {
      #pragma omp critical(stprcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  // not pairwise
  } else {
    const int n_i = targets.size();
    #pragma omp parallel for num_threads(ncores) if(ncores != 1)
    for(int i = 0; i < n_i; ++i) {
      const int starting_index_i = starting_index + i;
      P pre = targets[i];
      while(pre != start) {
        paths[starting_index_i].push_back(pre);
        pre = predecessor[pre];
      }
      paths[starting_index_i].push_back(start);
    }
    if(show_progress) {
      #pragma omp critical(stprcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  }
}

// without targets confirmed visited
template <typename G> // G: graph_to type
void stat_target_paths(const std::unordered_set<G>& graph_to_0, const std::vector<G>& predecessor, const G start, const std::vector<G>& targets,
  const int ncores, const int starting_index, const int n_targets, const int begin_target, const int exclude_index, const bool show_progress,
  std::vector<std::vector<G> >& paths) {
  
  // pairwise
  if(n_targets != -1) {
    #pragma omp parallel for num_threads(ncores) if(ncores != 1)
    for(int i = 0; i < n_targets; ++i) {
      const int starting_index_i = starting_index + i;
      G pre = targets[starting_index_i];
      // skip unvisited target
      if(predecessor[pre] == 0 && !graph_to_0.contains(pre)) {
        continue;
      }
      while(pre != start) {
        paths[starting_index_i].push_back(pre);
        pre = predecessor[pre];
      }
      paths[starting_index_i].push_back(start);
    }
    if(show_progress) {
      #pragma omp critical(stprcout)
      Rcpp::Rcout << std::string(n_targets, '=');
    }
  // no targets and directed
  } else if(exclude_index != -1) {
    int n_i = targets.size();
    #pragma omp parallel for num_threads(ncores) if(ncores != 1)
    for(int i = 0; i < n_i; ++i) {
      if(i != exclude_index) {
        G pre = targets[i];
        if(predecessor[pre] == 0 && !graph_to_0.contains(pre)) {
          continue;
        }
        int starting_index_i = (i < exclude_index) ? i : (i - 1);
        starting_index_i = starting_index + starting_index_i;
        while(pre != start) {
          paths[starting_index_i].push_back(pre);
          pre = predecessor[pre];
        }
        paths[starting_index_i].push_back(start);
      }
    }
    --n_i;
    if(show_progress) {
      #pragma omp critical(stprcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  // no targets and not directed
  } else if(begin_target != -1) {
    const int n_i = targets.size() - begin_target;
    #pragma omp parallel for num_threads(ncores) if(ncores != 1)
    for(int i = 0; i < n_i; ++i) {
      G pre = targets[begin_target + i];
      if(predecessor[pre] == 0 && !graph_to_0.contains(pre)) {
        continue;
      }
      const int starting_index_i = starting_index + i;
      while(pre != start) {
        paths[starting_index_i].push_back(pre);
        pre = predecessor[pre];
      }
      paths[starting_index_i].push_back(start);
    }
    if(show_progress) {
      #pragma omp critical(stprcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  // not pairwise
  } else {
    const int n_i = targets.size();
    #pragma omp parallel for num_threads(ncores) if(ncores != 1)
    for(int i = 0; i < n_i; ++i) {
      G pre = targets[i];
      if(predecessor[pre] == 0 && !graph_to_0.contains(pre)) {
        continue;
      }
      const int starting_index_i = starting_index + i;
      while(pre != start) {
        paths[starting_index_i].push_back(pre);
        pre = predecessor[pre];
      }
      paths[starting_index_i].push_back(start);
    }
    if(show_progress) {
      #pragma omp critical(stprcout)
      Rcpp::Rcout << std::string(n_i, '=');
    }
  }
}

#endif
