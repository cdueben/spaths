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
#include <iostream>
#include "paths_one_call_wweights.h"
#include "paths_multiple_calls_wweights.h"
#include "show_progress.h"
#include "paths_wweights.h"

// static paths with precomputed weights
// functions are overloaded with int and unsigned short int cell numbers and double, float, int, and unsigned short int distances
// void paths_wweights
// void paths_wweights
// void paths_wweights
// void paths_wweights
// void paths_wweights
// void paths_wweights
// void paths_wweights
// void paths_wweights

void paths_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<double> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const std::unordered_set<int>& graph_to_0, const bool show_progress, const bool bar,
  std::vector<std::vector<int> >& paths, std::vector<double>& distances) {
  // cases
  // 1: pairwise (check starting_indices for the respective targets)
  // 2: no targets and directed (all starts to all starts except self)
  // 3: no targets and not directed (check the starting indices)
  // 4: not pairwise (n_targets in each iteration)
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, upd_rst, bar, affected_paths,
        paths, distances, 0);
    } else {
      paths_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, bar, paths, distances);
    }
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
        n_targets -= starting_index;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index, n_targets);
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_starts_1;
          const int exclude_index = i;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, -1, exclude_index);
        }
      // no targets and not directed
      } else {
        --n_starts;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<float> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const std::unordered_set<int>& graph_to_0, const bool show_progress, const bool bar,
  std::vector<std::vector<int> >& paths, std::vector<float>& distances) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, upd_rst, bar, affected_paths,
        paths, distances, 0);
    } else {
      paths_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, bar, paths, distances);
    }
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
        n_targets -= starting_index;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index, n_targets);
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_starts_1;
          const int exclude_index = i;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, -1, exclude_index);
        }
      // no targets and not directed
      } else {
        --n_starts;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<int> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const std::unordered_set<int>& graph_to_0, const bool show_progress, const bool bar,
  std::vector<std::vector<int> >& paths, std::vector<int>& distances) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, upd_rst, bar, affected_paths,
        paths, distances, 0);
    } else {
      paths_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, bar, paths, distances);
    }
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
        n_targets -= starting_index;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index, n_targets);
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_starts_1;
          const int exclude_index = i;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, -1, exclude_index);
        }
      // no targets and not directed
      } else {
        --n_starts;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const std::unordered_set<int>& graph_to_0, const bool show_progress, const bool bar,
  std::vector<std::vector<int> >& paths, std::vector<unsigned short int>& distances) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, upd_rst, bar, affected_paths,
        paths, distances, 0);
    } else {
      paths_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, bar, paths, distances);
    }
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
        n_targets -= starting_index;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index, n_targets);
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_starts_1;
          const int exclude_index = i;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, -1, exclude_index);
        }
      // no targets and not directed
      } else {
        --n_starts;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<double> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const int ncores, const bool directed, const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices,
  const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, const bool bar, std::vector<std::vector<unsigned short int> >& paths,
  std::vector<double>& distances) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, upd_rst, bar, affected_paths,
        paths, distances, 0);
    } else {
      paths_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, bar, paths, distances);
    }
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
        n_targets -= starting_index;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index, n_targets);
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_starts_1;
          const int exclude_index = i;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, -1, exclude_index);
        }
      // no targets and not directed
      } else {
        --n_starts;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<float> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const int ncores, const bool directed, const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices,
  const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, const bool bar, std::vector<std::vector<unsigned short int> >& paths,
  std::vector<float>& distances) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, upd_rst, bar, affected_paths,
        paths, distances, 0);
    } else {
      paths_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, bar, paths, distances);
    }
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
        n_targets -= starting_index;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index, n_targets);
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_starts_1;
          const int exclude_index = i;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, -1, exclude_index);
        }
      // no targets and not directed
      } else {
        --n_starts;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<int> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const int ncores, const bool directed, const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices,
  const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, const bool bar, std::vector<std::vector<unsigned short int> >& paths,
  std::vector<int>& distances) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, upd_rst, bar, affected_paths,
        paths, distances, 0);
    } else {
      paths_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, bar, paths, distances);
    }
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
        n_targets -= starting_index;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index, n_targets);
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_starts_1;
          const int exclude_index = i;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, -1, exclude_index);
        }
      // no targets and not directed
      } else {
        --n_starts;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const int ncores, const bool directed, const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices,
  const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, const bool bar, std::vector<std::vector<unsigned short int> >& paths,
  std::vector<unsigned short int>& distances) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, upd_rst, bar, affected_paths,
        paths, distances, 0);
    } else {
      paths_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, ncores, graph_to_0, bar, paths, distances);
    }
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
        n_targets -= starting_index;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index, n_targets);
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_starts_1;
          const int exclude_index = i;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, -1, exclude_index);
        }
      // no targets and not directed
      } else {
        --n_starts;
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
            paths, distances, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, early_stopping, 1, graph_to_0, upd_rst, bar, affected_paths,
          paths, distances, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}
