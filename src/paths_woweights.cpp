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
#include "paths_one_call_woweights.h"
#include "paths_multiple_calls_woweights.h"
#include "show_progress.h"
#include "paths_woweights.h"

// paths without precomputed weights
// functions are overloaded with int and unsigned short int adjacency lists and double, float, int, and unsigned short int distances
// void paths_woweights
// void paths_woweights
// void paths_woweights
// void paths_woweights
// void paths_woweights
// void paths_woweights
// void paths_woweights
// void paths_woweights

void paths_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices, const std::unordered_set<int>& graph_to_0,
  const bool show_progress, const bool bar, std::vector<std::vector<int> >& paths, std::vector<double>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    if(upd_rst_defined) {
      paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0,
        upd_rst, bar, affected_paths, paths, distances, radius2, 0);
    } else {
      paths_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0, bar,
        paths, distances, radius2);
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
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index, n_targets);
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        const int begin_target = i + 1;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, early_stopping, haversine, 1, graph_to_0, upd_rst,
          bar, affected_paths, paths, distances, radius2, starting_index, -1, begin_target);
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices, const std::unordered_set<int>& graph_to_0,
  const bool show_progress, const bool bar, std::vector<std::vector<int> >& paths, std::vector<float>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    if(upd_rst_defined) {
      paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0,
        upd_rst, bar, affected_paths, paths, distances, radius2, 0);
    } else {
      paths_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0, bar,
        paths, distances, radius2);
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
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index, n_targets);
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        const int begin_target = i + 1;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, early_stopping, haversine, 1, graph_to_0, upd_rst,
          bar, affected_paths, paths, distances, radius2, starting_index, -1, begin_target);
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices, const std::unordered_set<int>& graph_to_0,
  const bool show_progress, const bool bar, std::vector<std::vector<int> >& paths, std::vector<int>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    if(upd_rst_defined) {
      paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0,
        upd_rst, bar, affected_paths, paths, distances, radius2, 0);
    } else {
      paths_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0, bar,
        paths, distances, radius2);
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
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index, n_targets);
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        const int begin_target = i + 1;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, early_stopping, haversine, 1, graph_to_0, upd_rst,
          bar, affected_paths, paths, distances, radius2, starting_index, -1, begin_target);
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices, const std::unordered_set<int>& graph_to_0,
  const bool show_progress, const bool bar, std::vector<std::vector<int> >& paths, std::vector<unsigned short int>& distances,
  const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    if(upd_rst_defined) {
      paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0,
        upd_rst, bar, affected_paths, paths, distances, radius2, 0);
    } else {
      paths_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0, bar,
        paths, distances, radius2);
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
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index, n_targets);
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        const int begin_target = i + 1;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, early_stopping, haversine, 1, graph_to_0, upd_rst,
          bar, affected_paths, paths, distances, radius2, starting_index, -1, begin_target);
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const bool upd_rst_defined,
  const std::vector<int>& starting_indices, const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, const bool bar,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<double>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    if(upd_rst_defined) {
      paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0,
        upd_rst, bar, affected_paths, paths, distances, radius2, 0);
    } else {
      paths_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0, bar,
        paths, distances, radius2);
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
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index, n_targets);
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        const int begin_target = i + 1;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, early_stopping, haversine, 1, graph_to_0, upd_rst,
          bar, affected_paths, paths, distances, radius2, starting_index, -1, begin_target);
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const bool upd_rst_defined,
  const std::vector<int>& starting_indices, const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, const bool bar,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<float>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    if(upd_rst_defined) {
      paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0,
        upd_rst, bar, affected_paths, paths, distances, radius2, 0);
    } else {
      paths_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0, bar,
        paths, distances, radius2);
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
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index, n_targets);
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        const int begin_target = i + 1;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, early_stopping, haversine, 1, graph_to_0, upd_rst,
          bar, affected_paths, paths, distances, radius2, starting_index, -1, begin_target);
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const bool upd_rst_defined,
  const std::vector<int>& starting_indices, const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, const bool bar,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<int>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    if(upd_rst_defined) {
      paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0,
        upd_rst, bar, affected_paths, paths, distances, radius2, 0);
    } else {
      paths_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0, bar,
        paths, distances, radius2);
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
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index, n_targets);
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        const int begin_target = i + 1;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, early_stopping, haversine, 1, graph_to_0, upd_rst,
          bar, affected_paths, paths, distances, radius2, starting_index, -1, begin_target);
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void paths_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const bool upd_rst_defined,
  const std::vector<int>& starting_indices, const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, const bool bar,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<unsigned short int>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(paths.size(), upd_rst_defined, true, bar);
  }
  
  if(n_starts == 1) {
    if(upd_rst_defined) {
      paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0,
        upd_rst, bar, affected_paths, paths, distances, radius2, 0);
    } else {
      paths_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, ncores, graph_to_0, bar,
        paths, distances, radius2);
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
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index, n_targets);
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = starting_indices[i];
        const int begin_target = i + 1;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, early_stopping, haversine, 1, graph_to_0, upd_rst,
          bar, affected_paths, paths, distances, radius2, starting_index, -1, begin_target);
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int i = 0; i < n_starts; ++i) {
        const int starting_index = i * n_targets;
        paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, early_stopping, haversine, 1, graph_to_0,
          upd_rst, bar, affected_paths, paths, distances, radius2, starting_index);
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}
