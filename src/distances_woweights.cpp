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
#include "distances_one_call_woweights.h"
#include "distances_multiple_calls_woweights.h"
#include "show_progress.h"
#include "distances_woweights.h"

// distances without precomputed weights
// dists_one_call_woweights: this function is only called once and marks vertices as visited by modifying an input vector
// unlike in distances_wweights these functions do not differ in parameter types, which requires varying, not overloaded function names
// void dists_woweights
// void dists_woweights
// void dists_woweights
// void dists_woweights
// void dists_woweights
// void dists_woweights
// void dists_woweights
// void dists_woweights

void dists_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<double>& distances,
  const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    dists_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, bar, distances, radius2);
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void dists_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<float>& distances,
  const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    dists_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, bar, distances, radius2);
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void dists_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<int>& distances,
  const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    dists_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, bar, distances, radius2);
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void dists_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<unsigned short int>& distances,
  const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    dists_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, bar, distances, radius2);
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<double>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    dists_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, bar, distances, radius2);
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<float>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    dists_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, bar, distances, radius2);
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<int>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    dists_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, bar, distances, radius2);
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

void dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<unsigned short int>& distances, const double radius2) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<unsigned short int> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    dists_one_call_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[0], targets, early_stopping, haversine, bar, distances, radius2);
  } else {
    // pairwise
    if(pairwise) {
      const int n_starts_1 = n_starts - 1;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, n_targets);
        }
      }
    // no targets
    } else if(targets.empty()) {
      --n_starts;
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          const int begin_target = i + 1;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], starts, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index, -1, begin_target);
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, true, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts[i], targets, false, haversine, upd_rst, bar, affected_paths,
            distances, radius2, starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}
