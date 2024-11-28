#ifndef DISTANCESWWEIGHTS_H
#define DISTANCESWWEIGHTS_H

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
#include "distances_one_call_wweights.h"
#include "distances_multiple_calls_wweights.h"
#include "show_progress.h"

// distances with precomputed weights
// functions are overloaded with int and unsigned short int graph_to adjacency lists and double, float, int, and unsigned short int distances
template <typename G, typename D> // G: graph_to type, D: distances type
void dists_wweights(std::vector<std::vector<G> >& graph_to, std::vector<std::vector<D> >& graph_weights, const std::size_t n_cells,
  const std::vector<G>& starts, const std::vector<G>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<D>& distances) {
  
  int n_starts = starts.size();
  const std::vector<int> affected_paths;
  const std::unordered_set<G> upd_rst;
  
  if(show_progress) {
    stat_show_progress_header(distances.size(), false, false, bar);
  }
  
  if(n_starts == 1) {
    // directed graphs rely on the visited vector to mark visited cells
    if(upd_rst_defined || directed) {
      dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, upd_rst, bar, affected_paths, distances, 0);
    } else {
      dists_one_call_wweights(graph_to, graph_weights, n_cells, starts[0], targets, early_stopping, bar, distances);
    }
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
          dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, true, upd_rst, bar, affected_paths, distances, starting_index,
            n_targets);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = starting_indices[i];
          int n_targets = (i == n_starts_1) ? targets.size() : starting_indices[i + 1];
          n_targets -= starting_index;
          dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, false, upd_rst, bar, affected_paths, distances,
            starting_index, n_targets);
        }
      }
    } else if(targets.empty()) {
      // no targets and directed
      if(directed) {
        const int n_starts_1 = n_starts - 1;
        if(early_stopping) {
          #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
          for(int i = 0; i < n_starts; ++i) {
            const int starting_index = i * n_starts_1;
            const int exclude_index = i;
            dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, true, upd_rst, bar, affected_paths, distances,
              starting_index, -1, -1, exclude_index);
          }
        } else {
          #pragma omp parallel for num_threads(ncores) if(ncores != 1)
          for(int i = 0; i < n_starts; ++i) {
            const int starting_index = i * n_starts_1;
            const int exclude_index = i;
            dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, false, upd_rst, bar, affected_paths, distances,
              starting_index, -1, -1, exclude_index);
          }
        }
      // no targets and not directed
      } else {
        --n_starts;
        if(early_stopping) {
          #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
          for(int i = 0; i < n_starts; ++i) {
            const int starting_index = starting_indices[i];
            const int begin_target = i + 1;
            dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, true, upd_rst, bar, affected_paths, distances,
              starting_index, -1, begin_target);
          }
        } else {
          #pragma omp parallel for num_threads(ncores) if(ncores != 1)
          for(int i = 0; i < n_starts; ++i) {
            const int starting_index = starting_indices[i];
            const int begin_target = i + 1;
            dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], starts, false, upd_rst, bar, affected_paths, distances,
              starting_index, -1, begin_target);
          }
        }
      }
    // not pairwise
    } else {
      const int n_targets = targets.size();
      if(early_stopping) {
        #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, true, upd_rst, bar, affected_paths, distances,
            starting_index);
        }
      } else {
        #pragma omp parallel for num_threads(ncores) if(ncores != 1)
        for(int i = 0; i < n_starts; ++i) {
          const int starting_index = i * n_targets;
          dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, starts[i], targets, false, upd_rst, bar, affected_paths, distances,
            starting_index);
        }
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

#endif
