#ifndef UPDDISTANCESWWEIGHTS_H
#define UPDDISTANCESWWEIGHTS_H

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
#include <unordered_map>
#include <cstddef>
#include <limits>
#include <ranges>
#include <iostream>
#include "paths_wweights.h"
#include "upd_affected_paths.h"
#include "repeat_distances.h"
#include "upd_starts_targets_map.h"
#include "distances_multiple_calls_wweights.h"
#include "show_progress.h"

// distances with grid updating
// functions are overloaded with int and unsigned short int graph_to adjacency lists and double, float, int, and unsigned short int distances
template <typename G, typename D> // G: graph_to type, D: distances type
void upd_dists_wweights(const std::vector<std::vector<G> >& graph_to, std::vector<std::vector<D> >& graph_weights, const std::size_t n_cells,
  const std::vector<G>& starts, const std::vector<G>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool directed,
  const bool early_stopping, const int ncores, const bool par_lvl_upd, const std::vector<std::unordered_set<G> >& upd_rst, const bool show_progress,
  const int bar_limit, std::vector<D>& distances) {
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<G> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
    
    {
      // static paths and distances
      std::vector<std::vector<G> > paths(n_paths);
      paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
        show_progress, bar, paths, distances);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_i(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<G, std::vector<G> > u_starts_targets;
        std::unordered_map<G, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, upd_rst[u], false, u_affected_paths[s], distances,
              starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<G> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const G s = keys[i];
            dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, upd_rst[u], false, u_affected_paths[s],
              distances, starting_index);
          }
        }
      }
      if(bar) {
        #pragma omp critical(urcout)
        Rcpp::Rcout << '=';
      }
    }
  } else {
    std::vector<std::vector<unsigned short int> > affected_paths(n_upd_rst);
    
    {
      // static paths and distances
      std::vector<std::vector<G> > paths(n_paths);
      paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
        show_progress, bar, paths, distances);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<G, std::vector<G> > u_starts_targets;
        std::unordered_map<G, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, upd_rst[u], false, u_affected_paths[s], distances,
              starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<G> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const G s = keys[i];
            dists_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, upd_rst[u], false, u_affected_paths[s],
              distances, starting_index);
          }
        }
      }
      if(bar) {
        #pragma omp critical(urcout)
        Rcpp::Rcout << '=';
      }
    }
  }
  if(bar) {
    Rcpp::Rcout << '|' << std::endl;
  }
}

#endif
