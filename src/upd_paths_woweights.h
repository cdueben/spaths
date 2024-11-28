#ifndef UPDPATHSWOWEIGHTS_H
#define UPDPATHSWOWEIGHTS_H

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
#include "paths_woweights.h"
#include "upd_affected_paths.h"
#include "upd_starts_targets_map.h"
#include "paths_multiple_calls_woweights.h"
#include "repeat_distances.h"
#include "show_progress.h"
#include "coordinates.h"

// paths without precomputed weights and grid updating
// functions are overloaded with int and unsigned int cell numbers and double, float, int, and unsigned short int distances
template <typename G, typename D> // G: graph_to type, D: distances type
void upd_paths_woweights(std::vector<std::vector<G> >& graph_to, Rcpp::List& coords, const std::vector<G>& starts, const std::vector<G>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine, const int ncores, const bool par_lvl_upd,
  std::vector<std::unordered_set<G> >& upd_rst, const bool show_progress, const int bar_limit, std::vector<std::vector<G> >& static_paths,
  std::vector<std::vector<std::vector<G> > >& upd_paths, std::vector<D>& distances, const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords, false);
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<G> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
    graph_to_0, show_progress, bar, static_paths, distances, radius2);
  
  if(show_progress) {
    bar = n_upd_rst <= bar_limit;
    upd_show_progress_header(n_upd_rst, true, bar);
  }
  
  // repeat the static distances
  repeat_distances(distances, n_upd_rst);
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  const int ap_n_upd_rst = (distances.empty()) ? 0 : n_upd_rst;
  
  if(n_paths * (ap_n_upd_rst + 1) > ((int) USHRT_MAX)) {
    // updated paths of affected connections
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      std::vector<int> affected_paths = upd_affected_paths_i(upd_rst[u], static_paths, ncores_u);
      if(!affected_paths.empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<G, std::vector<G> > u_starts_targets;
        std::unordered_map<G, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, ncores_p, graph_to_0,
              upd_rst[u], false, u_affected_paths[s], upd_paths[u], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<G> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const G s = keys[i];
            paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, ncores_p,
              graph_to_0, upd_rst[u], false, u_affected_paths[s], upd_paths[u], distances, radius2, starting_index);
          }
        }
      }
      std::unordered_set<G>().swap(upd_rst[u]);
      if(bar) {
        #pragma omp critical(urcout)
        Rcpp::Rcout << '=';
      }
    }
  } else {
    // updated paths of affected connections
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      std::vector<unsigned short int> affected_paths = upd_affected_paths_u(upd_rst[u], static_paths, ncores_u);
      if(!affected_paths.empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<G, std::vector<G> > u_starts_targets;
        std::unordered_map<G, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, ncores_p, graph_to_0,
              upd_rst[u], false, u_affected_paths[s], upd_paths[u], distances, radius2, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<G> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const G s = keys[i];
            paths_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, ncores_p,
              graph_to_0, upd_rst[u], false, u_affected_paths[s], upd_paths[u], distances, radius2, starting_index);
          }
        }
      }
      std::unordered_set<G>().swap(upd_rst[u]);
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
