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
#include "coordinates.h"
#include "paths_woweights.h"
#include "upd_affected_paths.h"
#include "upd_starts_targets_map.h"
#include "distances_multiple_calls_woweights.h"
#include "repeat_distances.h"
#include "show_progress.h"
#include "upd_distances_woweights.h"

// distances without precomputed weights and with grid updating
// functions are overloaded with int and unsigned short int graph_to and double, float, int, and unsigned short int distances
// void upd_dists_woweights
// void upd_dists_woweights
// void upd_dists_woweights
// void upd_dists_woweights
// void upd_dists_woweights
// void upd_dists_woweights
// void upd_dists_woweights
// void upd_dists_woweights

void upd_dists_woweights(std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const std::vector<int>& starts, const std::vector<int>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine, const int ncores, const bool par_lvl_upd,
  const std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress, const int bar_limit, std::vector<double>& distances,
  const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords);
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
  
    {
      // static paths and distances
      std::vector<std::vector<int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);

        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
      std::vector<std::vector<int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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

void upd_dists_woweights(std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const std::vector<int>& starts, const std::vector<int>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine, const int ncores, const bool par_lvl_upd,
  const std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress, const int bar_limit, std::vector<float>& distances,
  const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords);
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
  
    {
      // static paths and distances
      std::vector<std::vector<int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_i(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
      
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);

        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
      std::vector<std::vector<int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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

void upd_dists_woweights(std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const std::vector<int>& starts, const std::vector<int>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine, const int ncores, const bool par_lvl_upd,
  const std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress, const int bar_limit, std::vector<int>& distances, const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords);
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
  
    {
      // static paths and distances
      std::vector<std::vector<int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_i(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
      
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);

        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
      std::vector<std::vector<int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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

void upd_dists_woweights(std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const std::vector<int>& starts, const std::vector<int>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine, const int ncores, const bool par_lvl_upd,
  const std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress, const int bar_limit, std::vector<unsigned short int>& distances,
  const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords);
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
  
    {
      // static paths and distances
      std::vector<std::vector<int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_i(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
      
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);

        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
      std::vector<std::vector<int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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

void upd_dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping,
  const bool haversine, const int ncores, const bool par_lvl_upd, const std::vector<std::unordered_set<unsigned short int> >& upd_rst,
  const bool show_progress, const int bar_limit, std::vector<double>& distances, const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords);
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<unsigned short int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
  
    {
      // static paths and distances
      std::vector<std::vector<unsigned short int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_i(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
      
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);

        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
      std::vector<std::vector<unsigned short int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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

void upd_dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping,
  const bool haversine, const int ncores, const bool par_lvl_upd, const std::vector<std::unordered_set<unsigned short int> >& upd_rst,
  const bool show_progress, const int bar_limit, std::vector<float>& distances, const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords);
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<unsigned short int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
  
    {
      // static paths and distances
      std::vector<std::vector<unsigned short int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_i(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
      
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);

        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
      std::vector<std::vector<unsigned short int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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

void upd_dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping,
  const bool haversine, const int ncores, const bool par_lvl_upd, const std::vector<std::unordered_set<unsigned short int> >& upd_rst,
  const bool show_progress, const int bar_limit, std::vector<int>& distances, const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords);
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<unsigned short int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
  
    {
      // static paths and distances
      std::vector<std::vector<unsigned short int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_i(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
      
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);

        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
      std::vector<std::vector<unsigned short int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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

void upd_dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping,
  const bool haversine, const int ncores, const bool par_lvl_upd, const std::vector<std::unordered_set<unsigned short int> >& upd_rst,
  const bool show_progress, const int bar_limit, std::vector<unsigned short int>& distances, const double radius2) {
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords);
  
  const int n_paths = distances.size();
  const int n_upd_rst = upd_rst.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<unsigned short int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  const int ncores_u = (par_lvl_upd) ? 1 : ncores;
  
  if(n_paths * (n_upd_rst + 1) > ((int) USHRT_MAX)) {
    std::vector<std::vector<int> > affected_paths(n_upd_rst);
  
    {
      // static paths and distances
      std::vector<std::vector<unsigned short int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_i(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
      
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);

        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
      std::vector<std::vector<unsigned short int> > paths(n_paths);
      paths_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, true, starting_indices,
        graph_to_0, show_progress, bar, paths, distances, radius2);
      
      // paths affected by upd_rst
      #pragma omp parallel for simd num_threads(ncores) schedule(dynamic) if(ncores != 1)
      for(int u = 0; u < n_upd_rst; ++u) {
        affected_paths[u] = upd_affected_paths_u(upd_rst[u], paths, ncores_u);
      }
    }
    
    if(show_progress) {
      bar = n_upd_rst <= bar_limit;
      upd_show_progress_header(n_upd_rst, false, bar);;
    }
    
    // repeat the static distances
    repeat_distances(distances, n_upd_rst);
    
    // updated distances of affected paths
    #pragma omp parallel for num_threads(ncores) schedule(dynamic) if(par_lvl_upd && ncores != 1)
    for(int u = 0; u < n_upd_rst; ++u) {
      if(!affected_paths[u].empty()) {
        
        // starts, targets, and paths addressed by iterations of the subsequent loop
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths[u], starts, targets, pairwise, false, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, t, early_stopping, haversine, upd_rst[u], false,
              u_affected_paths[s], distances, radius2, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            dists_multiple_calls_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, s, u_starts_targets[s], early_stopping, haversine, upd_rst[u],
              false, u_affected_paths[s], distances, radius2, starting_index);
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
