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
#include "upd_starts_targets_map.h"
#include "paths_multiple_calls_wweights.h"
#include "repeat_distances.h"
#include "show_progress.h"
#include "upd_paths_wweights.h"

// paths with precomputed weights and grid updating
// functions are overloaded with int and unsigned int cell numbers and double, float, int, and unsigned short int distances
// void upd_paths_wweights
// void upd_paths_wweights
// void upd_paths_wweights
// void upd_paths_wweights
// void upd_paths_wweights
// void upd_paths_wweights
// void upd_paths_wweights
// void upd_paths_wweights

void upd_paths_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<double> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool directed,
  const bool early_stopping, const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress,
  const int bar_limit, std::vector<std::vector<int> >& static_paths, std::vector<std::vector<std::vector<int> > >& upd_paths,
  std::vector<double>& distances) {
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
    show_progress, bar, static_paths, distances);
  
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          // canonical loop instead of range-based for loop to meet OpenMP requirement
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      // To do: test efficiency of deleting the unordered_set
      std::unordered_set<int>().swap(upd_rst[u]);
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<int>().swap(upd_rst[u]);
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

void upd_paths_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<double> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool directed, const bool early_stopping, const int ncores, const bool par_lvl_upd,
  std::vector<std::unordered_set<unsigned short int> >& upd_rst, const bool show_progress, const int bar_limit,
  std::vector<std::vector<unsigned short int> >& static_paths, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths,
  std::vector<double>& distances) {
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<unsigned short int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
    show_progress, bar, static_paths, distances);
  
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
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<unsigned short int>().swap(upd_rst[u]);
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
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<unsigned short int>().swap(upd_rst[u]);
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

void upd_paths_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<float> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool directed,
  const bool early_stopping, const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress,
  const int bar_limit, std::vector<std::vector<int> >& static_paths, std::vector<std::vector<std::vector<int> > >& upd_paths,
  std::vector<float>& distances) {
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
    show_progress, bar, static_paths, distances);
  
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<int>().swap(upd_rst[u]);
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<int>().swap(upd_rst[u]);
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

void upd_paths_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<float> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool directed, const bool early_stopping, const int ncores, const bool par_lvl_upd,
  std::vector<std::unordered_set<unsigned short int> >& upd_rst, const bool show_progress, const int bar_limit,
  std::vector<std::vector<unsigned short int> >& static_paths, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths,
  std::vector<float>& distances) {
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<unsigned short int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
    show_progress, bar, static_paths, distances);
  
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
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<unsigned short int>().swap(upd_rst[u]);
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
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<unsigned short int>().swap(upd_rst[u]);
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

void upd_paths_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<int> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool directed,
  const bool early_stopping, const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress,
  const int bar_limit, std::vector<std::vector<int> >& static_paths, std::vector<std::vector<std::vector<int> > >& upd_paths, std::vector<int>& distances) {
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
    show_progress, bar, static_paths, distances);
  
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<int>().swap(upd_rst[u]);
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<int>().swap(upd_rst[u]);
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

void upd_paths_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<int> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool directed, const bool early_stopping, const int ncores, const bool par_lvl_upd,
  std::vector<std::unordered_set<unsigned short int> >& upd_rst, const bool show_progress, const int bar_limit,
  std::vector<std::vector<unsigned short int> >& static_paths, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths,
  std::vector<int>& distances) {
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<unsigned short int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
    show_progress, bar, static_paths, distances);
  
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
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<unsigned short int>().swap(upd_rst[u]);
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
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<unsigned short int>().swap(upd_rst[u]);
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

void upd_paths_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const std::vector<int>& starts, const std::vector<int>& targets, const std::vector<int>& starting_indices, const bool pairwise,
  const bool directed, const bool early_stopping, const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<int> >& upd_rst,
  const bool show_progress, const int bar_limit, std::vector<std::vector<int> >& static_paths, std::vector<std::vector<std::vector<int> > >& upd_paths,
  std::vector<unsigned short int>& distances) {
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
    show_progress, bar, static_paths, distances);
  
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<int>().swap(upd_rst[u]);
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
        std::unordered_map<int, std::vector<int> > u_starts_targets;
        std::unordered_map<int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<int>().swap(upd_rst[u]);
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

void upd_paths_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const std::vector<int>& starting_indices, const bool pairwise, const bool directed, const bool early_stopping, const int ncores, const bool par_lvl_upd,
  std::vector<std::unordered_set<unsigned short int> >& upd_rst, const bool show_progress, const int bar_limit,
  std::vector<std::vector<unsigned short int> >& static_paths, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths,
  std::vector<unsigned short int>& distances) {
  
  const int n_upd_rst = upd_rst.size();
  const int n_paths = static_paths.size();
  bool bar = show_progress && (n_paths <= bar_limit);
  
  // set used to identify visited targets
  const std::unordered_set<unsigned short int> graph_to_0(graph_to[0].begin(), graph_to[0].end());
  
  // static paths and distances
  paths_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, true, starting_indices, graph_to_0,
    show_progress, bar, static_paths, distances);
  
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
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<unsigned short int>().swap(upd_rst[u]);
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
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_starts_targets;
        std::unordered_map<unsigned short int, std::vector<unsigned short int> > u_affected_paths;
        upd_st_map(affected_paths, starts, targets, pairwise, directed, starting_indices, u_starts_targets, u_affected_paths);
        
        const int ncores_p = (!par_lvl_upd && u_starts_targets.size() == 1) ? ncores : 1;
        const int starting_index = n_paths * (u + 1);
        
        if(ncores_u == 1 || ncores_p != 1) {
          for(const auto & [s, t] : u_starts_targets) {
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, t, early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        } else {
          auto keys_views = std::views::keys(u_starts_targets);
          const std::vector<unsigned short int> keys{keys_views.begin(), keys_views.end()};
          const std::size_t n_keys = keys.size();
          #pragma omp parallel for num_threads(ncores_u) schedule(dynamic)
          for(std::size_t i = 0; i < n_keys; ++i) {
            const unsigned short int s = keys[i];
            paths_multiple_calls_wweights(graph_to, graph_weights, n_cells, s, u_starts_targets[s], early_stopping, ncores_p, graph_to_0, upd_rst[u], bar,
              u_affected_paths[s], upd_paths[u], distances, starting_index);
          }
        }
      }
      std::unordered_set<unsigned short int>().swap(upd_rst[u]);
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
