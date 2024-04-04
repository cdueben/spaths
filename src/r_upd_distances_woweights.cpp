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
#include "starts_targets.h"
#include "convert_upd_rst.h"
#include "graph_to.h"
#include "graph_weights.h"
#include "upd_distances_wweights.h"
#include "upd_distances_woweights.h"

// distances without precomputed weights and with grid updating
// distances that are exported to R
// graphs are not constructed in parallel as the from and to R vectors are not necessarily thread-safe on all machines
// Rcpp::NumericVector r_upd_dists_woweights_d
// Rcpp::IntegerVector r_upd_dists_woweights_i

// [[Rcpp::export]]
Rcpp::NumericVector r_upd_dists_woweights_d(Rcpp::List& from_to, Rcpp::List& starts_targets, Rcpp::List& coords, const std::size_t n_cells,
  Rcpp::List& upd_rst_r, const bool haversine, const bool queen, const int ncores, const bool par_lvl_upd, const bool pairwise, const bool pre,
  const bool early_stopping, const bool int_path, const bool double_weights, const bool show_progress, const int bar_limit, const double radius2) {
  
  const int n_paths = compute_n_paths(starts_targets, false, pairwise);
  
  // with this structure, upd_rst_c, graph_to, and graph_weights go out of scope before the R object is assembled (alternative to .swap() etc.)
  if(double_weights) {
    std::vector<double> distances (n_paths);
    if(int_path) {
      const std::vector<int> starts = get_starts_i(starts_targets);
      const std::vector<int> targets = get_targets_i(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
      const std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
      // with precomputed weights
      if(pre) {
        std::vector<std::vector<double> > graph_weights = graph_weights_d(graph_to, coords, haversine, queen, ncores, true, radius2);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, false, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      // without precomputed weights
      } else {
        upd_dists_woweights(graph_to, coords, starts, targets, starting_indices, pairwise, early_stopping, haversine, ncores, par_lvl_upd, upd_rst_c,
          show_progress, bar_limit, distances, radius2);
      }
    } else {
      const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
      const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
      const std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
      if(pre) {
        std::vector<std::vector<double> > graph_weights = graph_weights_d(graph_to, coords, haversine, queen, ncores, true, radius2);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, false, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      } else {
        upd_dists_woweights(graph_to, coords, starts, targets, starting_indices, pairwise, early_stopping, haversine, ncores, par_lvl_upd, upd_rst_c,
          show_progress, bar_limit, distances, radius2);
      }
    }
    return Rcpp::wrap(distances);
  } else {
    std::vector<float> distances (n_paths);
    if(int_path) {
      const std::vector<int> starts = get_starts_i(starts_targets);
      const std::vector<int> targets = get_targets_i(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
      const std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
      // with precomputed weights
      if(pre) {
        std::vector<std::vector<float> > graph_weights = graph_weights_f(graph_to, coords, haversine, queen, ncores, true, radius2);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, false, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      // without precomputed weights
      } else {
        upd_dists_woweights(graph_to, coords, starts, targets, starting_indices, pairwise, early_stopping, haversine, ncores, par_lvl_upd, upd_rst_c,
          show_progress, bar_limit, distances, radius2);
      }
    } else {
      const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
      const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
      const std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
      if(pre) {
        std::vector<std::vector<float> > graph_weights = graph_weights_f(graph_to, coords, haversine, queen, ncores, true, radius2);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, false, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      } else {
        upd_dists_woweights(graph_to, coords, starts, targets, starting_indices, pairwise, early_stopping, haversine, ncores, par_lvl_upd, upd_rst_c,
          show_progress, bar_limit, distances, radius2);
      }
    }
    return Rcpp::wrap(distances);
  }
}

// [[Rcpp::export]]
Rcpp::IntegerVector r_upd_dists_woweights_i(Rcpp::List& from_to, Rcpp::List& starts_targets, Rcpp::List& coords, const std::size_t n_cells,
  Rcpp::List& upd_rst_r, const bool haversine, const bool queen, const int ncores, const bool par_lvl_upd, const bool pairwise, const bool pre,
  const bool early_stopping, const bool int_path, const bool signed_weights, const bool show_progress, const int bar_limit, const double radius2) {
  
  const int n_paths = compute_n_paths(starts_targets, false, pairwise);
  
  // with this structure, upd_rst_c, graph_to, and graph_weights go out of scope before the R object is assembled (alternative to .swap() etc.)
  if(signed_weights) {
    std::vector<int> distances (n_paths);
    if(int_path) {
      const std::vector<int> starts = get_starts_i(starts_targets);
      const std::vector<int> targets = get_targets_i(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
      const std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
      // with precomputed weights
      if(pre) {
        std::vector<std::vector<int> > graph_weights = graph_weights_i(graph_to, coords, haversine, queen, ncores, true, radius2);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, false, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      // without precomputed weights
      } else {
        upd_dists_woweights(graph_to, coords, starts, targets, starting_indices, pairwise, early_stopping, haversine, ncores, par_lvl_upd, upd_rst_c,
          show_progress, bar_limit, distances, radius2);
      }
    } else {
      const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
      const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
      const std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
      if(pre) {
        std::vector<std::vector<int> > graph_weights = graph_weights_i(graph_to, coords, haversine, queen, ncores, true, radius2);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, false, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      } else {
        upd_dists_woweights(graph_to, coords, starts, targets, starting_indices, pairwise, early_stopping, haversine, ncores, par_lvl_upd, upd_rst_c,
          show_progress, bar_limit, distances, radius2);
      }
    }
    return Rcpp::wrap(distances);
  } else {
    std::vector<unsigned short int> distances (n_paths);
    if(int_path) {
      const std::vector<int> starts = get_starts_i(starts_targets);
      const std::vector<int> targets = get_targets_i(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
      const std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
      // with precomputed weights
      if(pre) {
        std::vector<std::vector<unsigned short int> > graph_weights = graph_weights_u(graph_to, coords, haversine, queen, ncores, true, radius2);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, false, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      // without precomputed weights
      } else {
        upd_dists_woweights(graph_to, coords, starts, targets, starting_indices, pairwise, early_stopping, haversine, ncores, par_lvl_upd, upd_rst_c,
          show_progress, bar_limit, distances, radius2);
      }
    } else {
      const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
      const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
      const std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
      if(pre) {
        std::vector<std::vector<unsigned short int> > graph_weights = graph_weights_u(graph_to, coords, haversine, queen, ncores, true, radius2);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, false, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      } else {
        upd_dists_woweights(graph_to, coords, starts, targets, starting_indices, pairwise, early_stopping, haversine, ncores, par_lvl_upd, upd_rst_c,
          show_progress, bar_limit, distances, radius2);
      }
    }
    return Rcpp::wrap(distances);
  }
}
