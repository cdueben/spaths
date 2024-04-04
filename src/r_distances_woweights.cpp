// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include "graph_to.h"
#include "starts_targets.h"
#include "distances_wweights.h"
#include "distances_woweights.h"
#include "graph_weights.h"

// distances without precomputed weights and update_rst
// Rcpp::NumericVector r_dists_woweights_d
// Rcpp::IntegerVector r_dists_woweights_i

// [[Rcpp::export]]
Rcpp::NumericVector r_dists_woweights_d(Rcpp::List& from_to, Rcpp::List& starts_targets, Rcpp::List& coords, const std::size_t n_cells,
  const bool haversine, const bool queen, const int ncores, const bool pairwise, const bool pre, const bool early_stopping, const bool int_path,
  const bool double_weights, const bool show_progress, const int bar_limit, const double radius2) {
  
  const int n_paths = compute_n_paths(starts_targets, false, pairwise);
  const bool bar = show_progress && (n_paths <= bar_limit);
  
  // with this structure, upd_rst_c, graph_to, and graph_weights go out of scope before the R object is assembled (alternative to .swap() etc.)
  if(double_weights) {
    std::vector<double> distances (n_paths);
    if(int_path) {
      const std::vector<int> starts = get_starts_i(starts_targets);
      const std::vector<int> targets = get_targets_i(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
      // with precomputed weights
      if(pre) {
        std::vector<std::vector<double> > graph_weights = graph_weights_d(graph_to, coords, haversine, queen, ncores, true, radius2);
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, false, pairwise, false, starting_indices, show_progress,
          bar, distances);
      // without precomputed weights
      } else {
        const double xres = coords["xres"];
        const double yres = coords["yres"];
        const int ncol = coords["ncol"];
        const double ymax = coords["ymax"];
        const std::vector<int> cell_numbers = get_cell_numbers(coords);
        dists_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, starting_indices,
          show_progress, bar, distances, radius2);
      }
    } else {
      const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
      const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
      if(pre) {
        std::vector<std::vector<double> > graph_weights = graph_weights_d(graph_to, coords, haversine, queen, ncores, true, radius2);
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, false, pairwise, false, starting_indices, show_progress,
          bar, distances);
      } else {
        const double xres = coords["xres"];
        const double yres = coords["yres"];
        const int ncol = coords["ncol"];
        const double ymax = coords["ymax"];
        const std::vector<int> cell_numbers = get_cell_numbers(coords);
        dists_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, starting_indices,
          show_progress, bar, distances, radius2);
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
      // with precomputed weights
      if(pre) {
        std::vector<std::vector<float> > graph_weights = graph_weights_f(graph_to, coords, haversine, queen, ncores, true, radius2);
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, false, pairwise, false, starting_indices, show_progress,
          bar, distances);
      // without precomputed weights
      } else {
        const double xres = coords["xres"];
        const double yres = coords["yres"];
        const int ncol = coords["ncol"];
        const double ymax = coords["ymax"];
        const std::vector<int> cell_numbers = get_cell_numbers(coords);
        dists_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, starting_indices,
          show_progress, bar, distances, radius2);
      }
    } else {
      const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
      const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
      if(pre) {
        std::vector<std::vector<float> > graph_weights = graph_weights_f(graph_to, coords, haversine, queen, ncores, true, radius2);
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, false, pairwise, false, starting_indices, show_progress,
          bar, distances);
      } else {
        const double xres = coords["xres"];
        const double yres = coords["yres"];
        const int ncol = coords["ncol"];
        const double ymax = coords["ymax"];
        const std::vector<int> cell_numbers = get_cell_numbers(coords);
        dists_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, starting_indices,
          show_progress, bar, distances, radius2);
      }
    }
    return Rcpp::wrap(distances);
  }
}

// [[Rcpp::export]]
Rcpp::IntegerVector r_dists_woweights_i(Rcpp::List& from_to, Rcpp::List& starts_targets, Rcpp::List& coords, const std::size_t n_cells,
  const bool haversine, const bool queen, const int ncores, const bool pairwise, const bool pre, const bool early_stopping, const bool int_path,
  const bool signed_weights, const bool show_progress, const int bar_limit, const double radius2) {
  
  const int n_paths = compute_n_paths(starts_targets, false, pairwise);
  const bool bar = show_progress && (n_paths <= bar_limit);
  
  // with this structure, upd_rst_c, graph_to, and graph_weights go out of scope before the R object is assembled (alternative to .swap() etc.)
  if(signed_weights) {
    std::vector<int> distances (n_paths);
    if(int_path) {
      const std::vector<int> starts = get_starts_i(starts_targets);
      const std::vector<int> targets = get_targets_i(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
      // with precomputed weights
      if(pre) {
        std::vector<std::vector<int> > graph_weights = graph_weights_i(graph_to, coords, haversine, queen, ncores, true, radius2);
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, false, pairwise, false, starting_indices, show_progress,
          bar, distances);
      // without precomputed weights
      } else {
        const double xres = coords["xres"];
        const double yres = coords["yres"];
        const int ncol = coords["ncol"];
        const double ymax = coords["ymax"];
        const std::vector<int> cell_numbers = get_cell_numbers(coords);
        dists_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, starting_indices,
          show_progress, bar, distances, radius2);
      }
    } else {
      const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
      const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
      if(pre) {
        std::vector<std::vector<int> > graph_weights = graph_weights_i(graph_to, coords, haversine, queen, ncores, true, radius2);
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, false, pairwise, false, starting_indices, show_progress,
          bar, distances);
      } else {
        const double xres = coords["xres"];
        const double yres = coords["yres"];
        const int ncol = coords["ncol"];
        const double ymax = coords["ymax"];
        const std::vector<int> cell_numbers = get_cell_numbers(coords);
        dists_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, starting_indices,
          show_progress, bar, distances, radius2);
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
      // with precomputed weights
      if(pre) {
        std::vector<std::vector<unsigned short int> > graph_weights = graph_weights_u(graph_to, coords, haversine, queen, ncores, true, radius2);
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, false, pairwise, false, starting_indices, show_progress,
          bar, distances);
      // without precomputed weights
      } else {
        const double xres = coords["xres"];
        const double yres = coords["yres"];
        const int ncol = coords["ncol"];
        const double ymax = coords["ymax"];
        const std::vector<int> cell_numbers = get_cell_numbers(coords);
        dists_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, starting_indices,
          show_progress, bar, distances, radius2);
      }
    } else {
      const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
      const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
      const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), targets.empty(), pairwise);
      std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
      if(pre) {
        std::vector<std::vector<unsigned short int> > graph_weights = graph_weights_u(graph_to, coords, haversine, queen, ncores, true, radius2);
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, false, pairwise, false, starting_indices, show_progress,
          bar, distances);
      } else {
        const double xres = coords["xres"];
        const double yres = coords["yres"];
        const int ncol = coords["ncol"];
        const double ymax = coords["ymax"];
        const std::vector<int> cell_numbers = get_cell_numbers(coords);
        dists_woweights(graph_to, cell_numbers, ncol, xres, yres, ymax, starts, targets, early_stopping, haversine, ncores, pairwise, starting_indices,
          show_progress, bar, distances, radius2);
      }
    }
    return Rcpp::wrap(distances);
  }
}
