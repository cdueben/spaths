// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <cstddef>
#include <vector>
#include "starts_targets.h"
#include "graph_to.h"
#include "graph_weights.h"
#include "distances_wweights.h"

// distances with precomputed weights and without grid updating
// distances that are exported to R
// Rcpp::NumericVector r_dists_wweights_d
// Rcpp::IntegerVector r_dists_wweights_i

// [[Rcpp::export]]
Rcpp::NumericVector r_dists_wweights_d(Rcpp::List& from_to, Rcpp::List& starts_targets, const std::size_t n_cells, const bool early_stopping,
  const int ncores, const bool directed, const bool pairwise, const bool int_path, const bool double_weights, const bool show_progress,
  const int bar_limit) {
  
  const int n_paths = compute_n_paths(starts_targets, directed, pairwise);
  const bool bar = show_progress && (n_paths <= bar_limit);
  
  if(double_weights) {
    std::vector<double> distances (n_paths);
    {
      std::vector<std::vector<double> > graph_weights = graph_weights_d(from_to, n_cells);
      if(int_path) {
        const std::vector<int> starts = get_starts_i(starts_targets);
        const std::vector<int> targets = get_targets_i(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
        from_to["from"] = R_NilValue;
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, false, starting_indices,
          show_progress, bar, distances);
      } else {
        const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
        const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
        from_to["from"] = R_NilValue;
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, false, starting_indices,
          show_progress, bar, distances);
      }
    }
    return Rcpp::wrap(distances);
  } else {
    std::vector<float> distances (n_paths);
    {
      std::vector<std::vector<float> > graph_weights = graph_weights_f(from_to, n_cells);
      if(int_path) {
        const std::vector<int> starts = get_starts_i(starts_targets);
        const std::vector<int> targets = get_targets_i(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
        from_to["from"] = R_NilValue;
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, false, starting_indices,
          show_progress, bar, distances);
      } else {
        const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
        const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
        from_to["from"] = R_NilValue;
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, false, starting_indices,
          show_progress, bar, distances);
      }
    }
    return Rcpp::wrap(distances);
  }
}

// [[Rcpp::export]]
Rcpp::IntegerVector r_dists_wweights_i(Rcpp::List& from_to, Rcpp::List& starts_targets, const std::size_t n_cells, const bool early_stopping,
  const int ncores, const bool directed, const bool pairwise, const bool int_path, const bool signed_weights, const bool show_progress,
  const int bar_limit) {
  
  const int n_paths = compute_n_paths(starts_targets, directed, pairwise);
  const bool bar = show_progress && (n_paths <= bar_limit);
  
  if(signed_weights) {
    std::vector<int> distances (n_paths);
    {
      std::vector<std::vector<int> > graph_weights = graph_weights_i(from_to, n_cells);
      if(int_path) {
        const std::vector<int> starts = get_starts_i(starts_targets);
        const std::vector<int> targets = get_targets_i(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
        from_to["from"] = R_NilValue;
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, false, starting_indices,
          show_progress, bar, distances);
      } else {
        const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
        const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
        from_to["from"] = R_NilValue;
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, false, starting_indices,
          show_progress, bar, distances);
      }
    }
    return Rcpp::wrap(distances);
  } else {
    std::vector<unsigned short int> distances (n_paths);
    {
      std::vector<std::vector<unsigned short int> > graph_weights = graph_weights_u(from_to, n_cells);
      if(int_path) {
        const std::vector<int> starts = get_starts_i(starts_targets);
        const std::vector<int> targets = get_targets_i(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
        from_to["from"] = R_NilValue;
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, false, starting_indices,
          show_progress, bar, distances);
      } else {
        const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
        const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
        from_to["from"] = R_NilValue;
        dists_wweights(graph_to, graph_weights, n_cells, starts, targets, early_stopping, ncores, directed, pairwise, false, starting_indices,
          show_progress, bar, distances);
      }
    }
    return Rcpp::wrap(distances);
  }
}
