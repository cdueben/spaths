// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <cstddef>
#include "starts_targets.h"
#include "convert_upd_rst.h"
#include "graph_to.h"
#include "graph_weights.h"
#include "upd_distances_wweights.h"

// distances with precomputed weights and grid updating
// distances that are exported to R
// Rcpp::NumericVector r_upd_dists_wweights_d
// Rcpp::IntegerVector r_upd_dists_wweights_i

// [[Rcpp::export]]
Rcpp::NumericVector r_upd_dists_wweights_d(Rcpp::List& from_to, Rcpp::List& starts_targets, const std::size_t n_cells, Rcpp::List& upd_rst_r,
  const bool early_stopping, const int ncores, const bool pairwise, const bool directed, const bool par_lvl_upd, const bool int_path,
  const bool double_weights, const bool show_progress, const int bar_limit) {
  
  const int n_paths = compute_n_paths(starts_targets, directed, pairwise);
  
  // with this structure, upd_rst_c, graph_to, and graph_weights go out of scope before the R object is assembled (alternative to .swap() etc.)
  if(double_weights) {
    std::vector<double> distances (n_paths);
    {
      std::vector<std::vector<double> > graph_weights = graph_weights_d(from_to, n_cells);
      if(int_path) {
        const std::vector<int> starts = get_starts_i(starts_targets);
        const std::vector<int> targets = get_targets_i(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        const std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
        const std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      } else {
        const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
        const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        const std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
        const std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
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
        const std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
        const std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      } else {
        const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
        const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        const std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
        const std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      }
    }
    return Rcpp::wrap(distances);
  }
}

// [[Rcpp::export]]
Rcpp::IntegerVector r_upd_dists_wweights_i(Rcpp::List& from_to, Rcpp::List& starts_targets, const std::size_t n_cells, Rcpp::List& upd_rst_r,
  const bool early_stopping, const int ncores, const bool pairwise, const bool directed, const bool par_lvl_upd, const bool int_path,
  const bool signed_weights, const bool show_progress, const int bar_limit) {
  
  const int n_paths = compute_n_paths(starts_targets, directed, pairwise);
  
  if(signed_weights) {
    std::vector<int> distances (n_paths);
    {
      std::vector<std::vector<int> > graph_weights = graph_weights_i(from_to, n_cells);
      if(int_path) {
        const std::vector<int> starts = get_starts_i(starts_targets);
        const std::vector<int> targets = get_targets_i(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        const std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
        const std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      } else {
        const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
        const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        const std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
        const std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
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
        const std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
        const std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      } else {
        const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
        const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
        const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
        const std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
        const std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
        upd_dists_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
          upd_rst_c, show_progress, bar_limit, distances);
      }
    }
    return Rcpp::wrap(distances);
  }
}
