// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <cstddef>
#include "starts_targets.h"
#include "convert_upd_rst.h"
#include "graph_to.h"
#include "graph_weights.h"
#include "upd_paths_wweights.h"
#include "coordinates.h"

// paths with precomputed weights and grid updating
// paths that are exported to R
// Rcpp::List r_upd_paths_wweights

// [[Rcpp::export]]
Rcpp::List r_upd_paths_wweights(Rcpp::List& from_to, Rcpp::List& starts_targets, Rcpp::List& coords, const std::size_t n_cells, Rcpp::List& upd_rst_r,
  const bool early_stopping, const int ncores, const bool pairwise, const bool directed, const bool par_lvl_upd, const bool int_path,
  const bool numeric_weights, const bool double_weights, const bool signed_weights, const bool return_dists, const bool show_progress,
  const int bar_limit) {
  // from_to contains a from vector, a to vector, and a weights vector (numeric or integer)
  // coords contains a cell_numbers vector, a double_coords boolean, a ncol integer, and xmin, ymax, xres, and yres doubles or integers
  // starts_targets contains a starts vector, a targets vector, a n_paths_per_start vector, a n_starts integer, and a n_targets integer
  // directed denotes whether the graph is directed
  // par_lvl_upd denotes whether the code parallelizes at the upd_rst level
  
  const int n_upd_rst = upd_rst_r.size();
  const int n_paths = compute_n_paths(starts_targets, directed, pairwise);
  const int n_dists = (return_dists) ? n_paths : 0;
  
  Rcpp::List paths_r = Rcpp::List::create(Rcpp::Named("paths") = R_NilValue);
  
  std::vector<int> unconnected_indices;
  
  // with this structure, upd_rst_c, graph_to, and graph_weights go out of scope before the R object data is assembled (alternative to .swap() etc.)
  if(int_path) {
    std::vector<std::vector<std::vector<int> > > upd_paths(n_upd_rst, std::vector<std::vector<int> >(n_paths));
    {
      std::vector<std::vector<int> > static_paths(n_paths);
      if(numeric_weights) {
        if(double_weights) {
          std::vector<double> distances (n_dists);
          {
            const std::vector<int> starts = get_starts_i(starts_targets);
            const std::vector<int> targets = get_targets_i(starts_targets);
            const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
            std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
            const std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
            std::vector<std::vector<double> > graph_weights = graph_weights_d(from_to, n_cells);
            upd_paths_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
              upd_rst_c, show_progress, bar_limit, static_paths, upd_paths, distances);
          }
          if(return_dists) {
            paths_r["distances"] = Rcpp::wrap(distances);
          }
        } else {
          std::vector<float> distances (n_dists);
          {
            const std::vector<int> starts = get_starts_i(starts_targets);
            const std::vector<int> targets = get_targets_i(starts_targets);
            const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
            std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
            const std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
            std::vector<std::vector<float> > graph_weights = graph_weights_f(from_to, n_cells);
            upd_paths_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
              upd_rst_c, show_progress, bar_limit, static_paths, upd_paths, distances);
          }
          if(return_dists) {
            paths_r["distances"] = Rcpp::wrap(distances);
          }
        }
      } else {
        if(signed_weights) {
          std::vector<int> distances (n_dists);
          {
            const std::vector<int> starts = get_starts_i(starts_targets);
            const std::vector<int> targets = get_targets_i(starts_targets);
            const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
            std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
            const std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
            std::vector<std::vector<int> > graph_weights = graph_weights_i(from_to, n_cells);
            upd_paths_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
              upd_rst_c, show_progress, bar_limit, static_paths, upd_paths, distances);
          }
          if(return_dists) {
            paths_r["distances"] = Rcpp::wrap(distances);
          }
        } else {
          std::vector<unsigned short int> distances (n_dists);
          {
            const std::vector<int> starts = get_starts_i(starts_targets);
            const std::vector<int> targets = get_targets_i(starts_targets);
            const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
            std::vector<std::unordered_set<int> > upd_rst_c = convert_upd_rst_i(upd_rst_r);
            const std::vector<std::vector<int> > graph_to = graph_to_i(from_to, n_cells);
            std::vector<std::vector<unsigned short int> > graph_weights = graph_weights_u(from_to, n_cells);
            upd_paths_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
              upd_rst_c, show_progress, bar_limit, static_paths, upd_paths, distances);
          }
          if(return_dists) {
            paths_r["distances"] = Rcpp::wrap(distances);
          }
        }
      }
      paths_r["paths"] = coordinates(static_paths, coords, n_upd_rst, return_dists, unconnected_indices);
    }
    coordinates(paths_r["paths"], upd_paths, coords, return_dists, unconnected_indices);
  } else {
    std::vector<std::vector<std::vector<unsigned short int> > > upd_paths(n_upd_rst, std::vector<std::vector<unsigned short int> >(n_paths));
    {
      std::vector<std::vector<unsigned short int> > static_paths(n_paths);
      if(numeric_weights) {
        if(double_weights) {
          std::vector<double> distances (n_dists);
          {
            const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
            const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
            const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
            std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
            const std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
            std::vector<std::vector<double> > graph_weights = graph_weights_d(from_to, n_cells);
            upd_paths_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
              upd_rst_c, show_progress, bar_limit, static_paths, upd_paths, distances);
          }
          if(return_dists) {
            paths_r["distances"] = Rcpp::wrap(distances);
          }
        } else {
          std::vector<float> distances (n_dists);
          {
            const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
            const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
            const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
            std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
            const std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
            std::vector<std::vector<float> > graph_weights = graph_weights_f(from_to, n_cells);
            upd_paths_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
              upd_rst_c, show_progress, bar_limit, static_paths, upd_paths, distances);
          }
          if(return_dists) {
            paths_r["distances"] = Rcpp::wrap(distances);
          }
        }
      } else {
        if(signed_weights) {
          std::vector<int> distances (n_dists);
          {
            const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
            const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
            const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
            std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
            const std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
            std::vector<std::vector<int> > graph_weights = graph_weights_i(from_to, n_cells);
            upd_paths_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
              upd_rst_c, show_progress, bar_limit, static_paths, upd_paths, distances);
          }
          if(return_dists) {
            paths_r["distances"] = Rcpp::wrap(distances);
          }
        } else {
          std::vector<unsigned short int> distances (n_dists);
          {
            const std::vector<unsigned short int> starts = get_starts_u(starts_targets);
            const std::vector<unsigned short int> targets = get_targets_u(starts_targets);
            const std::vector<int> starting_indices = get_starting_indices_i(starts_targets, (int) starts.size(), (targets.empty() && !directed), pairwise);
            std::vector<std::unordered_set<unsigned short int> > upd_rst_c = convert_upd_rst_u(upd_rst_r);
            const std::vector<std::vector<unsigned short int> > graph_to = graph_to_u(from_to, n_cells);
            std::vector<std::vector<unsigned short int> > graph_weights = graph_weights_u(from_to, n_cells);
            upd_paths_wweights(graph_to, graph_weights, n_cells, starts, targets, starting_indices, pairwise, directed, early_stopping, ncores, par_lvl_upd,
              upd_rst_c, show_progress, bar_limit, static_paths, upd_paths, distances);
          }
          if(return_dists) {
            paths_r["distances"] = Rcpp::wrap(distances);
          }
        }
      }
      paths_r["paths"] = coordinates(static_paths, coords, n_upd_rst, return_dists, unconnected_indices);
    }
    coordinates(paths_r["paths"], upd_paths, coords, return_dists, unconnected_indices);
  }
  if(!return_dists) {
    paths_r["unconnected_indices"] = Rcpp::wrap(unconnected_indices);
  }
  
  return paths_r;
}
