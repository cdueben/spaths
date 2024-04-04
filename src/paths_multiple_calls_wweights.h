#ifndef PATHSMULTIPLECALLSWWEIGHTS_H
#define PATHSMULTIPLECALLSWWEIGHTS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <unordered_set>
#include <queue>
#include <utility>
#include <functional>
#include <cstddef>
#include "pair_types.h"
#include "targets_set.h"
#include "target_distances.h"
#include "target_paths.h"
#include "visited.h"

void paths_multiple_calls_wweights(const std::vector<std::vector<int> >& graph_to, const std::vector<std::vector<double> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<int> >& paths, std::vector<double>& distances, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<int> >& graph_to, const std::vector<std::vector<float> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<int> >& paths, std::vector<float>& distances, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<int> >& graph_to, const std::vector<std::vector<int> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<int> >& paths, std::vector<int>& distances, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<int> >& graph_to, const std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<int> >& paths, std::vector<unsigned short int>& distances, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<std::vector<double> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<int>& affected_paths, std::vector<std::vector<unsigned short int> >& paths, std::vector<double>& distances,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<std::vector<float> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<int>& affected_paths, std::vector<std::vector<unsigned short int> >& paths, std::vector<float>& distances,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<std::vector<int> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<int>& affected_paths, std::vector<std::vector<unsigned short int> >& paths, std::vector<int>& distances, const int starting_index = -1,
  const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<unsigned short int> >& graph_to,
  const std::vector<std::vector<unsigned short int> >& graph_weights, const std::size_t n_cells, const unsigned short int start,
  const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<unsigned short int>& distances, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<int> >& graph_to, const std::vector<std::vector<double> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths, std::vector<double>& distances,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<int> >& graph_to, const std::vector<std::vector<float> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths, std::vector<float>& distances,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<int> >& graph_to, const std::vector<std::vector<int> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths, std::vector<int>& distances, const int starting_index = -1,
  const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<int> >& graph_to, const std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths, std::vector<unsigned short int>& distances,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<std::vector<double> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<unsigned short int> >& paths, std::vector<double>& distances,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<std::vector<float> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<unsigned short int> >& paths, std::vector<float>& distances,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<std::vector<int> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<unsigned short int> >& paths, std::vector<int>& distances,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1, const int exclude_index = -1);
void paths_multiple_calls_wweights(const std::vector<std::vector<unsigned short int> >& graph_to,
  const std::vector<std::vector<unsigned short int> >& graph_weights, const std::size_t n_cells, const unsigned short int start,
  const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<unsigned short int>& distances, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1, const int exclude_index = -1);

#endif
