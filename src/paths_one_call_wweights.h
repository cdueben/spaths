#ifndef PATHSONECALLWWEIGHTS_H
#define PATHSONECALLWWEIGHTS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <utility>
#include <functional>
#include <cstddef>
#include "pair_types.h"
#include "target_distances.h"
#include "target_paths.h"

void paths_one_call_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<double> >& graph_weights, const std::size_t n_cells,
  const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores, const std::unordered_set<int>& graph_to_0,
  const bool show_progress, std::vector<std::vector<int> >& paths, std::vector<double>& distances);
void paths_one_call_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<float> >& graph_weights, const std::size_t n_cells,
  const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores, const std::unordered_set<int>& graph_to_0,
  const bool show_progress, std::vector<std::vector<int> >& paths, std::vector<float>& distances);
void paths_one_call_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<int> >& graph_weights, const std::size_t n_cells,
  const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores, const std::unordered_set<int>& graph_to_0,
  const bool show_progress, std::vector<std::vector<int> >& paths, std::vector<int>& distances);
void paths_one_call_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<int>& graph_to_0, const bool show_progress, std::vector<std::vector<int> >& paths, std::vector<unsigned short int>& distances);
void paths_one_call_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<double> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, std::vector<std::vector<unsigned short int> >& paths,
  std::vector<double>& distances);
void paths_one_call_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<float> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, std::vector<std::vector<unsigned short int> >& paths,
  std::vector<float>& distances);
void paths_one_call_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<int> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, std::vector<std::vector<unsigned short int> >& paths,
  std::vector<int>& distances);
void paths_one_call_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress, std::vector<std::vector<unsigned short int> >& paths,
  std::vector<unsigned short int>& distances);

#endif
