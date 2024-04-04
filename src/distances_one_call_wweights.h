#ifndef DISTANCESONECALLWWEIGHTS_H
#define DISTANCESONECALLWWEIGHTS_H

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
#include "target_distances.h"

void dists_one_call_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<double> >& graph_weights, const std::size_t n_cells,
  const int start, const std::vector<int>& targets, const bool early_stopping, const bool show_progress, std::vector<double>& distances);
void dists_one_call_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<float> >& graph_weights, const std::size_t n_cells,
  const int start, const std::vector<int>& targets, const bool early_stopping, const bool show_progress, std::vector<float>& distances);
void dists_one_call_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<int> >& graph_weights, const std::size_t n_cells,
  const int start, const std::vector<int>& targets, const bool early_stopping, const bool show_progress, std::vector<int>& distances);
void dists_one_call_wweights(const std::vector<std::vector<int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const int start, const std::vector<int>& targets, const bool early_stopping, const bool show_progress,
  std::vector<unsigned short int>& distances);
void dists_one_call_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<double> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const bool show_progress, std::vector<double>& distances);
void dists_one_call_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<float> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const bool show_progress, std::vector<float>& distances);
void dists_one_call_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<int> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const bool show_progress, std::vector<int>& distances);
void dists_one_call_wweights(const std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const unsigned short int start, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const bool show_progress, std::vector<unsigned short int>& distances);

#endif
