#ifndef DISTANCESWWEIGHTS_H
#define DISTANCESWWEIGHTS_H

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
#include <iostream>
#include "distances_one_call_wweights.h"
#include "distances_multiple_calls_wweights.h"
#include "show_progress.h"

void dists_wweights(std::vector<std::vector<int> >& graph_to, std::vector<std::vector<double> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<double>& distances);
void dists_wweights(std::vector<std::vector<int> >& graph_to, std::vector<std::vector<float> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<float>& distances);
void dists_wweights(std::vector<std::vector<int> >& graph_to, std::vector<std::vector<int> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<int>& distances);
void dists_wweights(std::vector<std::vector<int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights, const std::size_t n_cells,
  const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const int ncores, const bool directed, const bool pairwise,
  const bool upd_rst_defined, const std::vector<int>& starting_indices, const bool show_progress, const bool bar,
  std::vector<unsigned short int>& distances);
void dists_wweights(std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<double> >& graph_weights, const std::size_t n_cells,
  const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const bool directed, const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices, const bool show_progress, const bool bar,
  std::vector<double>& distances);
void dists_wweights(std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<float> >& graph_weights, const std::size_t n_cells,
  const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const bool directed, const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices, const bool show_progress, const bool bar,
  std::vector<float>& distances);
void dists_wweights(std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<int> >& graph_weights, const std::size_t n_cells,
  const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, const bool early_stopping, const int ncores,
  const bool directed, const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices, const bool show_progress, const bool bar,
  std::vector<int>& distances);
void dists_wweights(std::vector<std::vector<unsigned short int> >& graph_to, std::vector<std::vector<unsigned short int> >& graph_weights,
  const std::size_t n_cells, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const int ncores, const bool directed, const bool pairwise, const bool upd_rst_defined, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<unsigned short int>& distances);

#endif
