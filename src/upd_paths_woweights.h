#ifndef UPDPATHSWOWEIGHTS_H
#define UPDPATHSWOWEIGHTS_H

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
#include "paths_woweights.h"
#include "upd_affected_paths.h"
#include "upd_starts_targets_map.h"
#include "paths_multiple_calls_woweights.h"
#include "repeat_distances.h"
#include "show_progress.h"
#include "coordinates.h"

void upd_paths_woweights(std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const std::vector<int>& starts,
  const std::vector<int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine,
  const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress, const int bar_limit,
  std::vector<std::vector<int> >& static_paths, std::vector<std::vector<std::vector<int> > >& upd_paths, std::vector<double>& distances,
  const double radius2);
void upd_paths_woweights(std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const std::vector<int>& starts,
  const std::vector<int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine,
  const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress, const int bar_limit,
  std::vector<std::vector<int> >& static_paths, std::vector<std::vector<std::vector<int> > >& upd_paths, std::vector<float>& distances,
  const double radius2);
void upd_paths_woweights(std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const std::vector<int>& starts,
  const std::vector<int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine,
  const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress, const int bar_limit,
  std::vector<std::vector<int> >& static_paths, std::vector<std::vector<std::vector<int> > >& upd_paths, std::vector<int>& distances, const double radius2);
void upd_paths_woweights(std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const std::vector<int>& starts,
  const std::vector<int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping, const bool haversine,
  const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<int> >& upd_rst, const bool show_progress, const int bar_limit,
  std::vector<std::vector<int> >& static_paths, std::vector<std::vector<std::vector<int> > >& upd_paths, std::vector<unsigned short int>& distances,
  const double radius2);
void upd_paths_woweights(std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping,
  const bool haversine, const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<unsigned short int> >& upd_rst, const bool show_progress,
  const int bar_limit, std::vector<std::vector<unsigned short int> >& static_paths, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths,
  std::vector<double>& distances, const double radius2);
void upd_paths_woweights(std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping,
  const bool haversine, const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<unsigned short int> >& upd_rst, const bool show_progress,
  const int bar_limit, std::vector<std::vector<unsigned short int> >& static_paths, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths,
  std::vector<float>& distances, const double radius2);
void upd_paths_woweights(std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping,
  const bool haversine, const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<unsigned short int> >& upd_rst, const bool show_progress,
  const int bar_limit, std::vector<std::vector<unsigned short int> >& static_paths, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths,
  std::vector<int>& distances, const double radius2);
void upd_paths_woweights(std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const std::vector<int>& starting_indices, const bool pairwise, const bool early_stopping,
  const bool haversine, const int ncores, const bool par_lvl_upd, std::vector<std::unordered_set<unsigned short int> >& upd_rst, const bool show_progress,
  const int bar_limit, std::vector<std::vector<unsigned short int> >& static_paths, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths,
  std::vector<unsigned short int>& distances, const double radius2);

#endif
