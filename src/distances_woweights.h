#ifndef DISTANCESWOWEIGHTS_H
#define DISTANCESWOWEIGHTS_H

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
#include "distances_one_call_woweights.h"
#include "distances_multiple_calls_woweights.h"
#include "show_progress.h"

void dists_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<double>& distances,
  const double radius2);
void dists_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<float>& distances,
  const double radius2);
void dists_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<int>& distances,
  const double radius2);
void dists_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres, const double yres,
  const double ymax, const std::vector<int>& starts, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const bool pairwise, const std::vector<int>& starting_indices, const bool show_progress, const bool bar, std::vector<unsigned short int>& distances,
  const double radius2);
void dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<double>& distances, const double radius2);
void dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<float>& distances, const double radius2);
void dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<int>& distances, const double radius2);
void dists_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const bool pairwise, const std::vector<int>& starting_indices,
  const bool show_progress, const bool bar, std::vector<unsigned short int>& distances, const double radius2);

#endif
