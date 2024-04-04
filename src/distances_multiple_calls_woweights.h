#ifndef DISTANCESMULTIPLECALLSWOWEIGHTS_H
#define DISTANCESMULTIPLECALLSWOWEIGHTS_H

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
#include "individual_distances.h"
#include "visited.h"

void dists_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths, std::vector<double>& distances,
  const double radius2, const int starting_index = -1, const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths, std::vector<float>& distances,
  const double radius2, const int starting_index = -1, const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths, std::vector<int>& distances,
  const double radius2, const int starting_index = -1, const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping,
  const bool haversine, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<unsigned short int>& distances, const double radius2, const int starting_index = -1, const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers,
  const int ncol, const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<int>& affected_paths, std::vector<double>& distances, const double radius2, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<int>& affected_paths, std::vector<float>& distances, const double radius2, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers,
  const int ncol, const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<int>& affected_paths, std::vector<int>& distances, const double radius2, const int starting_index = -1, const int n_targets = -1,
  const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<int>& affected_paths, std::vector<unsigned short int>& distances, const double radius2, const int starting_index = -1,
  const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths, std::vector<double>& distances,
  const double radius2, const int starting_index = -1, const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths, std::vector<float>& distances,
  const double radius2, const int starting_index = -1, const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths, std::vector<int>& distances,
  const double radius2, const int starting_index = -1, const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping,
  const bool haversine, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths,
  std::vector<unsigned short int>& distances, const double radius2, const int starting_index = -1, const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers,
  const int ncol, const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<double>& distances, const double radius2, const int starting_index = -1,
  const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<float>& distances, const double radius2, const int starting_index = -1,
  const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers,
  const int ncol, const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<int>& distances, const double radius2, const int starting_index = -1,
  const int n_targets = -1, const int begin_target = -1);
void dists_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<unsigned short int>& distances, const double radius2, const int starting_index = -1,
  const int n_targets = -1, const int begin_target = -1);

#endif
