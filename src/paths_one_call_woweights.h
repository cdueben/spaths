#ifndef PATHSONECALLWOWEIGHTS_H
#define PATHSONECALLWOWEIGHTS_H

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
#include "individual_distances.h"
#include "target_distances.h"
#include "target_paths.h"

void paths_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const bool show_progress, std::vector<std::vector<int> >& paths, std::vector<double>& distances,
  const double radius2);
void paths_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const bool show_progress, std::vector<std::vector<int> >& paths, std::vector<float>& distances,
  const double radius2);
void paths_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const bool show_progress, std::vector<std::vector<int> >& paths, std::vector<int>& distances,
  const double radius2);
void paths_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const bool show_progress, std::vector<std::vector<int> >& paths, std::vector<unsigned short int>& distances,
  const double radius2);
void paths_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const int start, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<double>& distances, const double radius2);
void paths_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const int start, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<float>& distances, const double radius2);
void paths_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const int start, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<int>& distances, const double radius2);
void paths_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const int start, const std::vector<unsigned short int>& targets, const bool early_stopping,
  const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0, const bool show_progress,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<unsigned short int>& distances, const double radius2);

#endif
