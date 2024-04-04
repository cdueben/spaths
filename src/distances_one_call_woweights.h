#ifndef DISTANCESONECALLWOWEIGHTS_H
#define DISTANCESONECALLWOWEIGHTS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <unordered_set>
#include <queue>
#include <utility>
#include <functional>
#include "pair_types.h"
#include "target_distances.h"
#include "individual_distances.h"

void dists_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const bool show_progress, std::vector<double>& distances, const double radius2);
void dists_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const bool show_progress, std::vector<float>& distances, const double radius2);
void dists_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const bool show_progress, std::vector<int>& distances, const double radius2);
void dists_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const bool show_progress, std::vector<unsigned short int>& distances, const double radius2);
void dists_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const bool show_progress, std::vector<double>& distances, const double radius2);
void dists_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const bool show_progress, std::vector<float>& distances, const double radius2);
void dists_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const bool show_progress, std::vector<int>& distances, const double radius2);
void dists_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const bool show_progress, std::vector<unsigned short int>& distances, const double radius2);

#endif
