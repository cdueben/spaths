#ifndef GRAPHWEIGHTS_H
#define GRAPHWEIGHTS_H

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
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstddef>
#include "coordinates.h"
#include "structs.h"

std::vector<std::vector<double> > graph_weights_d(Rcpp::List& from_to, const std::size_t n_cells, const bool int_path, const bool from_to_r);
std::vector<std::vector<float> > graph_weights_f(Rcpp::List& from_to, const std::size_t n_cells, const bool int_path, const bool from_to_r);
std::vector<std::vector<int> > graph_weights_i(Rcpp::List& from_to, const std::size_t n_cells, const bool int_path, const bool from_to_r);
std::vector<std::vector<unsigned short int> > graph_weights_u(Rcpp::List& from_to, const std::size_t n_cells, const bool int_path, const bool from_to_r);
std::vector<std::vector<double> > graph_weights_d(const std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2);
std::vector<std::vector<float> > graph_weights_f(const std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const bool haversine, const bool queen,
  const int ncores, const bool rm_rvector, const double radius2);
std::vector<std::vector<int> > graph_weights_i(const std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const bool haversine, const bool queen,
  const int ncores, const bool rm_rvector, const double radius2);
std::vector<std::vector<unsigned short int> > graph_weights_u(const std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2);
std::vector<std::vector<double> > graph_weights_d(const std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2);
std::vector<std::vector<float> > graph_weights_f(const std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2);
std::vector<std::vector<int> > graph_weights_i(const std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2);
std::vector<std::vector<unsigned short int> > graph_weights_u(const std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords,
  const bool haversine, const bool queen, const int ncores, const bool rm_rvector, const double radius2);

#endif
