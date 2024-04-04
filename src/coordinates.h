#ifndef COORDINATES_H
#define COORDINATES_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <cstddef>

inline std::vector<int> get_cell_numbers(Rcpp::List& coords, const bool rm_rvector = true) {
  Rcpp::IntegerVector cell_numbers_r = coords["cell_numbers"];
  const std::vector<int> cell_numbers_c = Rcpp::as<std::vector<int> >(cell_numbers_r);
  if(rm_rvector) {
    coords["cell_numbers"] = R_NilValue;
  }
  return cell_numbers_c;
}

Rcpp::List coordinates(const std::vector<std::vector<int> >& paths, Rcpp::List& coords, const int n_upd_rst, const bool return_dists,
  std::vector<int>& unconnected_indices);
Rcpp::List coordinates(std::vector<std::vector<unsigned short int> >& paths, Rcpp::List& coords, const int n_upd_rst, const bool return_dists,
  std::vector<int>& unconnected_indices);
void coordinates(Rcpp::List paths_coordinates, std::vector<std::vector<std::vector<int> > >& upd_paths, Rcpp::List& coords, const bool return_dists,
  std::vector<int>& unconnected_indices);
void coordinates(Rcpp::List paths_coordinates, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths, Rcpp::List& coords,
  const bool return_dists, std::vector<int>& unconnected_indices);

#endif
