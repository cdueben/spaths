#ifndef GRAPHTO_H
#define GRAPHTO_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>

std::vector<std::vector<int> > graph_to_i(Rcpp::List& from_to, const std::size_t n_cells);
std::vector<std::vector<unsigned short int> > graph_to_u(Rcpp::List& from_to, const std::size_t n_cells);

#endif
