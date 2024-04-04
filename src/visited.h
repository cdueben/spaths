#ifndef VISITED_H
#define VISITED_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <cstddef>

std::vector<bool> create_visited(const std::size_t n_cells, const std::unordered_set<int>& upd_rst);
std::vector<bool> create_visited(const std::size_t n_cells, const std::unordered_set<unsigned short int>& upd_rst);

#endif
