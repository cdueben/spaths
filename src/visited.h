#ifndef VISITED_H
#define VISITED_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <cstddef>

// initialize vector of visited cells
// functions are overloaded with int and unsigned short int upd_rst
template <typename U> // U: upd_rst type
std::vector<bool> create_visited(const std::size_t n_cells, const std::unordered_set<U>& upd_rst) {
  std::vector<bool> visited (n_cells);
  if(!upd_rst.empty()) {
    for(const auto& i: upd_rst) {
      visited[i] = true;
    }
  }
  return visited;
}

#endif
