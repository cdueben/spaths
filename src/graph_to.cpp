// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include "graph_to.h"

// adjacency list of vertices
// std::vector<std::vector<int> > graph_to_i
// std::vector<std::vector<unsigned short int> > graph_to_u

std::vector<std::vector<int> > graph_to_i(Rcpp::List& from_to, const std::size_t n_cells) {
  Rcpp::IntegerVector from = from_to["from"];
  Rcpp::IntegerVector to = from_to["to"];
  
  const std::size_t n_edges = from.size();
  std::vector<std::vector<int> > graph_to(n_cells);
  
  for(std::size_t i = 0; i < n_edges; ++i) {
    graph_to[from[i]].push_back(to[i]);
  }
  
  from_to["to"] = R_NilValue;
  return graph_to;
}

std::vector<std::vector<unsigned short int> > graph_to_u(Rcpp::List& from_to, const std::size_t n_cells) {
  Rcpp::IntegerVector from = from_to["from"];
  Rcpp::IntegerVector to = from_to["to"];
  
  const std::size_t n_edges = from.size();
  std::vector<std::vector<unsigned short int> > graph_to(n_cells);
  
  for(std::size_t i = 0; i < n_edges; ++i) {
    graph_to[from[i]].push_back(to[i]);
  }
  
  from_to["to"] = R_NilValue;
  return graph_to;
}
