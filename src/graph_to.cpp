// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include <type_traits>
#include "structs.h"
#include "graph_to.h"

// adjacency list of vertices
// std::vector<std::vector<int> > graph_to_i
// std::vector<std::vector<unsigned short int> > graph_to_u

template <typename T>
std::vector<std::vector<T> > graph_t(Rcpp::List& from_to, const std::size_t n_cells, const T t, const bool from_to_r) {
  
  std::vector<std::vector<T> > graph_to(n_cells);
  
  if(from_to_r) {
    Rcpp::IntegerVector from = from_to["from"];
    Rcpp::IntegerVector to = from_to["to"];
    
    const std::size_t n_edges = from.size();
    
    for(std::size_t i = 0; i < n_edges; ++i) {
      graph_to[from[i]].push_back(to[i]);
    }
    
    from_to["to"] = R_NilValue;
  } else {
    if constexpr (std::is_same_v<T, int>) {
      Rcpp::XPtr<From_To_I> ft = from_to["from_to"];
      
      const std::size_t n_edges = ft->from.size();
      
      for(std::size_t i = 0; i < n_edges; ++i) {
        graph_to[ft->from[i]].push_back(ft->to[i]);
      }
      
      std::vector<T>().swap(ft->to);
    } else {
      Rcpp::XPtr<From_To_U> ft = from_to["from_to"];
      
      const std::size_t n_edges = ft->from.size();
      
      for(std::size_t i = 0; i < n_edges; ++i) {
        graph_to[ft->from[i]].push_back(ft->to[i]);
      }
      
      std::vector<T>().swap(ft->to);
    }
  }
  
  return graph_to;
}

std::vector<std::vector<int> > graph_to_i(Rcpp::List& from_to, const std::size_t n_cells, const bool from_to_r) {
  constexpr int t {};
  return graph_t(from_to, n_cells, t, from_to_r);
}

std::vector<std::vector<unsigned short int> > graph_to_u(Rcpp::List& from_to, const std::size_t n_cells, const bool from_to_r) {
  constexpr unsigned short int t {};
  return graph_t(from_to, n_cells, t, from_to_r);
}
