// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <cstddef>
#include "convert_upd_rst.h"

// convert upd_rst to unordered_sets, which are faster to search for C++ than R vectors are
// std::vector<std::unordered_set<int> > convert_upd_rst_i
// std::vector<std::unordered_set<unsigned short int> > convert_upd_rst_u

std::vector<std::unordered_set<int> > convert_upd_rst_i(Rcpp::List& upd_rst_r) {
  const std::size_t n_upd_rst = upd_rst_r.size();
  std::vector<std::unordered_set<int> > upd_rst_c(n_upd_rst);
  
  for(std::size_t u = 0; u < n_upd_rst; ++u) {
    Rcpp::IntegerVector upd_rst_u = upd_rst_r[u];
    upd_rst_c[u].insert(upd_rst_u.begin(), upd_rst_u.end());
    upd_rst_r[u] = R_NilValue;
  }
  
  return upd_rst_c;
}

std::vector<std::unordered_set<unsigned short int> > convert_upd_rst_u(Rcpp::List& upd_rst_r) {
  const std::size_t n_upd_rst = upd_rst_r.size();
  std::vector<std::unordered_set<unsigned short int> > upd_rst_c(n_upd_rst);
  
  for(std::size_t u = 0; u < n_upd_rst; ++u) {
    Rcpp::IntegerVector upd_rst_u = upd_rst_r[u];
    upd_rst_c[u].insert(upd_rst_u.begin(), upd_rst_u.end());
    upd_rst_r[u] = R_NilValue;
  }
  
  return upd_rst_c;
}
