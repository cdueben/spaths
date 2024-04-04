#ifndef CONVERTUPDRST_H
#define CONVERTUPDRST_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_set>
#include <cstddef>

std::vector<std::unordered_set<int> > convert_upd_rst_i(Rcpp::List& upd_rst_r);
std::vector<std::unordered_set<unsigned short int> > convert_upd_rst_u(Rcpp::List& upd_rst_r);

#endif
