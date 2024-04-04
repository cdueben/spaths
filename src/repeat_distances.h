#ifndef REPEATDISTANCES_H
#define REPEATDISTANCES_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <iterator>

// repeat the static distances
void repeat_distances(std::vector<double>& distances, const int n_upd_rst);
void repeat_distances(std::vector<float>& distances, const int n_upd_rst);
void repeat_distances(std::vector<int>& distances, const int n_upd_rst);
void repeat_distances(std::vector<unsigned short int>& distances, const int n_upd_rst);

#endif
