// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include "repeat_distances.h"

// repeat static distances
// void repeat_distances
// void repeat_distances
// void repeat_distances
// void repeat_distances

void repeat_distances(std::vector<double>& distances, const int n_upd_rst) {
  const int n_paths = distances.size();
  if(n_paths != 0) {
    distances.reserve(n_paths * (n_upd_rst + 1));
    for(int u = 0; u < n_upd_rst; ++u) {
      std::copy_n(distances.begin(), n_paths, std::back_inserter(distances));
    }
  }
}

void repeat_distances(std::vector<float>& distances, const int n_upd_rst) {
  const int n_paths = distances.size();
  if(n_paths != 0) {
    distances.reserve(n_paths * (n_upd_rst + 1));
    for(int u = 0; u < n_upd_rst; ++u) {
      std::copy_n(distances.begin(), n_paths, std::back_inserter(distances));
    }
  }
}

void repeat_distances(std::vector<int>& distances, const int n_upd_rst) {
  const int n_paths = distances.size();
  if(n_paths != 0) {
    distances.reserve(n_paths * (n_upd_rst + 1));
    for(int u = 0; u < n_upd_rst; ++u) {
      std::copy_n(distances.begin(), n_paths, std::back_inserter(distances));
    }
  }
}

void repeat_distances(std::vector<unsigned short int>& distances, const int n_upd_rst) {
  const int n_paths = distances.size();
  if(n_paths != 0) {
    distances.reserve(n_paths * (n_upd_rst + 1));
    for(int u = 0; u < n_upd_rst; ++u) {
      std::copy_n(distances.begin(), n_paths, std::back_inserter(distances));
    }
  }
}
