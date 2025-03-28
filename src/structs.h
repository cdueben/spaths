#ifndef STRUCTS_H
#define STRUCTS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>

struct From_To_I {
  std::vector<int> from;
  std::vector<int> to;
};

struct From_To_U {
  std::vector<unsigned short int> from;
  std::vector<unsigned short int> to;
};

#endif
