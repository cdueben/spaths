// This file makes RcppExports.cpp aware of some of the headers used in the package
// Without it, the package does not compile because Rcpp is not aware of, e.g., the From_To structs

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include "../../src/structs.h"
