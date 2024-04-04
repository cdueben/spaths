// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <cstddef>

// destination names, when destinations are not specified

// [[Rcpp::export]]
Rcpp::CharacterVector destination_names_character(Rcpp::CharacterVector& origin_names) {
  const std::size_t n_names = origin_names.size();
  Rcpp::CharacterVector destination_names (n_names * (n_names - 1) / 2);
  int position {0};
  const std::size_t n_names_1 = n_names - 1;
  for(std::size_t i = 0; i < n_names_1; ++i) {
    for(std::size_t j = i + 1; j < n_names; ++j) {
      destination_names[position] = origin_names[j];
      ++position;
    }
  }
  return destination_names;
}

// [[Rcpp::export]]
Rcpp::NumericVector destination_names_numeric(Rcpp::NumericVector& origin_names) {
  const std::size_t n_names = origin_names.size();
  Rcpp::NumericVector destination_names (n_names * (n_names - 1) / 2);
  int position {0};
  const std::size_t n_names_1 = n_names - 1;
  for(std::size_t i = 1; i < n_names_1; ++i) {
    for(std::size_t j = i + 1; j < n_names; ++j) {
      destination_names[position] = origin_names[j];
      ++position;
    }
  }
  return destination_names;
}

// [[Rcpp::export]]
Rcpp::IntegerVector destination_names_integer(Rcpp::IntegerVector& origin_names) {
  const std::size_t n_names = origin_names.size();
  Rcpp::IntegerVector destination_names (n_names * (n_names - 1) / 2);
  int position {0};
  const std::size_t n_names_1 = n_names - 1;
  for(std::size_t i = 1; i < n_names_1; ++i) {
    for(std::size_t j = i + 1; j < n_names; ++j) {
      destination_names[position] = origin_names[j];
      ++position;
    }
  }
  return destination_names;
}

// [[Rcpp::export]]
Rcpp::IntegerVector destination_names_auto(const int n_names) {
  Rcpp::IntegerVector destination_names (n_names * (n_names - 1) / 2);
  int position {0};
  const int n_names_1 = n_names + 1;
  for(int i = 1; i < n_names; ++i) {
    for(int j = i + 1; j < n_names_1; ++j) {
      destination_names[position] = j;
      ++position;
    }
  }
  return destination_names;
}
