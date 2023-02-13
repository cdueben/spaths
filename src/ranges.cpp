#include <Rcpp.h>
// [[Rcpp::plugins(cpp17)]]

// Row ranges
// [[Rcpp::export]]
Rcpp::IntegerVector define_ranges(Rcpp::IntegerVector &n_rows, int avg_rows, int ncores) {
  Rcpp::IntegerVector rows_per_element(ncores);
  int i {0};
  int n_unique = n_rows.length();
  for(int nc {0}; nc < ncores; nc++) {
    int nr {0};
    while(nr < avg_rows && i < n_unique) {
      nr += n_rows[i];
      i++;
    }
    rows_per_element[nc] = nr;
  }
  if(i != n_unique) {
    int nr {0};
    while(i < n_unique) {
      nr += n_rows[i];
      i++;
    }
    rows_per_element[ncores - 1] += nr;
  }
  return rows_per_element;
}

// Origins vector
// [[Rcpp::export]]
Rcpp::IntegerVector list_origins(int n_origins) {
  Rcpp::IntegerVector origins(n_origins * (n_origins - 1) / 2);
  int i {0};
  for(int r {1}; r < n_origins; r++) {
    int j = r + 1;
    for(int o {1}; o < j; o++) {
      origins[i] = o;
      i++;
    }
  }
  return origins;
}
