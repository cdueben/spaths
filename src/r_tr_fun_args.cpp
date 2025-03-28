// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads()  1
#define omp_get_thread_num()   0
#define omp_get_max_threads()  1
#define omp_get_thread_limit() 1
#define omp_get_num_procs()    1
#endif
#include <vector>
#include <cstddef>
#include <string>
#include <type_traits>
#include "structs.h"

// this file assembles the inputs for transition functions receiving pointers instead if native R objects

template <typename T>
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_haversine_queen(Rcpp::XPtr<T> from_to, Rcpp::IntegerVector& row_number, Rcpp::IntegerVector& col_number,
  Rcpp::NumericVector& d_horizontal, const double d_vertical, Rcpp::NumericVector& d_diagonal, const int ncores) {
  
  const std::size_t nedges = from_to->from.size();
  std::vector<double>* distance = new std::vector<double>(nedges);
  
  #pragma omp simd
  for(std::size_t i = 0; i < nedges; ++i) {
    if(row_number[from_to->from[i]] == row_number[from_to->to[i]]) {
      (*distance)[i] = d_horizontal[row_number[from_to->from[i]]];
    } else if(col_number[from_to->from[i]] == col_number[from_to->to[i]]) {
      (*distance)[i] = d_vertical;
    } else {
      (*distance)[i] = d_diagonal[row_number[from_to->from[i]]];
    }
  }
  
  Rcpp::XPtr<std::vector<double> > p(distance);
  return p;
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_haversine_queen_i(Rcpp::XPtr<From_To_I> from_to, Rcpp::IntegerVector& row_number,
  Rcpp::IntegerVector& col_number, Rcpp::NumericVector& d_horizontal, const double d_vertical, Rcpp::NumericVector& d_diagonal, const int ncores) {
  
  return tr_fun_args_d_haversine_queen(from_to, row_number, col_number, d_horizontal, d_vertical, d_diagonal, ncores);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_haversine_queen_u(Rcpp::XPtr<From_To_U> from_to, Rcpp::IntegerVector& row_number,
  Rcpp::IntegerVector& col_number, Rcpp::NumericVector& d_horizontal, const double d_vertical, Rcpp::NumericVector& d_diagonal, const int ncores) {
  
  return tr_fun_args_d_haversine_queen(from_to, row_number, col_number, d_horizontal, d_vertical, d_diagonal, ncores);
}

template <typename T>
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_haversine_rook(Rcpp::XPtr<T> from_to, Rcpp::IntegerVector& row_number, Rcpp::NumericVector& d_horizontal,
  const double d_vertical, const int ncores) {
  
  const std::size_t nedges = from_to->from.size();
  std::vector<double>* distance = new std::vector<double>(nedges);
  
  #pragma omp simd
  for(std::size_t i = 0; i < nedges; ++i) {
    if(row_number[from_to->from[i]] == row_number[from_to->to[i]]) {
      (*distance)[i] = d_horizontal[row_number[from_to->from[i]]];
    } else {
      (*distance)[i] = d_vertical;
    }
  }
  
  Rcpp::XPtr<std::vector<double> > p(distance);
  return p;
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_haversine_rook_i(Rcpp::XPtr<From_To_I> from_to, Rcpp::IntegerVector& row_number,
  Rcpp::NumericVector& d_horizontal, const double d_vertical, const int ncores) {
  
  return tr_fun_args_d_haversine_rook(from_to, row_number, d_horizontal, d_vertical, ncores);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_haversine_rook_u(Rcpp::XPtr<From_To_U> from_to, Rcpp::IntegerVector& row_number,
  Rcpp::NumericVector& d_horizontal, const double d_vertical, const int ncores) {
  
  return tr_fun_args_d_haversine_rook(from_to, row_number, d_horizontal, d_vertical, ncores);
}

template <typename T>
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_euclidean_queen(Rcpp::XPtr<T> from_to, Rcpp::IntegerVector& row_number,
  Rcpp::IntegerVector& col_number, const double rst_xres, const double rst_yres, const double d_vertical, const double d_diagonal, const int ncores) {
  
  const std::size_t nedges = from_to->from.size();
  std::vector<double>* distance = new std::vector<double>(nedges);
  
  #pragma omp simd
  for(std::size_t i = 0; i < nedges; ++i) {
    if(row_number[from_to->from[i]] == row_number[from_to->to[i]]) {
      (*distance)[i] = rst_xres;
    } else if(col_number[from_to->from[i]] == col_number[from_to->to[i]]) {
      (*distance)[i] = rst_yres;
    } else {
      (*distance)[i] = d_diagonal;
    }
  }
  
  Rcpp::XPtr<std::vector<double> > p(distance);
  return p;
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_euclidean_queen_i(Rcpp::XPtr<From_To_I> from_to, Rcpp::IntegerVector& row_number,
  Rcpp::IntegerVector& col_number, const double rst_xres, const double rst_yres, const double d_vertical, const double d_diagonal, const int ncores) {
  
  return tr_fun_args_d_euclidean_queen(from_to, row_number, col_number, rst_xres, rst_yres, d_vertical, d_diagonal, ncores);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_euclidean_queen_u(Rcpp::XPtr<From_To_U> from_to, Rcpp::IntegerVector& row_number,
  Rcpp::IntegerVector& col_number, const double rst_xres, const double rst_yres, const double d_vertical, const double d_diagonal, const int ncores) {
  
  return tr_fun_args_d_euclidean_queen(from_to, row_number, col_number, rst_xres, rst_yres, d_vertical, d_diagonal, ncores);
}

template <typename T>
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_euclidean_rook(Rcpp::XPtr<T> from_to, Rcpp::IntegerVector& row_number, const double rst_xres,
  const double rst_yres, const int ncores) {
  
  const std::size_t nedges = from_to->from.size();
  std::vector<double>* distance = new std::vector<double>(nedges);
  
  #pragma omp simd
  for(std::size_t i = 0; i < nedges; ++i) {
    if(row_number[from_to->from[i]] == row_number[from_to->to[i]]) {
      (*distance)[i] = rst_xres;
    } else {
      (*distance)[i] = rst_yres;
    }
  }
  
  Rcpp::XPtr<std::vector<double> > p(distance);
  return p;
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_euclidean_rook_i(Rcpp::XPtr<From_To_I> from_to, Rcpp::IntegerVector& row_number, const double rst_xres,
  const double rst_yres, const int ncores) {
  
  return tr_fun_args_d_euclidean_rook(from_to, row_number, rst_xres, rst_yres, ncores);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_d_euclidean_rook_u(Rcpp::XPtr<From_To_U> from_to, Rcpp::IntegerVector& row_number, const double rst_xres,
  const double rst_yres, const int ncores) {
  
  return tr_fun_args_d_euclidean_rook(from_to, row_number, rst_xres, rst_yres, ncores);
}

template <typename T>
Rcpp::XPtr<std::vector<double> > tr_fun_args_coords(Rcpp::XPtr<T> from_to, Rcpp::NumericVector& col_row_number, const bool is_from) {
  
  const std::size_t nedges = from_to->from.size();
  std::vector<double>* c = new std::vector<double>(nedges);
  
  if(is_from) {
    #pragma omp simd
    for(std::size_t i = 0; i < nedges; ++i) {
      (*c)[i] = col_row_number[from_to->from[i]];
    }
  } else {
    #pragma omp simd
    for(std::size_t i = 0; i < nedges; ++i) {
      (*c)[i] = col_row_number[from_to->to[i]];
    }
  }
  
  Rcpp::XPtr<std::vector<double> > p(c);
  return p;
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_coords_i(Rcpp::XPtr<From_To_I> from_to, Rcpp::NumericVector& col_row_number, const bool is_from) {
  
  return tr_fun_args_coords(from_to, col_row_number, is_from);
}

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<double> > tr_fun_args_coords_u(Rcpp::XPtr<From_To_U> from_to, Rcpp::NumericVector& col_row_number, const bool is_from) {
  
  return tr_fun_args_coords(from_to, col_row_number, is_from);
}

template <typename C, typename V>
void tr_fun_args_v(const std::vector<C>& cell, Rcpp::DataFrame& crd, const std::string col_name, const std::size_t nedges, std::vector<V>* cpp_vec) {
  
  if constexpr (std::is_same_v<V, int>) {
    Rcpp::IntegerVector crd_col = crd[col_name];
    #pragma omp simd
    for(std::size_t i = 0; i < nedges; ++i) {
      (*cpp_vec)[i] = crd_col[cell[i]];
    }
  } else if constexpr (std::is_same_v<V, double>) {
    Rcpp::NumericVector crd_col = crd[col_name];
    #pragma omp simd
    for(std::size_t i = 0; i < nedges; ++i) {
      (*cpp_vec)[i] = crd_col[cell[i]];
    }
  } else if constexpr (std::is_same_v<V, std::string>) {
    Rcpp::CharacterVector crd_col = crd[col_name];
    #pragma omp simd
    for(std::size_t i = 0; i < nedges; ++i) {
      (*cpp_vec)[i] = crd_col[cell[i]];
    }
  } else if constexpr (std::is_same_v<V, bool>) {
    Rcpp::LogicalVector crd_col = crd[col_name];
    #pragma omp simd
    for(std::size_t i = 0; i < nedges; ++i) {
      (*cpp_vec)[i] = crd_col[cell[i]];
    }
  }
}

// Extract layer values corresponding to adjacency list cells
// [[Rcpp::export]]
Rcpp::List tr_fun_args_v_i(Rcpp::XPtr<From_To_I> from_to, Rcpp::DataFrame& crd, Rcpp::CharacterVector& v_vars, Rcpp::CharacterVector& v_vars_classes,
  const bool is_from) {
  
  const std::size_t nedges = from_to->from.size();
  const std::size_t nvars = v_vars.size();
  
  Rcpp::List v;
  
  for(std::size_t j = 0; j < nvars; ++j) {
    const std::string col_name = Rcpp::as<std::string>(v_vars[j]);
    if(v_vars_classes[j] == "integer") {
      std::vector<int>* cpp_vec = new std::vector<int>(nedges);
      tr_fun_args_v((is_from) ? from_to->from : from_to->to, crd, col_name, nedges, cpp_vec);
      Rcpp::XPtr<std::vector<int> > p(cpp_vec);
      v[col_name] = p;
    } else if(v_vars_classes[j] == "numeric") {
      std::vector<double>* cpp_vec = new std::vector<double>(nedges);
      tr_fun_args_v((is_from) ? from_to->from : from_to->to, crd, col_name, nedges, cpp_vec);
      Rcpp::XPtr<std::vector<double> > p(cpp_vec);
      v[col_name] = p;
    } else if(v_vars_classes[j] == "character") {
      std::vector<std::string>* cpp_vec = new std::vector<std::string>(nedges);
      tr_fun_args_v((is_from) ? from_to->from : from_to->to, crd, col_name, nedges, cpp_vec);
      Rcpp::XPtr<std::vector<std::string> > p(cpp_vec);
      v[col_name] = p;
    } else {
      std::vector<bool>* cpp_vec = new std::vector<bool>(nedges);
      tr_fun_args_v((is_from) ? from_to->from : from_to->to, crd, col_name, nedges, cpp_vec);
      Rcpp::XPtr<std::vector<bool> > p(cpp_vec);
      v[col_name] = p;
    }
  }
  
  return v;
}

// [[Rcpp::export]]
Rcpp::List tr_fun_args_v_u(Rcpp::XPtr<From_To_U> from_to, Rcpp::DataFrame& crd, Rcpp::CharacterVector& v_vars, Rcpp::CharacterVector& v_vars_classes,
  const bool is_from) {
  
  const std::size_t nedges = from_to->from.size();
  const std::size_t nvars = v_vars.size();
  
  Rcpp::List v;
  
  for(std::size_t j = 0; j < nvars; ++j) {
    const std::string col_name = Rcpp::as<std::string>(v_vars[j]);
    if(v_vars_classes[j] == "integer") {
      std::vector<int>* cpp_vec = new std::vector<int>(nedges);
      tr_fun_args_v((is_from) ? from_to->from : from_to->to, crd, col_name, nedges, cpp_vec);
      Rcpp::XPtr<std::vector<int> > p(cpp_vec);
      v[col_name] = p;
    } else if(v_vars_classes[j] == "numeric") {
      Rcpp::NumericVector crd_col = crd[col_name];
      std::vector<double>* cpp_vec = new std::vector<double>(nedges);
      tr_fun_args_v((is_from) ? from_to->from : from_to->to, crd, col_name, nedges, cpp_vec);
      Rcpp::XPtr<std::vector<double> > p(cpp_vec);
      v[col_name] = p;
    } else if(v_vars_classes[j] == "character") {
      Rcpp::CharacterVector crd_col = crd[col_name];
      std::vector<std::string>* cpp_vec = new std::vector<std::string>(nedges);
      tr_fun_args_v((is_from) ? from_to->from : from_to->to, crd, col_name, nedges, cpp_vec);
      Rcpp::XPtr<std::vector<std::string> > p(cpp_vec);
      v[col_name] = p;
    } else {
      Rcpp::LogicalVector crd_col = crd[col_name];
      std::vector<bool>* cpp_vec = new std::vector<bool>(nedges);
      tr_fun_args_v((is_from) ? from_to->from : from_to->to, crd, col_name, nedges, cpp_vec);
      Rcpp::XPtr<std::vector<bool> > p(cpp_vec);
      v[col_name] = p;
    }
  }
  
  return v;
}
