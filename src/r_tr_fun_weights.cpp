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
#include <algorithm>
#include <cmath>
#include <type_traits>
#include "structs.h"

// check the output of a user-defined C++ tr_fun transition function
template <typename T>
void check_weights(const std::size_t from_to_size, Rcpp::XPtr<std::vector<T> > weights) {
  if(from_to_size != weights->size()) {
    Rcpp::stop("The vector returned by your tr_fun transition function does not have the correct length");
  }
  if constexpr (!std::is_same_v<T, unsigned short int>) {
    if(*std::min_element(weights->begin(), weights->end()) < 0) {
      Rcpp::stop("tr_fun must not return negative values");
    }
  }
  if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
    bool non_finite = false;
    const std::size_t weights_size = weights->size();
    #pragma omp parallel for simd reduction(||:non_finite)
    for(std::size_t i = 0; i < weights_size; ++i) {
      if(!std::isfinite((*weights)[i])) {
          non_finite = true;
      }
    }
    if(non_finite) {
      Rcpp::stop("tr_fun must exclusively return finite values");
    }
  }
}

// [[Rcpp::export]]
void check_weights_i_d(Rcpp::XPtr<From_To_I> from_to, Rcpp::XPtr<std::vector<double> > weights) {
  check_weights(from_to->from.size(), weights);
}

// [[Rcpp::export]]
void check_weights_i_f(Rcpp::XPtr<From_To_I> from_to, Rcpp::XPtr<std::vector<float> > weights) {
  check_weights(from_to->from.size(), weights);
}

// [[Rcpp::export]]
void check_weights_i_i(Rcpp::XPtr<From_To_I> from_to, Rcpp::XPtr<std::vector<int> > weights) {
  check_weights(from_to->from.size(), weights);
}

// [[Rcpp::export]]
void check_weights_i_u(Rcpp::XPtr<From_To_I> from_to, Rcpp::XPtr<std::vector<unsigned short int> > weights) {
  check_weights(from_to->from.size(), weights);
}

// [[Rcpp::export]]
void check_weights_u_d(Rcpp::XPtr<From_To_U> from_to, Rcpp::XPtr<std::vector<double> > weights) {
  check_weights(from_to->from.size(), weights);
}

// [[Rcpp::export]]
void check_weights_u_f(Rcpp::XPtr<From_To_U> from_to, Rcpp::XPtr<std::vector<float> > weights) {
  check_weights(from_to->from.size(), weights);
}

// [[Rcpp::export]]
void check_weights_u_i(Rcpp::XPtr<From_To_U> from_to, Rcpp::XPtr<std::vector<int> > weights) {
  check_weights(from_to->from.size(), weights);
}

// [[Rcpp::export]]
void check_weights_u_u(Rcpp::XPtr<From_To_U> from_to, Rcpp::XPtr<std::vector<unsigned short int> > weights) {
  check_weights(from_to->from.size(), weights);
}
