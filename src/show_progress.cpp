// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <string>
#include <iostream>
#include "show_progress.h"

// progress print

void stat_show_progress_header(const int n_paths, const bool upd_rst_defined, const bool paths, const bool bar) {
  const std::string upd = (upd_rst_defined) ? "static " : "";
  const std::string data_type = (paths) ? "paths" : "distances";
  Rcpp::Rcout << "Starting " << upd << data_type << " calculation" << std::endl;
  if(bar) {
    Rcpp::Rcout << '|' << std::string(n_paths, '-') << '|' << std::endl << '|';
  }
}

void upd_show_progress_header(const int n_upd_rst, const bool paths, const bool bar) {
  const std::string data_type = (paths) ? "paths" : "distances";
  Rcpp::Rcout << "Starting updated " << data_type << " calculation" << std::endl;
  if(bar) {
    Rcpp::Rcout << '|' << std::string(n_upd_rst, '-') << '|' << std::endl << '|';
  }
}
