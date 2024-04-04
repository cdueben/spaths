#ifndef UPDSHOWPROGRESS_H
#define UPDSHOWPROGRESS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <string>
#include <iostream>

void stat_show_progress_header(const int n_paths, const bool upd_rst_defined, const bool paths, const bool bar);
void upd_show_progress_header(const int n_upd_rst, const bool paths, const bool bar);

#endif
