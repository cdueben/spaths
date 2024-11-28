#ifndef COORDINATES_H
#define COORDINATES_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <cstddef>
#include <type_traits>

inline std::vector<int> get_cell_numbers(Rcpp::List& coords, const bool rm_rvector = true) {
  Rcpp::IntegerVector cell_numbers_r = coords["cell_numbers"];
  const std::vector<int> cell_numbers_c = Rcpp::as<std::vector<int> >(cell_numbers_r);
  if(rm_rvector) {
    coords["cell_numbers"] = R_NilValue;
  }
  return cell_numbers_c;
}

// coordinates of static paths
template <typename T> // T: paths data type (int, unsigned short int)
Rcpp::List coordinates(std::vector<std::vector<T> >& paths, Rcpp::List& coords, const int n_upd_rst, const bool return_dists,
  std::vector<int>& unconnected_indices) {
  // xmin is xmin + 0.5 * xres
  // ymax is ymax - 0.5 * yres
  
  const int n_paths = paths.size();
  Rcpp::List paths_coordinates(n_paths * (n_upd_rst + 1));
  
  Rcpp::IntegerVector cell_numbers = coords["cell_numbers"];
  const bool double_coords = coords["double_coords"];
  const T ncol = coords["ncol"];
  
  if(double_coords) {
    const double xmin = coords["xmin"];
    const double ymax = coords["ymax"];
    const double xres = coords["xres"];
    const double yres = coords["yres"];
    
    Rcpp::NumericVector unconnected = Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    unconnected.attr("dim") = Rcpp::Dimension(1, 2);
    
    for(int i = 0; i < n_paths; ++i) {
      const int path_length = paths[i].size();
      if(path_length == 0) {
        paths_coordinates[i] = unconnected;
        if(!return_dists) {
          unconnected_indices.push_back(i);
        }
      } else {
        Rcpp::NumericVector paths_coordinates_i (path_length * 2);
        for(int j = path_length - 1; j > -1; --j) {
          const T cell_number = cell_numbers[paths[i][j]];
          const T row = cell_number / ncol;
          paths_coordinates_i[j] = xmin + (cell_number - row * ncol) * xres;       // longitude
          paths_coordinates_i[path_length + j] = ymax - row * yres;                // latitude
        }
        paths_coordinates_i.attr("dim") = Rcpp::Dimension(path_length, 2);
        paths_coordinates[i] = paths_coordinates_i;
      }
    }
  } else {
    const int xmin = coords["xmin"];
    const int ymax = coords["ymax"];
    const int xres = coords["xres"];
    const int yres = coords["yres"];
    
    Rcpp::IntegerVector unconnected = Rcpp::IntegerVector::create(NA_INTEGER, NA_INTEGER);
    unconnected.attr("dim") = Rcpp::Dimension(1, 2);
    
    for(int i = 0; i < n_paths; ++i) {
      const int path_length = paths[i].size();
      if(path_length == 0) {
        paths_coordinates[i] = unconnected;
        if(!return_dists) {
          unconnected_indices.push_back(i);
        }
      } else {
        Rcpp::IntegerVector paths_coordinates_i (path_length * 2);
        for(int j = path_length - 1; j > -1; --j) {
          const T cell_number = cell_numbers[paths[i][j]];
          const T row = cell_number / ncol;
          paths_coordinates_i[j] = xmin + (cell_number - row * ncol) * xres;       // longitude
          paths_coordinates_i[path_length + j] = ymax - row * yres;                // latitude
        }
        paths_coordinates_i.attr("dim") = Rcpp::Dimension(path_length, 2);
        paths_coordinates[i] = paths_coordinates_i;
      }
    }
  }
  
  return paths_coordinates;
}

// coordinates of upd_rst paths
template <typename T> // T: upd_paths data type (int, unsigned short int)
void coordinates(Rcpp::List paths_coordinates, std::vector<std::vector<std::vector<T> > >& upd_paths, Rcpp::List& coords, const bool return_dists,
  std::vector<int>& unconnected_indices) {
  // xmin is xmin + 0.5 * xres
  // ymax is ymax - 0.5 * yres
  
  const int n_upd_rst = upd_paths.size();
  const int n_paths = upd_paths[0].size();
  constexpr T inf = (std::is_same_v<T, int>) ? -1 : std::numeric_limits<unsigned short int>::max();
  
  Rcpp::IntegerVector cell_numbers = coords["cell_numbers"];
  const bool double_coords = coords["double_coords"];
  const T ncol = coords["ncol"];
  
  if(double_coords) {
    const double xmin = coords["xmin"];
    const double ymax = coords["ymax"];
    const double xres = coords["xres"];
    const double yres = coords["yres"];
    
    Rcpp::NumericVector unconnected = Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    unconnected.attr("dim") = Rcpp::Dimension(1, 2);
    
    for(int u = 0; u < n_upd_rst; ++u) {
      const int starting_index = (u + 1) * n_paths;
      for(int i = 0; i < n_paths; ++i) {
        const int path_length = upd_paths[u][i].size();
        if(path_length == 0) {
          paths_coordinates[starting_index + i] = paths_coordinates[i];
        } else if(upd_paths[u][i][0] == inf) {
          paths_coordinates[starting_index + i] = unconnected;
          if(!return_dists) {
            unconnected_indices.push_back(starting_index + i);
          }
        } else {
          Rcpp::NumericVector paths_coordinates_i (path_length * 2);
          for(int j = path_length - 1; j > -1; --j) {
            const T cell_number = cell_numbers[upd_paths[u][i][j]];
            const T row = cell_number / ncol;
            paths_coordinates_i[j] = xmin + (cell_number - row * ncol) * xres;       // longitude
            paths_coordinates_i[path_length + j] = ymax - row * yres;                // latitude
          }
          paths_coordinates_i.attr("dim") = Rcpp::Dimension(path_length, 2);
          paths_coordinates[starting_index + i] = paths_coordinates_i;
        }
      }
    }
  } else {
    const int xmin = coords["xmin"];
    const int ymax = coords["ymax"];
    const int xres = coords["xres"];
    const int yres = coords["yres"];
    
    Rcpp::IntegerVector unconnected = Rcpp::IntegerVector::create(NA_INTEGER, NA_INTEGER);
    unconnected.attr("dim") = Rcpp::Dimension(1, 2);
    
    for(int u = 0; u < n_upd_rst; ++u) {
      const int starting_index = (u + 1) * n_paths;
      for(int i = 0; i < n_paths; ++i) {
        const int path_length = upd_paths[u][i].size();
        if(path_length == 0) {
          paths_coordinates[starting_index + i] = paths_coordinates[i];
        } else if(upd_paths[u][i][0] == inf) {
          paths_coordinates[starting_index + i] = unconnected;
          if(!return_dists) {
            unconnected_indices.push_back(starting_index + i);
          }
        } else {
          Rcpp::IntegerVector paths_coordinates_i (path_length * 2);
          for(int j = path_length - 1; j > -1; --j) {
            const T cell_number = cell_numbers[upd_paths[u][i][j]];
            const T row = cell_number / ncol;
            paths_coordinates_i[j] = xmin + (cell_number - row * ncol) * xres;       // longitude
            paths_coordinates_i[path_length + j] = ymax - row * yres;                // latitude
          }
          paths_coordinates_i.attr("dim") = Rcpp::Dimension(path_length, 2);
          paths_coordinates[starting_index + i] = paths_coordinates_i;
        }
      }
    }
  }
  
  coords["cell_numbers"] = R_NilValue;
}

#endif
