// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <cstddef>
#include "coordinates.h"

// path coordinates
// inline std::vector<int> get_cell_numbers
// Rcpp::List coordinates
// Rcpp::List coordinates
// void coordinates
// void coordinates

// coordinates of static paths
Rcpp::List coordinates(const std::vector<std::vector<int> >& paths, Rcpp::List& coords, const int n_upd_rst, const bool return_dists,
  std::vector<int>& unconnected_indices) {
  // xmin is xmin + 0.5 * xres
  // ymax is ymax - 0.5 * yres
  
  const int n_paths = paths.size();
  Rcpp::List paths_coordinates(n_paths * (n_upd_rst + 1));
  
  Rcpp::IntegerVector cell_numbers = coords["cell_numbers"];
  const bool double_coords = coords["double_coords"];
  const int ncol = coords["ncol"];
  
  if(double_coords) {
    const double xmin = coords["xmin"];
    const double ymax = coords["ymax"];
    const double xres = coords["xres"];
    const double yres = coords["yres"];
    
    Rcpp::NumericVector unconnected = Rcpp::NumericVector::create(NA_REAL, NA_REAL);
    unconnected.attr("dim") = Rcpp::Dimension(1, 2);
    
    // this loop is not multi-threaded as Rcpp objects are not thread-safe
    // std::vectors with parallelization turn out to be slower in this case because of the subsequent conversion to an R vector
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
          const int cell_number = cell_numbers[paths[i][j]];
          const int row = cell_number / ncol;
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
          const int cell_number = cell_numbers[paths[i][j]];
          const int row = cell_number / ncol;
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

Rcpp::List coordinates(std::vector<std::vector<unsigned short int> >& paths, Rcpp::List& coords, const int n_upd_rst, const bool return_dists,
  std::vector<int>& unconnected_indices) {
  // xmin is xmin + 0.5 * xres
  // ymax is ymax - 0.5 * yres
  
  const int n_paths = paths.size();
  Rcpp::List paths_coordinates(n_paths * (n_upd_rst + 1));
  
  Rcpp::IntegerVector cell_numbers = coords["cell_numbers"];
  const bool double_coords = coords["double_coords"];
  const unsigned short int ncol = coords["ncol"];
  
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
          const unsigned short int cell_number = cell_numbers[paths[i][j]];
          const unsigned short int row = cell_number / ncol;
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
          const unsigned short int cell_number = cell_numbers[paths[i][j]];
          const unsigned short int row = cell_number / ncol;
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
void coordinates(Rcpp::List paths_coordinates, std::vector<std::vector<std::vector<int> > >& upd_paths, Rcpp::List& coords, const bool return_dists,
  std::vector<int>& unconnected_indices) {
  // xmin is xmin + 0.5 * xres
  // ymax is ymax - 0.5 * yres
  
  const int n_upd_rst = upd_paths.size();
  const int n_paths = upd_paths[0].size();
  
  Rcpp::IntegerVector cell_numbers = coords["cell_numbers"];
  const bool double_coords = coords["double_coords"];
  const int ncol = coords["ncol"];
  
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
        } else if(upd_paths[u][i][0] == -1) {
          paths_coordinates[starting_index + i] = unconnected;
          if(!return_dists) {
            unconnected_indices.push_back(starting_index + i);
          }
        } else {
          Rcpp::NumericVector paths_coordinates_i (path_length * 2);
          for(int j = path_length - 1; j > -1; --j) {
            const int cell_number = cell_numbers[upd_paths[u][i][j]];
            const int row = cell_number / ncol;
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
        } else if(upd_paths[u][i][0] == -1) {
          paths_coordinates[starting_index + i] = unconnected;
          if(!return_dists) {
            unconnected_indices.push_back(starting_index + i);
          }
        } else {
          Rcpp::IntegerVector paths_coordinates_i (path_length * 2);
          for(int j = path_length - 1; j > -1; --j) {
            const int cell_number = cell_numbers[upd_paths[u][i][j]];
            const int row = cell_number / ncol;
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

void coordinates(Rcpp::List paths_coordinates, std::vector<std::vector<std::vector<unsigned short int> > >& upd_paths, Rcpp::List& coords,
  const bool return_dists, std::vector<int>& unconnected_indices) {
  // xmin is xmin + 0.5 * xres
  // ymax is ymax - 0.5 * yres
  
  const int n_upd_rst = upd_paths.size();
  const int n_paths = upd_paths[0].size();
  const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
  
  Rcpp::IntegerVector cell_numbers = coords["cell_numbers"];
  const bool double_coords = coords["double_coords"];
  const unsigned short int ncol = coords["ncol"];
  
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
            const unsigned short int cell_number = cell_numbers[upd_paths[u][i][j]];
            const unsigned short int row = cell_number / ncol;
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
            const unsigned short int cell_number = cell_numbers[upd_paths[u][i][j]];
            const unsigned short int row = cell_number / ncol;
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
