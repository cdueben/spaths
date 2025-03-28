// [[Rcpp::plugins(cpp20)]]
// [[Rcpp::plugins(openmp)]]
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
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstddef>
#include "coordinates.h"
#include "structs.h"
#include "graph_weights.h"

// adjacency list of edge weights
// functions are overloaded with and without precomputed weights and with int and unsigned short int graph_to adjacency lists
// std::vector<std::vector<double> > graph_weights_d
// std::vector<std::vector<float> > graph_weights_f
// std::vector<std::vector<int> > graph_weights_i
// std::vector<std::vector<unsigned short int> > graph_weights_u
// std::vector<std::vector<double> > graph_weights_d
// std::vector<std::vector<float> > graph_weights_f
// std::vector<std::vector<int> > graph_weights_i
// std::vector<std::vector<unsigned short int> > graph_weights_u
// std::vector<std::vector<double> > graph_weights_d
// std::vector<std::vector<float> > graph_weights_f
// std::vector<std::vector<int> > graph_weights_i
// std::vector<std::vector<unsigned short int> > graph_weights_u

template <typename T>
std::vector<std::vector<T> > graph_w(Rcpp::List& from_to, const std::size_t n_cells, const T t, const bool int_path, const bool from_to_r) {
  
  std::vector<std::vector<T> > graph_weights(n_cells);
  
  if(from_to_r) {
    Rcpp::IntegerVector from = from_to["from"];
    Rcpp::NumericVector weights = from_to["weights"];
    
    const std::size_t n_edges = from.size();
    
    for(std::size_t i = 0; i < n_edges; ++i) {
      graph_weights[from[i]].push_back(weights[i]);
    }
  } else {
    Rcpp::XPtr<std::vector<T> > weights = from_to["weights"];
    
    if(int_path) {
      Rcpp::XPtr<From_To_I> ft = from_to["from_to"];
      
      const std::size_t n_edges = ft->from.size();
      
      for(std::size_t i = 0; i < n_edges; ++i) {
        graph_weights[ft->from[i]].push_back((*weights)[i]);
      }
    } else {
      Rcpp::XPtr<From_To_U> ft = from_to["from_to"];
      
      const std::size_t n_edges = ft->from.size();
      
      for(std::size_t i = 0; i < n_edges; ++i) {
        graph_weights[ft->from[i]].push_back((*weights)[i]);
      }
    }
  }
  from_to["weights"] = R_NilValue;
  
  return graph_weights;
}

std::vector<std::vector<double> > graph_weights_d(Rcpp::List& from_to, const std::size_t n_cells, const bool int_path, const bool from_to_r) {
  constexpr double t {};
  return graph_w(from_to, n_cells, t, int_path, from_to_r);
}

std::vector<std::vector<float> > graph_weights_f(Rcpp::List& from_to, const std::size_t n_cells, const bool int_path, const bool from_to_r) {
  constexpr float t {};
  return graph_w(from_to, n_cells, t, int_path, from_to_r);
}

std::vector<std::vector<int> > graph_weights_i(Rcpp::List& from_to, const std::size_t n_cells, const bool int_path, const bool from_to_r) {
  constexpr int t {};
  return graph_w(from_to, n_cells, t, int_path, from_to_r);
}

std::vector<std::vector<unsigned short int> > graph_weights_u(Rcpp::List& from_to, const std::size_t n_cells, const bool int_path, const bool from_to_r) {
  constexpr unsigned short int t {};
  return graph_w(from_to, n_cells, t, int_path, from_to_r);
}

template <typename G, typename D> // G: graph_to type, D: distances type
std::vector<std::vector<D> > graph_weights(const std::vector<std::vector<G> >& graph_to, Rcpp::List& coords, const bool haversine, const bool queen,
  const int ncores, const bool rm_rvector, const double radius2, const D t) {
  // ymax is ymax - 0.5 * yres
  
  const double xres = coords["xres"];
  const double yres = coords["yres"];
  const int nrow = coords["nrow"];
  const int ncol = coords["ncol"];
  const double ymax = coords["ymax"];
  const std::vector<int> cell_numbers = get_cell_numbers(coords, rm_rvector);

  const std::size_t n_cells = cell_numbers.size();
  std::vector<std::vector<D> > dist(n_cells);                                   // Output adjacency list
  
  // Haversine distances
  if(haversine) {
    // Cell coordinates
    std::vector<int> row (n_cells);
    std::vector<int> col (n_cells);
    #pragma omp parallel for simd num_threads(ncores)
    for(std::size_t i = 0; i < n_cells; ++i) {
      row[i] = cell_numbers[i] / ncol;
      col[i] = cell_numbers[i] - row[i] * ncol;
    }

    constexpr double pi180 = 0.0174532925199433;                                // Pi / 180
    const double yres2 = std::sin(yres * pi180 / 2.0);
    
    // Distances between neighboring cells
    std::vector<D> d_horizontal(nrow);                                          // Intialization of vector on distances between horizontal neighbors
    
    // Distance equation components that are constant across iterations
    const double xres2 = std::sin(xres * pi180 / 2.0);
    
    // Queen case contiguity
    if(queen) {
      const double yres3 = std::pow(yres2, 2);
      
      // Distances between neighboring cells
      D d_vertical; // Distance between vertical neighbors
      if constexpr (std::is_same_v<D, double> || std::is_same_v<D, float>) {
        d_vertical = radius2 * std::atan2(yres2, std::sqrt(1.0 - yres3));
      } else {
        d_vertical = radius2 * std::atan2(yres2, std::sqrt(1.0 - yres3)) + 0.5;
      }
      std::vector<D> d_diagonal(nrow);                                          // Intialization of vector on distances between diagonal neighbors
      
      // Distance equation components that are constant across iterations
      const double xres3 = std::pow(xres2, 2);
      
      // Distances between neighboring cells
      #pragma omp parallel for simd num_threads(ncores)
      for(int r = 0; r < nrow; r++) {
        const double yres4 = std::cos((ymax - yres * r) * pi180);
        double d = yres4 * xres2;
        if constexpr (std::is_same_v<D, double> || std::is_same_v<D, float>) {
          d_horizontal[r] = radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2)));
        } else {
          d_horizontal[r] = radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2))) + 0.5;
        }
        d = yres3 + yres4 * std::cos((ymax - yres * (r + 1)) * pi180) * xres3;
        if constexpr (std::is_same_v<D, double> || std::is_same_v<D, float>) {
          d_diagonal[r] = radius2 * std::atan2(std::sqrt(d), std::sqrt(1.0 - d));
        } else {
          d_diagonal[r] = radius2 * std::atan2(std::sqrt(d), std::sqrt(1.0 - d)) + 0.5;
        }
      }
    
      // Filling in edge weights
      #pragma omp parallel for num_threads(ncores)
      for(std::size_t i = 0; i < n_cells; i++) {
        for(const G & to_i : graph_to[i]) {
          if(col[i] == col[to_i]) {
            dist[i].push_back(d_vertical);
          } else if(row[i] == row[to_i]) {
            dist[i].push_back(d_horizontal[row[i]]);
          } else {
            dist[i].push_back(d_diagonal[std::min(row[i], row[to_i])]);
          }
        }
      }
    // Rook case contiguity
    } else {
      // Distances between neighboring cells
      D d_vertical; // Distance between vertical neighbors
      if constexpr (std::is_same_v<D, double> || std::is_same_v<D, float>) {
        d_vertical = radius2 * std::atan2(yres2, std::sqrt(1.0 - std::pow(yres2, 2)));
      } else {
        d_vertical = radius2 * std::atan2(yres2, std::sqrt(1.0 - std::pow(yres2, 2))) + 0.5;
      }
      
      #pragma omp parallel for simd num_threads(ncores)
      for(int r = 0; r < nrow; r++) {
        const double d = std::cos((ymax - yres * r) * pi180) * xres2;
        if constexpr (std::is_same_v<D, double> || std::is_same_v<D, float>) {
          d_horizontal[r] = radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2)));
        } else {
          d_horizontal[r] = radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2))) + 0.5;
        }
      }
      
      // Filling in edge weights
      #pragma omp parallel for num_threads(ncores)
      for(std::size_t i = 0; i < n_cells; i++) {
        for(const G & to_i : graph_to[i]) {
          if(col[i] == col[to_i]) {
            dist[i].push_back(d_vertical);
          } else {
            dist[i].push_back(d_horizontal[row[i]]);
          }
        }
      }
    }
  // Euclidean distances
  } else {
    // Queen case contiguity
    if(queen) {
      // Cell coordinates
      std::vector<int> row (n_cells);
      std::vector<int> col (n_cells);
      #pragma omp parallel for simd num_threads(ncores)
      for(std::size_t i = 0; i < n_cells; ++i) {
        row[i] = cell_numbers[i] / ncol;
        col[i] = cell_numbers[i] - row[i] * ncol;
      }

      D d_diagonal; // Distance between diagonal neighbors
      if constexpr (std::is_same_v<D, double> || std::is_same_v<D, float>) {
        d_diagonal = std::sqrt(std::pow(xres, 2) + std::pow(yres, 2));
      } else {
        d_diagonal = std::sqrt(std::pow(xres, 2) + std::pow(yres, 2)) + 0.5;
      }
      
      // Filling in edge weights
      #pragma omp parallel for num_threads(ncores)
      for(std::size_t i = 0; i < n_cells; i++) {
        for(const G & to_i : graph_to[i]) {
          if(col[i] == col[to_i]) {
            dist[i].push_back(yres);
          } else if(row[i] == row[to_i]) {
            dist[i].push_back(xres);
          } else {
            dist[i].push_back(d_diagonal);
          }
        }
      }
    // Rook case contiguity
    } else {
      // Cell coordinates
      std::vector<int> col (n_cells);
      #pragma omp parallel for simd num_threads(ncores)
      for(std::size_t i = 0; i < n_cells; ++i) {
        const int row = cell_numbers[i] / ncol;
        col[i] = cell_numbers[i] - row * ncol;
      }

      // Filling in edge weights
      #pragma omp parallel for num_threads(ncores)
      for(std::size_t i = 0; i < n_cells; i++) {
        for(const G & to_i : graph_to[i]) {
          if(col[i] == col[to_i]) {
            dist[i].push_back(yres);
          } else {
            dist[i].push_back(xres);
          }
        }
      }
    }
  }
  return dist;
}

std::vector<std::vector<double> > graph_weights_d(const std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2) {
  constexpr double t {};
  return graph_weights(graph_to, coords, haversine, queen, ncores, rm_rvector, radius2, t);
}

std::vector<std::vector<float> > graph_weights_f(const std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const bool haversine, const bool queen,
  const int ncores, const bool rm_rvector, const double radius2) {
  constexpr float t {};
  return graph_weights(graph_to, coords, haversine, queen, ncores, rm_rvector, radius2, t);
}

std::vector<std::vector<int> > graph_weights_i(const std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const bool haversine, const bool queen,
  const int ncores, const bool rm_rvector, const double radius2) {
  constexpr int t {};
  return graph_weights(graph_to, coords, haversine, queen, ncores, rm_rvector, radius2, t);
}

std::vector<std::vector<unsigned short int> > graph_weights_u(const std::vector<std::vector<int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2) {
  constexpr unsigned short int t {};
  return graph_weights(graph_to, coords, haversine, queen, ncores, rm_rvector, radius2, t);
}

std::vector<std::vector<double> > graph_weights_d(const std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2) {
  constexpr double t {};
  return graph_weights(graph_to, coords, haversine, queen, ncores, rm_rvector, radius2, t);
}

std::vector<std::vector<float> > graph_weights_f(const std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2) {
  constexpr float t {};
  return graph_weights(graph_to, coords, haversine, queen, ncores, rm_rvector, radius2, t);
}

std::vector<std::vector<int> > graph_weights_i(const std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords, const bool haversine,
  const bool queen, const int ncores, const bool rm_rvector, const double radius2) {
  constexpr int t {};
  return graph_weights(graph_to, coords, haversine, queen, ncores, rm_rvector, radius2, t);
}

std::vector<std::vector<unsigned short int> > graph_weights_u(const std::vector<std::vector<unsigned short int> >& graph_to, Rcpp::List& coords,
  const bool haversine, const bool queen, const int ncores, const bool rm_rvector, const double radius2) {
  constexpr unsigned short int t {};
  return graph_weights(graph_to, coords, haversine, queen, ncores, rm_rvector, radius2, t);
}
