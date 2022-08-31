#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads()  1
#define omp_get_thread_num()   0
#define omp_get_max_threads()  1
#define omp_get_thread_limit() 1
#define omp_get_num_procs()    1
#endif
// [[Rcpp::plugins(openmp)]]
#include <cmath>

// Numeric distances with queen case contiguity
// [[Rcpp::export]]
Rcpp::NumericVector dists_queen_d(Rcpp::NumericVector &lat1, Rcpp::NumericVector &lon1, Rcpp::NumericVector &lat2, Rcpp::NumericVector &lon2, double yres,
  double xres, unsigned int nrow, double ymin, bool haversine, unsigned short int ncores, double radius2 = 12742020.0) {
  size_t nedges = lat1.length();                                                // Number of edges
  Rcpp::NumericVector dist(nedges);                                             // Output distance vector
  
  // Haversine distances
  if(haversine) {
    double pi180 = 0.0174532925199433;                                          // Pi / 180
    double yres2 = std::sin(yres * pi180 / 2.0);
    double yres3 = std::pow(yres2, 2);
    
    // Distances between neighboring cells
    double d_vertical = radius2 * std::atan2(yres2, std::sqrt(1.0 - yres3));    // Distance between vertical neighbors
    double d_horizontal[nrow];                                                  // Intialization of array on distances between horizontal neighbors
    double d_diagonal[nrow];                                                    // Intialization of array on distances between diagonal neighbors
    
    // Distance equation components that are constant across iterations
    ymin += 0.5 * yres;                                                         // Centroid latitude of the grid's bottom row
    double xres2 = std::sin(xres * pi180 / 2.0);
    double xres3 = std::pow(xres2, 2);
    
    // Distances between neighboring cells
    #pragma omp parallel for num_threads(ncores) if(ncores > 1)
    for(unsigned int r = 0; r < nrow; r++) {
      double yres4 = std::cos((ymin + yres * r) * pi180);
      double d = yres4 * xres2;
      d_horizontal[r] = radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2)));
      d = yres3 + yres4 * std::cos((ymin + yres * (r + 1)) * pi180) * xres3;
      d_diagonal[r] = radius2 * std::atan2(std::sqrt(d), std::sqrt(1.0 - d));
    }
  
    // Filling output distance vector using multiple cores
    if(ncores > 1) {
      size_t v = std::floor(nedges / ncores);                                   // Basic vector length per core
      size_t m = nedges % ncores;                                               // Additional elements in first vectors, in case the edges are not evenly distributed across cores
      size_t s {};                                                              // The starting element per core
      size_t e {};                                                              // The first element of the next iteration
      #pragma omp parallel for private(s, e) num_threads(ncores)
      for(unsigned short int i = 0; i < ncores; i++) {
        if(i < m) {
          s = (v + 1) * i;
          e = s + v + 1;
        } else {
          s = v * i + m;
          e = s + v;
        }
        for(size_t j = s; j < e; j++) {
          if(lon1[j] == lon2[j]) {
            dist[j] = d_vertical;
          } else if(lat1[j] == lat2[j]) {
            dist[j] = d_horizontal[(int)((lat1[j] - ymin) / yres + 0.5)];
          } else {
            dist[j] = d_diagonal[(int)((std::fmin(lat1[j], lat2[j]) - ymin) / yres + 0.5)];
          }
        }
      }
    // Filling output distance vector using single core
    } else {
      for(size_t i = 0; i < nedges; i++) {
        if(lon1[i] == lon2[i]) {
          dist[i] = d_vertical;
        } else if(lat1[i] == lat2[i]) {
          dist[i] = d_horizontal[(int)((lat1[i] - ymin) / yres + 0.5)];
        } else {
          dist[i] = d_diagonal[(int)((std::fmin(lat1[i], lat2[i]) - ymin) / yres + 0.5)];
        }
      }
    }
  // Euclidean distances
  } else {
    double d_diagonal = std::sqrt(std::pow(xres, 2) + std::pow(yres, 2));       // Distance between diagonal neighbors
    
    // Filling output distance vector using multiple cores
    if(ncores > 1) {
      size_t v = std::floor(nedges / ncores);                                   // Basic vector length per core
      size_t m = nedges % ncores;                                               // Additional elements in first vectors, in case the edges are not evenly distributed across cores
      size_t s {};                                                              // The starting element per core
      size_t e {};                                                              // The first element of the next iteration
      #pragma omp parallel for private(s, e) num_threads(ncores)
      for(unsigned short int i = 0; i < ncores; i++) {
        if(i < m) {
          s = (v + 1) * i;
          e = s + v + 1;
        } else {
          s = v * i + m;
          e = s + v;
        }
        for(size_t j = s; j < e; j++) {
          if(lon1[j] == lon2[j]) {
            dist[j] = yres;
          } else if(lat1[j] == lat2[j]) {
            dist[j] = xres;
          } else {
            dist[j] = d_diagonal;
          }
        }
      }
    // Filling output distance vector using single core
    } else {
      for(size_t i = 0; i < nedges; i++) {
        if(lon1[i] == lon2[i]) {
          dist[i] = yres;
        } else if(lat1[i] == lat2[i]) {
          dist[i] = xres;
        } else {
          dist[i] = d_diagonal;
        }
      }
    }
  }
  return dist;
}

// Numeric distances with rook case contiguity
// [[Rcpp::export]]
Rcpp::NumericVector dists_rook_d(Rcpp::NumericVector &lat1, Rcpp::NumericVector &lon1, Rcpp::NumericVector &lat2, Rcpp::NumericVector &lon2, double yres,
  double xres, unsigned int nrow, double ymin, bool haversine, unsigned short int ncores, double radius2 = 12742020.0) {
  size_t nedges = lat1.length();                                                // Number of edges
  Rcpp::NumericVector dist(nedges);                                             // Output distance vector
  
  // Haversine distances
  if(haversine) {
    double pi180 = 0.0174532925199433;                                          // Pi / 180
    double yres2 = std::sin(yres * pi180 / 2.0);
    
    // Distances between neighboring cells
    double d_vertical = radius2 * std::atan2(yres2, std::sqrt(1.0 - std::pow(yres2, 2))); // Distance between vertical neighbors
    double d_horizontal[nrow];                                                  // Intialization of array on distances between horizontal neighbors
    
    // Distance equation components that are constant across iterations
    ymin += 0.5 * yres;                                                         // Centroid latitude of the grid's bottom row
    double xres2 = std::sin(xres * pi180 / 2.0);
    
    // Distances between neighboring cells
    #pragma omp parallel for num_threads(ncores) if(ncores > 1)
    for(unsigned int r = 0; r < nrow; r++) {
      double d = std::cos((ymin + yres * r) * pi180) * xres2;
      d_horizontal[r] = radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2)));
    }
    
    // Filling output distance vector using multiple cores
    if(ncores > 1) {
      size_t v = std::floor(nedges / ncores);                                   // Basic vector length per core
      size_t m = nedges % ncores;                                               // Additional elements in first vectors, in case the edges are not evenly distributed across cores
      size_t s {};                                                              // The starting element per core
      size_t e {};                                                              // The first element of the next iteration
      #pragma omp parallel for private(s, e) num_threads(ncores)
      for(unsigned short int i = 0; i < ncores; i++) {
        if(i < m) {
          s = (v + 1) * i;
          e = s + v + 1;
        } else {
          s = v * i + m;
          e = s + v;
        }
        for(size_t j = s; j < e; j++) {
          if(lon1[j] == lon2[j]) {
            dist[j] = d_vertical;
          } else {
            dist[j] = d_horizontal[(int)((lat1[j] - ymin) / yres + 0.5)];
          }
        }
      }
    // Filling output distance vector using single core
    } else {
      for(size_t i = 0; i < nedges; i++) {
        if(lon1[i] == lon2[i]) {
          dist[i] = d_vertical;
        } else {
          dist[i] = d_horizontal[(int)((lat1[i] - ymin) / yres + 0.5)];
        }
      }
    }
  // Euclidean distances
  } else {
    // Filling output distance vector using multiple cores
    if(ncores > 1) {
      size_t v = std::floor(nedges / ncores);                                   // Basic vector length per core
      size_t m = nedges % ncores;                                               // Additional elements in first vectors, in case the edges are not evenly distributed across cores
      size_t s {};                                                              // The starting element per core
      size_t e {};                                                              // The first element of the next iteration
      #pragma omp parallel for private(s, e) num_threads(ncores)
      for(unsigned short int i = 0; i < ncores; i++) {
        if(i < m) {
          s = (v + 1) * i;
          e = s + v + 1;
        } else {
          s = v * i + m;
          e = s + v;
        }
        for(size_t j = s; j < e; j++) {
          if(lon1[j] == lon2[j]) {
            dist[j] = yres;
          } else {
            dist[j] = xres;
          }
        }
      }
    // Filling output distance vector using single core
    } else {
      for(size_t i = 0; i < nedges; i++) {
        if(lon1[i] == lon2[i]) {
          dist[i] = yres;
        } else {
          dist[i] = xres;
        }
      }
    }
  }
  return dist;
}

// Integer distances with queen case contiguity
// [[Rcpp::export]]
Rcpp::IntegerVector dists_queen_i(Rcpp::NumericVector &lat1, Rcpp::NumericVector &lon1, Rcpp::NumericVector &lat2, Rcpp::NumericVector &lon2, double yres,
  double xres, unsigned int nrow, double ymin, bool haversine, unsigned short int ncores, double radius2 = 12742020.0) {
  size_t nedges = lat1.length();                                                // Number of edges
  Rcpp::IntegerVector dist(nedges);                                             // Output distance vector
  
  // Haversine distances
  if(haversine) {
    double pi180 = 0.0174532925199433;                                          // Pi / 180
    double yres2 = std::sin(yres * pi180 / 2.0);
    double yres3 = std::pow(yres2, 2);
    
    // Distances between neighboring cells
    int d_vertical = (int)(radius2 * std::atan2(yres2, std::sqrt(1.0 - yres3)) + 0.5); // Distance between vertical neighbors
    int d_horizontal[nrow];                                                     // Intialization of array on distances between horizontal neighbors
    int d_diagonal[nrow];                                                       // Intialization of array on distances between diagonal neighbors
    
    // Distance equation components that are constant across iterations
    ymin += 0.5 * yres;                                                         // Centroid latitude of the grid's bottom row
    double xres2 = std::sin(xres * pi180 / 2.0);
    double xres3 = std::pow(xres2, 2);
    
    // Distances between neighboring cells
    #pragma omp parallel for num_threads(ncores) if(ncores > 1)
    for(unsigned int r = 0; r < nrow; r++) {
      double yres4 = std::cos((ymin + yres * r) * pi180);
      double d = yres4 * xres2;
      d_horizontal[r] = (int)(radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2))) + 0.5);
      d = yres3 + yres4 * std::cos((ymin + yres * (r + 1)) * pi180) * xres3;
      d_diagonal[r] = (int)(radius2 * std::atan2(std::sqrt(d), std::sqrt(1.0 - d)) + 0.5);
    }
    
    // Filling output distance vector using multiple cores
    if(ncores > 1) {
      size_t v = std::floor(nedges / ncores);                                   // Basic vector length per core
      size_t m = nedges % ncores;                                               // Additional elements in first vectors, in case the edges are not evenly distributed across cores
      size_t s {};                                                              // The starting element per core
      size_t e {};                                                              // The first element of the next iteration
      #pragma omp parallel for private(s, e) num_threads(ncores)
      for(unsigned short int i = 0; i < ncores; i++) {
        if(i < m) {
          s = (v + 1) * i;
          e = s + v + 1;
        } else {
          s = v * i + m;
          e = s + v;
        }
        for(size_t j = s; j < e; j++) {
          if(lon1[j] == lon2[j]) {
            dist[j] = d_vertical;
          } else if(lat1[j] == lat2[j]) {
            dist[j] = d_horizontal[(int)((lat1[j] - ymin) / yres + 0.5)];
          } else {
            dist[j] = d_diagonal[(int)((std::fmin(lat1[j], lat2[j]) - ymin) / yres + 0.5)];
          }
        }
      }
    // Filling output distance vector using single core
    } else {
      for(size_t i = 0; i < nedges; i++) {
        if(lon1[i] == lon2[i]) {
          dist[i] = d_vertical;
        } else if(lat1[i] == lat2[i]) {
          dist[i] = d_horizontal[(int)((lat1[i] - ymin) / yres + 0.5)];
        } else {
          dist[i] = d_diagonal[(int)((std::fmin(lat1[i], lat2[i]) - ymin) / yres + 0.5)];
        }
      }
    }
  // Euclidean distances
  } else {
    int d_vertical = (int)(yres + 0.5);
    int d_horizontal = (int)(xres + 0.5);
    int d_diagonal = (int)(std::sqrt(std::pow(xres, 2) + std::pow(yres, 2)) + 0.5); // Distance between diagonal neighbors
    
    // Filling output distance vector using multiple cores
    if(ncores > 1) {
      size_t v = std::floor(nedges / ncores);                                   // Basic vector length per core
      size_t m = nedges % ncores;                                               // Additional elements in first vectors, in case the edges are not evenly distributed across cores
      size_t s {};                                                              // The starting element per core
      size_t e {};                                                              // The first element of the next iteration
      #pragma omp parallel for private(s, e) num_threads(ncores)
      for(unsigned short int i = 0; i < ncores; i++) {
        if(i < m) {
          s = (v + 1) * i;
          e = s + v + 1;
        } else {
          s = v * i + m;
          e = s + v;
        }
        for(size_t j = s; j < e; j++) {
          if(lon1[j] == lon2[j]) {
            dist[j] = d_vertical;
          } else if(lat1[j] == lat2[j]) {
            dist[j] = d_horizontal;
          } else {
            dist[j] = d_diagonal;
          }
        }
      }
    // Filling output distance vector using single core
    } else {
      for(size_t i = 0; i < nedges; i++) {
        if(lon1[i] == lon2[i]) {
          dist[i] = d_vertical;
        } else if(lat1[i] == lat2[i]) {
          dist[i] = d_horizontal;
        } else {
          dist[i] = d_diagonal;
        }
      }
    }
  }
  return dist;
}

// Haversine integer distances with rook case contiguity
// [[Rcpp::export]]
Rcpp::IntegerVector dists_rook_i(Rcpp::NumericVector &lat1, Rcpp::NumericVector &lon1, Rcpp::NumericVector &lat2, Rcpp::NumericVector &lon2, double yres,
  double xres, unsigned int nrow, double ymin, bool haversine, unsigned short int ncores, double radius2 = 12742020.0) {
  size_t nedges = lat1.length();                                                // Number of edges
  Rcpp::IntegerVector dist(nedges);                                             // Output distance vector
  
  // Haversine distances
  if(haversine) {
    double pi180 = 0.0174532925199433;                                          // Pi / 180
    double yres2 = std::sin(yres * pi180 / 2.0);
    
    // Distances between neighboring cells
    int d_vertical = (int)(radius2 * std::atan2(yres2, std::sqrt(1.0 - std::pow(yres2, 2))) + 0.5); // Distance between vertical neighbors
    int d_horizontal[nrow];                                                     // Intialization of array on distances between horizontal neighbors
    
    // Distance equation components that are constant across iterations
    ymin += 0.5 * yres;                                                         // Centroid latitude of the grid's bottom row
    double xres2 = std::sin(xres * pi180 / 2.0);
    
    // Distances between neighboring cells
    #pragma omp parallel for num_threads(ncores) if(ncores > 1)
    for(unsigned int r = 0; r < nrow; r++) {
      double d = std::cos((ymin + yres * r) * pi180) * xres2;
      d_horizontal[r] = (int)(radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2))) + 0.5);
    }
    
    // Filling output distance vector using multiple cores
    if(ncores > 1) {
      size_t v = std::floor(nedges / ncores);                                   // Basic vector length per core
      size_t m = nedges % ncores;                                               // Additional elements in first vectors, in case the edges are not evenly distributed across cores
      size_t s {};                                                              // The starting element per core
      size_t e {};                                                              // The first element of the next iteration
      #pragma omp parallel for private(s, e) num_threads(ncores)
      for(unsigned short int i = 0; i < ncores; i++) {
        if(i < m) {
          s = (v + 1) * i;
          e = s + v + 1;
        } else {
          s = v * i + m;
          e = s + v;
        }
        for(size_t j = s; j < e; j++) {
          if(lon1[j] == lon2[j]) {
            dist[j] = d_vertical;
          } else {
            dist[j] = d_horizontal[(int)((lat1[j] - ymin) / yres + 0.5)];
          }
        }
      }
    // Filling output distance vector using single core
    } else {
      for(size_t i = 0; i < nedges; i++) {
        if(lon1[i] == lon2[i]) {
          dist[i] = d_vertical;
        } else {
          dist[i] = d_horizontal[(int)((lat1[i] - ymin) / yres + 0.5)];
        }
      }
    }
  // Euclidean distances
  } else {
    int d_vertical = (int)(yres + 0.5);
    int d_horizontal = (int)(xres + 0.5);
    
    // Filling output distance vector using multiple cores
    if(ncores > 1) {
      size_t v = std::floor(nedges / ncores);                                   // Basic vector length per core
      size_t m = nedges % ncores;                                               // Additional elements in first vectors, in case the edges are not evenly distributed across cores
      size_t s {};                                                              // The starting element per core
      size_t e {};                                                              // The first element of the next iteration
      #pragma omp parallel for private(s, e) num_threads(ncores)
      for(unsigned short int i = 0; i < ncores; i++) {
        if(i < m) {
          s = (v + 1) * i;
          e = s + v + 1;
        } else {
          s = v * i + m;
          e = s + v;
        }
        for(size_t j = s; j < e; j++) {
          if(lon1[j] == lon2[j]) {
            dist[j] = d_vertical;
          } else {
            dist[j] = d_horizontal;
          }
        }
      }
    // Filling output distance vector using single core
    } else {
      for(size_t i = 0; i < nedges; i++) {
        if(lon1[i] == lon2[i]) {
          dist[i] = d_vertical;
        } else {
          dist[i] = d_horizontal;
        }
      }
    }
  }
  return dist;
}
