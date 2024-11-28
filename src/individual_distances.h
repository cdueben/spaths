#ifndef INDIVIDUALDISTANCES_H
#define INDIVIDUALDISTANCES_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <cmath>
#include <algorithm>

// haversine and euclidean distances called while running Dijkstra's algorithm
// all functions are defined in the header file, because they are inline

// double haversine distance
inline double haversine_dist_d(const int from_cell_number, const int to_cell_number, const int ncol, const double xres, const double yres,
  const double ymax, double radius2) {
  
  constexpr double pi180 = 0.0174532925199433;                                  // Pi / 180
  const int from_row = from_cell_number / ncol;
  const int to_row = to_cell_number / ncol;
  
  if(from_row == to_row) {                                                      // horizontal
    const double d = std::cos((ymax - yres * from_row) * pi180) * std::sin(xres * pi180 / 2.0);
    return radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2)));
  } else if((from_cell_number - from_row * ncol) == (to_cell_number - to_row * ncol)) { // vertical
    const double d = std::sin(yres * pi180 / 2.0);
    return radius2 * std::atan2(d, std::sqrt(1.0 - std::pow(d, 2)));
  } else {                                                                      // diagonal
    const double d = std::pow(std::sin(yres * pi180 / 2.0), 2) + std::cos((ymax - yres * from_row) * pi180) * std::cos((ymax - yres * to_row) * pi180) * 
      std::pow(std::sin(xres * pi180 / 2.0), 2);
    return radius2 * std::atan2(std::sqrt(d), std::sqrt(1.0 - d));
  }
}

// double euclidean distance
inline double euclidean_dist_d(const int from_cell_number, const int to_cell_number, const int ncol, const double xres, const double yres) {
  const int from_row = from_cell_number / ncol;
  const int to_row = to_cell_number / ncol;
  
  if(from_row == to_row) {                                                      // horizontal
    return xres;
  } else if((from_cell_number - from_row * ncol) == (to_cell_number - to_row * ncol)) { // vertical
    return yres;
  } else {                                                                      // diagonal
    return std::sqrt(std::pow(xres, 2) + std::pow(yres, 2));
  }
}

// float haversine distance
inline float haversine_dist_f(const int from_cell_number, const int to_cell_number, const int ncol, const double xres, const double yres, const double ymax,
  double radius2) {
  return (float)(haversine_dist_d(from_cell_number, to_cell_number, ncol, xres, yres, ymax, radius2));
}

// float euclidean distance
inline float euclidean_dist_f(const int from_cell_number, const int to_cell_number, const int ncol, const double xres, const double yres) {
  return (float)(euclidean_dist_d(from_cell_number, to_cell_number, ncol, xres, yres));
}

// integer haversine distance
inline int haversine_dist_i(const int from_cell_number, const int to_cell_number, const int ncol, const double xres, const double yres, const double ymax,
  double radius2) {
  return (int)(haversine_dist_d(from_cell_number, to_cell_number, ncol, xres, yres, ymax, radius2) + 0.5);
}

// integer euclidean distance
inline int euclidean_dist_i(const int from_cell_number, const int to_cell_number, const int ncol, const double xres, const double yres) {
  return (int)(euclidean_dist_d(from_cell_number, to_cell_number, ncol, xres, yres) + 0.5);
}

// unsigned short integer haversine distance
inline unsigned short int haversine_dist_u(const int from_cell_number, const int to_cell_number, const int ncol, const double xres, const double yres,
  const double ymax, double radius2) {
  return (unsigned short int)(haversine_dist_d(from_cell_number, to_cell_number, ncol, xres, yres, ymax, radius2) + 0.5);
}

// unsigned short integer euclidean distance
inline unsigned short int euclidean_dist_u(const int from_cell_number, const int to_cell_number, const int ncol, const double xres, const double yres) {
  return (unsigned short int)(euclidean_dist_d(from_cell_number, to_cell_number, ncol, xres, yres) + 0.5);
}

#endif
