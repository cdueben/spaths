#ifndef DISTANCESMULTIPLECALLSWOWEIGHTS_H
#define DISTANCESMULTIPLECALLSWOWEIGHTS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <unordered_set>
#include <queue>
#include <utility>
#include <functional>
#include <cstddef>
#include <type_traits>
#include "targets_set.h"
#include "target_distances.h"
#include "individual_distances.h"
#include "visited.h"

// distances without precomputed edge weights
// functions are overloaded with int and unsigned short int graph_to adjacency lists, double, float, int, and unsigned short int distances, and int and
// unsigned short int affected_paths
template <typename G, typename D, typename A> // G: graph_to type, D: distances type, A: affected_paths type
void dists_multiple_calls_woweights(const std::vector<std::vector<G> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const G start, const std::vector<G>& targets, const bool early_stopping, const bool haversine,
  const std::unordered_set<G>& upd_rst, const bool show_progress, const std::vector<A>& affected_paths, std::vector<D>& distances, const double radius2,
  const int starting_index = -1, const int n_targets = -1, const int begin_target = -1) {
  
  const std::size_t n_cells = cell_numbers.size();
  
  // distances to all cells are initialized as infinite
  constexpr D inf = (std::is_same_v<D, double> || std::is_same_v<D, float>) ? std::numeric_limits<D>::infinity() : std::numeric_limits<D>::max();
  std::vector<D> vertex_distance (n_cells, inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visited cells
    std::priority_queue<std::pair<D, G>, std::vector<std::pair<D, G> >, std::greater<std::pair<D, G> > > pq;
    constexpr D start_distance = (std::is_same_v<D, double> || std::is_same_v<D, float>) ? 0.0 : 0;
    pq.push(std::make_pair(start_distance, start));
    vertex_distance[start] = start_distance;
  
    // indicator of visited cells
    std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      const std::unordered_set<G> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target, -1);
      int targets_not_found = targets_set.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        G current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!visited[current]) {
          for(const G & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!visited[to_i]) {
              D distance_i;
              if constexpr (std::is_same_v<D, double>) {
                distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres, ymax,
                  radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              } else if constexpr (std::is_same_v<D, float>) {
                distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres, ymax,
                  radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              } else if constexpr (std::is_same_v<D, int>) {
                distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres, ymax,
                  radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              } else {
                distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres, ymax,
                  radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              }
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          if(targets_set.contains(current)) {
            --targets_not_found;
            if(targets_not_found == 0) {
              break;
            }
          }
          visited[current] = true;
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        G current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!visited[current]) {
          for(const G & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!visited[to_i]) {
              D distance_i;
              if constexpr (std::is_same_v<D, double>) {
                distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres, ymax,
                  radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              } else if constexpr (std::is_same_v<D, float>) {
                distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres, ymax,
                  radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              } else if constexpr (std::is_same_v<D, int>) {
                distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres, ymax,
                  radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              } else {
                distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres, ymax,
                  radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              }
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          visited[current] = true;
        }
      }
    }
  }
  // distances to the target cells
  target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, show_progress, distances);
}

#endif
