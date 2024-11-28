#ifndef PATHSONECALLWOWEIGHTS_H
#define PATHSONECALLWOWEIGHTS_H

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
#include "individual_distances.h"
#include "target_distances.h"
#include "target_paths.h"

// paths without precomputed weights and with Dijkstra's algorithm run once on an undirected graph
// visited cells are marked in the graph_to vector
template <typename G, typename D> // G: graph_to type, D: distances type
void paths_one_call_woweights(std::vector<std::vector<G> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const G start, const std::vector<G>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<G>& graph_to_0, const bool show_progress, std::vector<std::vector<G> >& paths, std::vector<D>& distances, const double radius2) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<G> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    constexpr D inf = (std::is_same_v<D, double> || std::is_same_v<D, float>) ? std::numeric_limits<D>::infinity() : std::numeric_limits<D>::max();
    std::vector<D> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<std::pair<D, G>, std::vector<std::pair<D, G> >, std::greater<std::pair<D, G> > > pq;
      constexpr D start_distance = (std::is_same_v<D, double> || std::is_same_v<D, float>) ? 0.0 : 0;
      pq.push(std::make_pair(start_distance, start));
      vertex_distance[start] = start_distance;
      
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        std::unordered_set<G> targets_set(targets.begin(), targets.end());
        int targets_not_found = targets.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          G current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!graph_to[current].empty()) {
            for(const G & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!graph_to[to_i].empty()) {
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
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            if(targets_set.contains(current)) {
              --targets_not_found;
              if(targets_not_found == 0) {
                all_visited = true;
                break;
              }
            }
            graph_to[current].clear();
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          G current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!graph_to[current].empty()) {
            for(const G & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!graph_to[to_i].empty()) {
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
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            graph_to[current].clear();
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, false, distances);
    }
  }
  const std::vector<int> affected_paths;
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, 0, -1, -1, -1, show_progress, paths);
}

#endif
