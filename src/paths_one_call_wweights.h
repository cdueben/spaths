#ifndef PATHSONECALLWWEIGHTS_H
#define PATHSONECALLWWEIGHTS_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <utility>
#include <functional>
#include <cstddef>
#include <type_traits>
#include "target_distances.h"
#include "target_paths.h"

// paths with precomputed edge weights and Dijkstra's algorithm run once on an undirected graph
// visited cells are marked in the graph_weights vector
template <typename G, typename D> // G: graph_to type, D: distances type
void paths_one_call_wweights(const std::vector<std::vector<G> >& graph_to, std::vector<std::vector<D> >& graph_weights, const std::size_t n_cells,
  const G start, const std::vector<G>& targets, const bool early_stopping, const int ncores, const std::unordered_set<G>& graph_to_0,
  const bool show_progress, std::vector<std::vector<G> >& paths, std::vector<D>& distances) {
  
  // ids of the predecessor cells on the shortest paths
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
        const std::unordered_set<G> targets_set(targets.begin(), targets.end());
        int targets_not_found = targets.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          G current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!graph_weights[current].empty()) {
            const std::size_t neighbors = graph_to[current].size();
            for(std::size_t i = 0; i < neighbors; ++i) {
              const G to_i = graph_to[current][i];
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!graph_weights[to_i].empty()) {
                const D distance_i = vertex_distance[current] + graph_weights[current][i];
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
            graph_weights[current].clear();
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          G current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!graph_weights[current].empty()) {
            const std::size_t neighbors = graph_to[current].size();
            for(std::size_t i = 0; i < neighbors; ++i) {
              const G to_i = graph_to[current][i];
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!graph_weights[to_i].empty()) {
                const D distance_i = vertex_distance[current] + graph_weights[current][i];
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            graph_weights[current].clear();
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
