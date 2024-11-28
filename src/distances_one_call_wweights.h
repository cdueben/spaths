#ifndef DISTANCESONECALLWWEIGHTS_H
#define DISTANCESONECALLWWEIGHTS_H

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
#include "target_distances.h"

// distances with precomputed edge weights and Dijkstra's algorithm run once on an undirected graph
// visited cells are marked in the graph_weights vector
// functions are overloaded with int and unsigned short int graph_to adjacency lists and double, float, int, and unsigned short int distances
template <typename G, typename D> // G: graph_to type, D: distances type
void dists_one_call_wweights(const std::vector<std::vector<G> >& graph_to, std::vector<std::vector<D> >& graph_weights, const std::size_t n_cells,
  const G start, const std::vector<G>& targets, const bool early_stopping, const bool show_progress, std::vector<D>& distances) {
  
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
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_weights[current].clear();
        }
      }
    }
  }
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}

#endif
