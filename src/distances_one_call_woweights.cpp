// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <unordered_set>
#include <queue>
#include <utility>
#include <functional>
#include "pair_types.h"
#include "target_distances.h"
#include "individual_distances.h"
#include "distances_one_call_woweights.h"

// distances without precomputed weights and with Dijkstra's algorithm run once on an undirected graph
// visited cells are marked in the graph_to vector
// void dists_one_call_woweights
// void dists_one_call_woweights
// void dists_one_call_woweights
// void dists_one_call_woweights
// void dists_one_call_woweights
// void dists_one_call_woweights
// void dists_one_call_woweights
// void dists_one_call_woweights

void dists_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const bool show_progress, std::vector<double>& distances, const double radius2) {

  // distances to all cells are initialized as infinite
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> vertex_distance (cell_numbers.size(), inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visted cells
    std::priority_queue<dPair, std::vector<dPair>, std::greater<dPair> > pq;
    pq.push(std::make_pair(0.0, start));
    vertex_distance[start] = 0.0;
  
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      std::unordered_set<int> targets_set(targets.begin(), targets.end());
      int targets_not_found = targets.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
          graph_to[current].clear();
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_to[current].clear();
        }
      }
    }
  }
  
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}

void dists_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const bool show_progress, std::vector<float>& distances, const double radius2) {

  // distances to all cells are initialized as infinite
  const float inf = std::numeric_limits<float>::infinity();
  std::vector<float> vertex_distance (cell_numbers.size(), inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visted cells
    std::priority_queue<fPair, std::vector<fPair>, std::greater<fPair> > pq;
    pq.push(std::make_pair(0.0, start));
    vertex_distance[start] = 0.0;
  
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      std::unordered_set<int> targets_set(targets.begin(), targets.end());
      int targets_not_found = targets.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
          graph_to[current].clear();
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_to[current].clear();
        }
      }
    }
  }
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}

void dists_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const bool show_progress, std::vector<int>& distances, const double radius2) {

  // distances to all cells are initialized as infinite
  const int inf = std::numeric_limits<int>::max();
  std::vector<int> vertex_distance (cell_numbers.size(), inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visted cells
    std::priority_queue<iPair, std::vector<iPair>, std::greater<iPair> > pq;
    pq.push(std::make_pair(0, start));
    vertex_distance[start] = 0;
  
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      std::unordered_set<int> targets_set(targets.begin(), targets.end());
      int targets_not_found = targets.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
          graph_to[current].clear();
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_to[current].clear();
        }
      }
    }
  }
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}

void dists_one_call_woweights(std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine,
  const bool show_progress, std::vector<unsigned short int>& distances, const double radius2) {

  // distances to all cells are initialized as infinite
  const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
  std::vector<unsigned short int> vertex_distance (cell_numbers.size(), inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visted cells
    std::priority_queue<uPair, std::vector<uPair>, std::greater<uPair> > pq;
    pq.push(std::make_pair(0, start));
    vertex_distance[start] = 0;
  
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      std::unordered_set<int> targets_set(targets.begin(), targets.end());
      int targets_not_found = targets.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current], cell_numbers[to_i],
                ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
          graph_to[current].clear();
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current], cell_numbers[to_i],
                ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_to[current].clear();
        }
      }
    }
  }
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}

void dists_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const bool show_progress, std::vector<double>& distances, const double radius2) {

  // distances to all cells are initialized as infinite
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> vertex_distance (cell_numbers.size(), inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visted cells
    std::priority_queue<duPair, std::vector<duPair>, std::greater<duPair> > pq;
    pq.push(std::make_pair(0.0, start));
    vertex_distance[start] = 0.0;
  
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      std::unordered_set<unsigned short int> targets_set(targets.begin(), targets.end());
      int targets_not_found = targets.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        unsigned short int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const unsigned short int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
          graph_to[current].clear();
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        unsigned short int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const unsigned short int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_to[current].clear();
        }
      }
    }
  }
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}

void dists_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const bool show_progress, std::vector<float>& distances, const double radius2) {

  // distances to all cells are initialized as infinite
  const float inf = std::numeric_limits<float>::infinity();
  std::vector<float> vertex_distance (cell_numbers.size(), inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visted cells
    std::priority_queue<fuPair, std::vector<fuPair>, std::greater<fuPair> > pq;
    pq.push(std::make_pair(0.0, start));
    vertex_distance[start] = 0.0;
  
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      std::unordered_set<unsigned short int> targets_set(targets.begin(), targets.end());
      int targets_not_found = targets.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        unsigned short int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const unsigned short int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
          graph_to[current].clear();
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        unsigned short int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const unsigned short int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_to[current].clear();
        }
      }
    }
  }
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}

void dists_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const bool show_progress, std::vector<int>& distances, const double radius2) {

  // distances to all cells are initialized as infinite
  const int inf = std::numeric_limits<int>::max();
  std::vector<int> vertex_distance (cell_numbers.size(), inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visted cells
    std::priority_queue<iuPair, std::vector<iuPair>, std::greater<iuPair> > pq;
    pq.push(std::make_pair(0, start));
    vertex_distance[start] = 0;
  
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      std::unordered_set<unsigned short int> targets_set(targets.begin(), targets.end());
      int targets_not_found = targets.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        unsigned short int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const unsigned short int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
          graph_to[current].clear();
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        unsigned short int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const unsigned short int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_to[current].clear();
        }
      }
    }
  }
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}

void dists_one_call_woweights(std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const bool show_progress, std::vector<unsigned short int>& distances,
  const double radius2) {

  // distances to all cells are initialized as infinite
  const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
  std::vector<unsigned short int> vertex_distance (cell_numbers.size(), inf);
  
  {
    // priority queue stores distances of already assessed, but not yet visted cells
    std::priority_queue<uuPair, std::vector<uuPair>, std::greater<uuPair> > pq;
    pq.push(std::make_pair(0, start));
    vertex_distance[start] = 0;
  
    // if algorithm is set to stop once all targets have been visited
    if(early_stopping) {
      std::unordered_set<unsigned short int> targets_set(targets.begin(), targets.end());
      int targets_not_found = targets.size();
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        unsigned short int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const unsigned short int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current], cell_numbers[to_i],
                ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
          graph_to[current].clear();
        }
      }
    } else {
      // loop until all queued cells have been visited
      while(!pq.empty()) {
        // select cell with lowest total distance from priority queue
        unsigned short int current = pq.top().second;
        pq.pop();
        
        // check cell, if it has not been visited yet
        if(!graph_to[current].empty()) {
          for(const unsigned short int & to_i : graph_to[current]) {
            // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
            if(!graph_to[to_i].empty()) {
              const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current], cell_numbers[to_i],
                ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
              // if the path going through the current cell is of a lower distance than the previous total distance
              if(distance_i < vertex_distance[to_i]) {
                vertex_distance[to_i] = distance_i;
                pq.push(std::make_pair(distance_i, to_i));
              }
            }
          }
          graph_to[current].clear();
        }
      }
    }
  }
  // distances to the target cells
  stat_target_distances(vertex_distance, targets, 0, -1, -1, -1, show_progress, distances);
}
