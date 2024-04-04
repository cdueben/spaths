// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <limits>
#include <unordered_set>
#include <queue>
#include <utility>
#include <functional>
#include <cstddef>
#include "pair_types.h"
#include "targets_set.h"
#include "individual_distances.h"
#include "target_distances.h"
#include "target_paths.h"
#include "visited.h"
#include "paths_multiple_calls_woweights.h"

// paths without precomputed weights
// functions are overloaded with int and unsigned int paths, double, float, int, and unsigned int distances, and int and unsigned short int affected paths
// void paths_multiple_calls_woweights
// void paths_multiple_calls_woweights
// void paths_multiple_calls_woweights
// void paths_multiple_calls_woweights
// void paths_multiple_calls_woweights
// void paths_multiple_calls_woweights
// void paths_multiple_calls_woweights
// void paths_multiple_calls_woweights

void paths_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<int> >& paths, std::vector<double>& distances, const double radius2, const int starting_index, const int n_targets,
  const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<dPair, std::vector<dPair>, std::greater<dPair> > pq;
      pq.push(std::make_pair(0.0, start));
      vertex_distance[start] = 0.0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol,
                  xres, yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol,
                  xres, yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<int> >& paths, std::vector<float>& distances, const double radius2, const int starting_index,
  const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const float inf = std::numeric_limits<float>::infinity();
    std::vector<float> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<fPair, std::vector<fPair>, std::greater<fPair> > pq;
      pq.push(std::make_pair(0.0, start));
      vertex_distance[start] = 0.0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target,
          -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<int> >& paths, std::vector<int>& distances, const double radius2, const int starting_index,
  const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const int inf = std::numeric_limits<int>::max();
    std::vector<int> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<iPair, std::vector<iPair>, std::greater<iPair> > pq;
      pq.push(std::make_pair(0, start));
      vertex_distance[start] = 0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target,
          -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<int> >& paths, std::vector<unsigned short int>& distances, const double radius2, const int starting_index,
  const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
    std::vector<unsigned short int> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<uPair, std::vector<uPair>, std::greater<uPair> > pq;
      pq.push(std::make_pair(0, start));
      vertex_distance[start] = 0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target,
          -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current],
                  cell_numbers[to_i], ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current],
                  cell_numbers[to_i], ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<double>& distances, const double radius2, const int starting_index,
  const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<unsigned short int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<duPair, std::vector<duPair>, std::greater<duPair> > pq;
      pq.push(std::make_pair(0.0, start));
      vertex_distance[start] = 0.0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<unsigned short int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets,
          begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol,
                  xres, yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol,
                  xres, yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<float>& distances, const double radius2, const int starting_index,
  const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<unsigned short int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const float inf = std::numeric_limits<float>::infinity();
    std::vector<float> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<fuPair, std::vector<fuPair>, std::greater<fuPair> > pq;
      pq.push(std::make_pair(0.0, start));
      vertex_distance[start] = 0.0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<unsigned short int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets,
          begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<int>& distances, const double radius2, const int starting_index,
  const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<unsigned short int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const int inf = std::numeric_limits<int>::max();
    std::vector<int> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<iuPair, std::vector<iuPair>, std::greater<iuPair> > pq;
      pq.push(std::make_pair(0, start));
      vertex_distance[start] = 0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<unsigned short int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets,
          begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<unsigned short int>& distances, const double radius2,
  const int starting_index, const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<unsigned short int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
    std::vector<unsigned short int> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<uuPair, std::vector<uuPair>, std::greater<uuPair> > pq;
      pq.push(std::make_pair(0, start));
      vertex_distance[start] = 0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<unsigned short int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets,
          begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current],
                  cell_numbers[to_i], ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current],
                  cell_numbers[to_i], ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths, std::vector<double>& distances, const double radius2,
  const int starting_index, const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<dPair, std::vector<dPair>, std::greater<dPair> > pq;
      pq.push(std::make_pair(0.0, start));
      vertex_distance[start] = 0.0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol,
                  xres, yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol,
                  xres, yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths, std::vector<float>& distances, const double radius2,
  const int starting_index, const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const float inf = std::numeric_limits<float>::infinity();
    std::vector<float> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<fPair, std::vector<fPair>, std::greater<fPair> > pq;
      pq.push(std::make_pair(0.0, start));
      vertex_distance[start] = 0.0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target,
          -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths, std::vector<int>& distances, const double radius2,
  const int starting_index, const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const int inf = std::numeric_limits<int>::max();
    std::vector<int> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<iPair, std::vector<iPair>, std::greater<iPair> > pq;
      pq.push(std::make_pair(0, start));
      vertex_distance[start] = 0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target,
          -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol, const double xres,
  const double yres, const double ymax, const int start, const std::vector<int>& targets, const bool early_stopping, const bool haversine, const int ncores,
  const std::unordered_set<int>& graph_to_0, const std::unordered_set<int>& upd_rst, const bool show_progress,
  const std::vector<unsigned short int>& affected_paths, std::vector<std::vector<int> >& paths, std::vector<unsigned short int>& distances,
  const double radius2, const int starting_index, const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
    std::vector<unsigned short int> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<uPair, std::vector<uPair>, std::greater<uPair> > pq;
      pq.push(std::make_pair(0, start));
      vertex_distance[start] = 0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets, begin_target,
          -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current],
                  cell_numbers[to_i], ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current],
                  cell_numbers[to_i], ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<double>& distances, const double radius2, const int starting_index, const int n_targets,
  const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<unsigned short int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<duPair, std::vector<duPair>, std::greater<duPair> > pq;
      pq.push(std::make_pair(0.0, start));
      vertex_distance[start] = 0.0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<unsigned short int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets,
          begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol,
                  xres, yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const double distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_d(cell_numbers[current], cell_numbers[to_i], ncol,
                  xres, yres, ymax, radius2) : euclidean_dist_d(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<float>& distances, const double radius2, const int starting_index, const int n_targets,
  const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<unsigned short int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const float inf = std::numeric_limits<float>::infinity();
    std::vector<float> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<fuPair, std::vector<fuPair>, std::greater<fuPair> > pq;
      pq.push(std::make_pair(0.0, start));
      vertex_distance[start] = 0.0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<unsigned short int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets,
          begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const float distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_f(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<int>& distances, const double radius2, const int starting_index, const int n_targets,
  const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<unsigned short int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const int inf = std::numeric_limits<int>::max();
    std::vector<int> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<iuPair, std::vector<iuPair>, std::greater<iuPair> > pq;
      pq.push(std::make_pair(0, start));
      vertex_distance[start] = 0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<unsigned short int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets,
          begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres,
                  yres, ymax, radius2) : euclidean_dist_i(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}

void paths_multiple_calls_woweights(const std::vector<std::vector<unsigned short int> >& graph_to, const std::vector<int>& cell_numbers, const int ncol,
  const double xres, const double yres, const double ymax, const unsigned short int start, const std::vector<unsigned short int>& targets,
  const bool early_stopping, const bool haversine, const int ncores, const std::unordered_set<unsigned short int>& graph_to_0,
  const std::unordered_set<unsigned short int>& upd_rst, const bool show_progress, const std::vector<unsigned short int>& affected_paths,
  std::vector<std::vector<unsigned short int> >& paths, std::vector<unsigned short int>& distances, const double radius2, const int starting_index,
  const int n_targets, const int begin_target) {
  
  // ids of the predecessor cells on the shortest paths
  const std::size_t n_cells = cell_numbers.size();
  std::vector<unsigned short int> predecessor (n_cells);
  
  // indicator of whether all targets have been visited
  bool all_visited {false};
  
  {
    // distances to all cells are initialized as infinite
    const unsigned short int inf = std::numeric_limits<unsigned short int>::max();
    std::vector<unsigned short int> vertex_distance (n_cells, inf);
    
    {
      // priority queue stores distances of already assessed, but not yet visted cells
      std::priority_queue<uuPair, std::vector<uuPair>, std::greater<uuPair> > pq;
      pq.push(std::make_pair(0, start));
      vertex_distance[start] = 0;
      
      // indicator of visited cells
      std::vector<bool> visited = create_visited(n_cells, upd_rst);
    
      // if algorithm is set to stop once all targets have been visited
      if(early_stopping) {
        const std::unordered_set<unsigned short int> targets_set = create_targets_set(targets, affected_paths.empty(), starting_index, n_targets,
          begin_target, -1);
        int targets_not_found = targets_set.size();
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current],
                  cell_numbers[to_i], ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
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
            visited[current] = true;
          }
        }
      } else {
        // loop until all queued cells have been visited
        while(!pq.empty()) {
          // select cell with lowest total distance from priority queue
          unsigned short int current = pq.top().second;
          pq.pop();
          
          // check cell, if it has not been visited yet
          if(!visited[current]) {
            for(const unsigned short int & to_i : graph_to[current]) {
              // check neighboring cell, if it has not been visited yet; the edge from the current cell can by definition not provide a lower distance
              if(!visited[to_i]) {
                const unsigned short int distance_i = vertex_distance[current] + ((haversine) ? haversine_dist_u(cell_numbers[current],
                  cell_numbers[to_i], ncol, xres, yres, ymax, radius2) : euclidean_dist_u(cell_numbers[current], cell_numbers[to_i], ncol, xres, yres));
                // if the path going through the current cell is of a lower distance than the previous total distance
                if(distance_i < vertex_distance[to_i]) {
                  vertex_distance[to_i] = distance_i;
                  predecessor[to_i] = current;
                  pq.push(std::make_pair(distance_i, to_i));
                }
              }
            }
            visited[current] = true;
          }
        }
      }
    }
    // distances of the target cells
    if(!distances.empty()) {
      target_distances(vertex_distance, affected_paths, targets, starting_index, n_targets, begin_target, -1, false, distances);
    }
  }
  
  // paths to the target cells
  target_paths(predecessor, start, targets, graph_to_0, affected_paths, all_visited, ncores, starting_index, n_targets, begin_target, -1, show_progress,
    paths);
}
