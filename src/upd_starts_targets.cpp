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
#include <vector>
#include <cstddef>
#include "upd_starts_targets.h"

// pairwise and no targets and not directed cases
// starts indices search uses binary search
void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<int>& affected_paths, const std::vector<int>& starts,
  const std::vector<int>& targets, std::vector<int>& upd_starts, std::vector<int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    bool found {false};
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        upd_starts[p] = starts[i];
        found = true;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    if(!found) {
      upd_starts[p] = starts[right];
    }
    upd_targets[p] = targets[p];
  }
}

void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<unsigned short int>& affected_paths,
  const std::vector<int>& starts, const std::vector<int>& targets, std::vector<int>& upd_starts, std::vector<int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    bool found {false};
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        upd_starts[p] = starts[i];
        found = true;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    if(!found) {
      upd_starts[p] = starts[right];
    }
    upd_targets[p] = targets[p];
  }
}

void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<int>& affected_paths,
  const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, std::vector<unsigned short int>& upd_starts,
  std::vector<unsigned short int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    bool found {false};
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        upd_starts[p] = starts[i];
        found = true;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    if(!found) {
      upd_starts[p] = starts[right];
    }
    upd_targets[p] = targets[p];
  }
}

void upd_starts_targets_pairwise(const std::vector<int>& starting_indices, const std::vector<unsigned short int>& affected_paths,
  const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets, std::vector<unsigned short int>& upd_starts,
  std::vector<unsigned short int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    bool found {false};
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        upd_starts[p] = starts[i];
        found = true;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    if(!found) {
      upd_starts[p] = starts[right];
    }
    upd_targets[p] = targets[p];
  }
}

// no targets and directed
void upd_starts_targets_no_targets_directed(const std::vector<int>& affected_paths, const std::vector<int>& starts, std::vector<int>& upd_starts,
  std::vector<int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starts_size_1 = starts.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    const int start = affected_paths[p] / starts_size_1;
    int target = affected_paths[p] % starts_size_1;
    if(start <= target) {
      ++target;
    }
    upd_starts[p] = starts[start];
    upd_targets[p] = starts[target];
    
  }
}

void upd_starts_targets_no_targets_directed(const std::vector<unsigned short int>& affected_paths, const std::vector<int>& starts,
  std::vector<int>& upd_starts, std::vector<int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const unsigned short int starts_size_1 = starts.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    const unsigned short int start = affected_paths[p] / starts_size_1;
    unsigned short int target = affected_paths[p] % starts_size_1;
    if(start <= target) {
      ++target;
    }
    upd_starts[p] = starts[start];
    upd_targets[p] = starts[target];
  }
}

void upd_starts_targets_no_targets_directed(const std::vector<int>& affected_paths, const std::vector<unsigned short int>& starts,
  std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starts_size_1 = starts.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    const int start = affected_paths[p] / starts_size_1;
    int target = affected_paths[p] % starts_size_1;
    if(start <= target) {
      ++target;
    }
    upd_starts[p] = starts[start];
    upd_targets[p] = starts[target];
  }
}

void upd_starts_targets_no_targets_directed(const std::vector<unsigned short int>& affected_paths, const std::vector<unsigned short int>& starts,
  std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const unsigned short int starts_size_1 = starts.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    const unsigned short int start = affected_paths[p] / starts_size_1;
    unsigned short int target = affected_paths[p] % starts_size_1;
    if(start <= target) {
      ++target;
    }
    upd_starts[p] = starts[start];
    upd_targets[p] = starts[target];
  }
}

// no targets and not directed
void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<int>& affected_paths,
  const std::vector<int>& starts, std::vector<int>& upd_starts, std::vector<int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        right = i;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    upd_starts[p] = starts[right];
    upd_targets[p] = starts[affected_path - starting_indices[right] + right + 1];
  }
}

void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<unsigned short int>& affected_paths,
  const std::vector<int>& starts, std::vector<int>& upd_starts, std::vector<int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        right = i;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    upd_starts[p] = starts[right];
    upd_targets[p] = starts[affected_path - starting_indices[right] + right + 1];
  }
}

void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<int>& affected_paths,
  const std::vector<unsigned short int>& starts, std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        right = i;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    upd_starts[p] = starts[right];
    upd_targets[p] = starts[affected_path - starting_indices[right] + right + 1];
  }
}

void upd_starts_targets_no_targets_not_directed(const std::vector<int>& starting_indices, const std::vector<unsigned short int>& affected_paths,
  const std::vector<unsigned short int>& starts, std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int starting_indices_size_1 = starting_indices.size() - 1;
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    int left {0};
    int right = starting_indices_size_1;
    const int affected_path = affected_paths[p];
    while(left <= right) {
      int i = left + (right - left) / 2;
      if(starting_indices[i] == affected_path) {
        right = i;
        break;
      }
      if(starting_indices[i] < affected_path) {
        left = i + 1;
      } else {
        right = i - 1;
      }
    }
    upd_starts[p] = starts[right];
    upd_targets[p] = starts[affected_path - starting_indices[right] + right + 1];
  }
}

// not pairwise
void upd_starts_targets_not_pairwise(const std::vector<int>& affected_paths, const std::vector<int>& starts, const std::vector<int>& targets,
  std::vector<int>& upd_starts, std::vector<int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int targets_size = targets.size();
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    upd_starts[p] = starts[affected_paths[p] / targets_size];
    upd_targets[p] = targets[affected_paths[p] % targets_size];
  }
}

void upd_starts_targets_not_pairwise(const std::vector<unsigned short int>& affected_paths, const std::vector<int>& starts, const std::vector<int>& targets,
  std::vector<int>& upd_starts, std::vector<int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const unsigned short int targets_size = targets.size();
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    upd_starts[p] = starts[affected_paths[p] / targets_size];
    upd_targets[p] = targets[affected_paths[p] % targets_size];
  }
}

void upd_starts_targets_not_pairwise(const std::vector<int>& affected_paths, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const int targets_size = targets.size();
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    upd_starts[p] = starts[affected_paths[p] / targets_size];
    upd_targets[p] = targets[affected_paths[p] % targets_size];
  }
}

void upd_starts_targets_not_pairwise(const std::vector<unsigned short int>& affected_paths, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, std::vector<unsigned short int>& upd_starts, std::vector<unsigned short int>& upd_targets) {
  
  const std::size_t n_affected_paths = affected_paths.size();
  const unsigned short int targets_size = targets.size();
  
  #pragma omp simd
  for(std::size_t p = 0; p < n_affected_paths; ++p) {
    upd_starts[p] = starts[affected_paths[p] / targets_size];
    upd_targets[p] = targets[affected_paths[p] % targets_size];
  }
}
