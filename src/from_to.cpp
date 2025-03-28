// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include "structs.h"
#include "from_to.h"

// construct adjacency list for non-NA cells and shift cell numbers
template <typename T>
void from_to_ptr(T* from_to, Rcpp::IntegerVector& cell_numbers, const bool queen, const bool global, const int rst_ncol, const int n_cells,
  const bool tr_fun_specified) {
  
  const int n_not_na_cells = cell_numbers.size();
  std::vector<unsigned short int> neighbors;
  if(queen) {
    neighbors = {1, 2, 3, 4, 5, 6, 7, 8};
  } else {
    neighbors = {2, 4, 5, 7};
  }
  
  for(const unsigned short int& n : neighbors) {
    if(n == 1) {
      if(global) {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] - rst_ncol - 1;
          if((to % rst_ncol) == 0) {
            to += rst_ncol;
          }
          if(to > 0) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      } else {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] - rst_ncol - 1;
          if(to > 0 && (to % rst_ncol) != 0) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      }
    } else if(n == 2) {
      for(int from = 0; from < n_not_na_cells; ++from) {
        int to = cell_numbers[from] - rst_ncol;
        if(to > 0) {
          to = shift_to(cell_numbers, to, n_not_na_cells);
          if(to != -1) {
            from_to->from.push_back(from);
            from_to->to.push_back(to);
          }
        }
      }
    } else if(n == 3) {
      if(global) {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] - rst_ncol + 1;
          if((to % rst_ncol) == 1) {
            to -= rst_ncol;
          }
          if(to > 0) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      } else {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] - rst_ncol + 1;
          if(to > 0 && (to % rst_ncol) != 1) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      }
    } else if(n == 4) {
      if(global) {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] - 1;
          if((to % rst_ncol) == 0) {
            to += rst_ncol;
          }
          to = shift_to(cell_numbers, to, n_not_na_cells);
          if(to != -1) {
            from_to->from.push_back(from);
            from_to->to.push_back(to);
          }
        }
      } else {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] - 1;
          if((to % rst_ncol) != 0) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      }
    } else if(n == 5) {
      if(global) {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] + 1;
          if((to % rst_ncol) == 1) {
            to -= rst_ncol;
          }
          to = shift_to(cell_numbers, to, n_not_na_cells);
          if(to != -1) {
            from_to->from.push_back(from);
            from_to->to.push_back(to);
          }
        }
      } else {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] + 1;
          if((to % rst_ncol) != 1) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      }
    } else if(n == 6) {
      if(global) {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] + rst_ncol - 1;
          if((to % rst_ncol) == 0) {
            to += rst_ncol;
          }
          if(to <= n_cells) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      } else {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] + rst_ncol - 1;
          if(to <= n_cells && (to % rst_ncol) != 0) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      }
    } else if(n == 7) {
      for(int from = 0; from < n_not_na_cells; ++from) {
        int to = cell_numbers[from] + rst_ncol;
        if(to <= n_cells) {
          to = shift_to(cell_numbers, to, n_not_na_cells);
          if(to != -1) {
            from_to->from.push_back(from);
            from_to->to.push_back(to);
          }
        }
      }
    } else if(n == 8) {
      if(global) {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] + rst_ncol + 1;
          if((to % rst_ncol) == 1) {
            to -= rst_ncol;
          }
          if(to <= n_cells) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      } else {
        for(int from = 0; from < n_not_na_cells; ++from) {
          int to = cell_numbers[from] + rst_ncol + 1;
          if(to <= n_cells && (to % rst_ncol) != 1) {
            to = shift_to(cell_numbers, to, n_not_na_cells);
            if(to != -1) {
              from_to->from.push_back(from);
              from_to->to.push_back(to);
            }
          }
        }
      }
    }
  }
}

// [[Rcpp::export]]
Rcpp::XPtr<From_To_I> from_to_ptr_i(Rcpp::IntegerVector& cell_numbers, const bool queen, const bool global, const int rst_ncol, const int n_cells,
  const bool tr_fun_specified, const std::size_t max_neighbors) {
  
  From_To_I* from_to = new From_To_I;
  from_to->from.reserve(max_neighbors);
  from_to->to.reserve(max_neighbors);
  
  from_to_ptr(from_to, cell_numbers, queen, global, rst_ncol, n_cells, tr_fun_specified);
  
  Rcpp::XPtr<From_To_I> p(from_to);
  return p;
}

// [[Rcpp::export]]
Rcpp::XPtr<From_To_U> from_to_ptr_u(Rcpp::IntegerVector& cell_numbers, const bool queen, const bool global, const int rst_ncol, const int n_cells,
  const bool tr_fun_specified, const std::size_t max_neighbors) {
  
  From_To_U* from_to = new From_To_U;
  from_to->from.reserve(max_neighbors);
  from_to->to.reserve(max_neighbors);
  
  from_to_ptr(from_to, cell_numbers, queen, global, rst_ncol, n_cells, tr_fun_specified);
  
  Rcpp::XPtr<From_To_U> p(from_to);
  return p;
}
