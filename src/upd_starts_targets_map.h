#ifndef UPDSTARTSTARGETSMAP_H
#define UPDSTARTSTARGETSMAP_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstddef>
#include "upd_starts_targets.h"

void upd_st_map(const std::vector<int>& affected_paths, const std::vector<int>& starts, const std::vector<int>& targets, const bool pairwise,
  const bool directed, const std::vector<int>& starting_indices, std::unordered_map<int, std::vector<int> >& u_starts_targets,
  std::unordered_map<int, std::vector<int> >& u_affected_paths);
void upd_st_map(const std::vector<unsigned short int>& affected_paths, const std::vector<int>& starts, const std::vector<int>& targets, const bool pairwise,
  const bool directed, const std::vector<int>& starting_indices, std::unordered_map<int, std::vector<int> >& u_starts_targets,
  std::unordered_map<int, std::vector<unsigned short int> >& u_affected_paths);
void upd_st_map(const std::vector<int>& affected_paths, const std::vector<unsigned short int>& starts, const std::vector<unsigned short int>& targets,
  const bool pairwise, const bool directed, const std::vector<int>& starting_indices,
  std::unordered_map<unsigned short int, std::vector<unsigned short int> >& u_starts_targets,
  std::unordered_map<unsigned short int, std::vector<int> >& u_affected_paths);
void upd_st_map(const std::vector<unsigned short int>& affected_paths, const std::vector<unsigned short int>& starts,
  const std::vector<unsigned short int>& targets, const bool pairwise, const bool directed, const std::vector<int>& starting_indices, 
  std::unordered_map<unsigned short int, std::vector<unsigned short int> >& u_starts_targets,
  std::unordered_map<unsigned short int, std::vector<unsigned short int> >& u_affected_paths);

#endif
