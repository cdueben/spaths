#ifndef PAIRTYPES_H
#define PAIRTYPES_H

// [[Rcpp::plugins(cpp20)]]
#include <Rcpp.h>
#include <utility>

typedef std::pair<double, int> dPair;
typedef std::pair<float, int> fPair;
typedef std::pair<int, int> iPair;
typedef std::pair<unsigned short int, int> uPair;
typedef std::pair<double, unsigned short int> duPair;
typedef std::pair<float, unsigned short int> fuPair;
typedef std::pair<int, unsigned short int> iuPair;
typedef std::pair<unsigned short int, unsigned short int> uuPair;

#endif
