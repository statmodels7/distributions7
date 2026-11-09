// R evaluates these expressions one operation at a time; clang (arm64
// macOS) would otherwise fuse a*b + c into one rounding
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
#include <Rcpp.h>
using namespace Rcpp;
