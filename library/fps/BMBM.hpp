#pragma once
#include "library/fps/Berlekamp-Massey.hpp"
#include "library/fps/kitamasa.hpp"
//2*k項必要
//O(k^2+klogklogn)
template <typename mint>
mint BMBM(long long n, const vector<mint> &s) {
  using fps = FormalPowerSeries<mint>;
  auto bm = BerlekampMassey<mint>(s);
  return kitamasa(n, fps{begin(bm), end(bm)}, fps{begin(s), end(s)});
}