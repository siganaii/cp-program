#pragma once
#include "library/fps/famous-series.hpp"

// f(x) = f_0 + f_1 * x + f_2 * x^2 + ... f_n x^nとなる
// このとき g(i) = sum_{0<=j<=i} f(j) の一般項を表す多項式を求める
template <typename mint>
FormalPowerSeries<mint> fualhuber(const FormalPowerSeries<mint>& f,
                                  Binomial<mint>& C) {
  using fps = FormalPowerSeries<mint>;
  int N = f.size();
  auto B = BernoulliEGF(N + 3, C);
  B[1] = -B[1];
  fps s(N), t(N);
  for (int i = 0; i < N; i++) s[i] = f[i] * C.fac(i), t[i] = B[i];
  fps u = s * t.rev();
  fps g(N + 1);
  g[0] = f[0];
  for (int i = 0; i < N; i++) g[i + 1] = u[i + N - 1] * C.finv(i + 1);
  return g;
}
