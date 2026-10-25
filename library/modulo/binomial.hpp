#pragma once
#include <limits>

// 素数 mod を前提とする。階乗・逆階乗を参照する添字と MAX は mod 未満とする。
template <typename T>
struct Binomial {
  vector<T> f, g, h;
  Binomial(int MAX = 0) {
    assert(T::get_mod() != 0 && "Binomial<mint>()");
    assert(MAX < T::get_mod());
    f.resize(1, T{1});
    g.resize(1, T{1});
    h.resize(1, T{1});
    while (MAX >= (int)f.size()) extend();
  }
  void extend() {
    int n = f.size();
    // 倍増で mod に達すると階乗がゼロになり、逆階乗を求められない。
    int m = min(2LL * n, (long long)T::get_mod());
    assert(n < m);
    f.resize(m);
    g.resize(m);
    h.resize(m);
    for (int i = n; i < m; i++) f[i] = f[i - 1] * T(i);
    g[m - 1] = f[m - 1].inverse();
    h[m - 1] = g[m - 1] * f[m - 2];
    for (int i = m - 2; i >= n; i--) {
      g[i] = g[i + 1] * T(i + 1);
      h[i] = g[i] * f[i - 1];
    }
  }
  T fac(int i) {
    if (i < 0) return T(0);
    assert(i < T::get_mod());
    while (i >= (int)f.size()) extend();
    return f[i];
  }

  T finv(int i) {
    if (i < 0) return T(0);
    assert(i < T::get_mod());
    while (i >= (int)g.size()) extend();
    return g[i];
  }

  T inv(int i) {
    // INT_MIN は正負を反転できないため、表を作らず直接逆元を求める。
    if (i == std::numeric_limits<int>::min()) return T(i).inverse();
    if (i < 0) return -inv(-i);
    if (i >= T::get_mod()) return T(i).inverse();
    while (i >= (int)h.size()) extend();
    return h[i];
  }

  T C(int n, int r) {
    if (n < 0 || n < r || r < 0) return T(0);
    return fac(n) * finv(n - r) * finv(r);
  }

  inline T operator()(int n, int r) { return C(n, r); }

  template <typename I>
  T multinomial(const vector<I>& r) {
    static_assert(is_integral<I>::value == true);
    int n = 0;
    for (auto& x : r) {
      if (x < 0) return T(0);
      n += x;
    }
    T res = fac(n);
    for (auto& x : r) res *= finv(x);
    return res;
  }

  template <typename I>
  T operator()(const vector<I>& r) {
    return multinomial(r);
  }

  T C_naive(int n, int r) {
    if (n < 0 || n < r || r < 0) return T(0);
    T ret = T(1);
    r = min(r, n - r);
    for (int i = 1; i <= r; ++i) ret *= inv(i) * (n--);
    return ret;
  }

  T P(int n, int r) {
    if (n < 0 || n < r || r < 0) return T(0);
    return fac(n) * finv(n - r);
  }

  T H(int n, int r) {
    if (n < 0 || r < 0) return T(0);
    return r == 0 ? 1 : C(n + r - 1, r);
  }
};
