#pragma once
#include <immintrin.h>
struct bit_vector {
  using u32 = uint32_t;
  using u64 = uint64_t;

  static constexpr u32 w = 64;
  vector<u64> block;
  vector<u32> count;
  u32 n, zeros;

  inline u32 get(u32 i) const { return u32(block[i / w] >> (i % w)) & 1u; }
  inline void set(u32 i) { block[i / w] |= 1ULL << (i % w); }

  bit_vector() {}
  bit_vector(int _n) { init(_n); }

  __attribute__((optimize("O3", "unroll-loops"))) void init(int _n) {
    n = zeros = _n;
    block.resize(n / w + 1, 0);
    count.resize(block.size(), 0);
  }

  __attribute__((target("popcnt"))) void build() {
    for (u32 i = 1; i < block.size(); ++i)
      count[i] = count[i - 1] + _mm_popcnt_u64(block[i - 1]);
    zeros = rank0(n);
  }

  inline u32 rank0(u32 i) const { return i - rank1(i); }

  __attribute__((target("bmi2,popcnt"))) inline u32 rank1(u32 i) const {
    return count[i / w] + _mm_popcnt_u64(_bzhi_u64(block[i / w], i % w));
  }
};

// 負の数が出てこないようにする
template <typename T, bool EnableSum = false>
struct WaveletMatrix {
  using u32 = uint32_t;
  using u64 = uint64_t;

  int n, lg;
  T max_value = 0;
  vector<T> a;
  vector<bit_vector> bv;
  vector<vector<T>> pref;
  vector<T> orig_pref;

  WaveletMatrix(u32 _n) : n(_n), a(n) {}
  WaveletMatrix(const vector<T>& _a) : n(_a.size()), a(_a) { build(); }

  __attribute__((optimize("O3"))) void build() {
    max_value = a.empty() ? T(0) : *max_element(begin(a), end(a));
    lg = a.empty() ? 0 : __lg(max<T>(max_value, 1)) + 1;

    if constexpr (EnableSum) {
      orig_pref.assign(n + 1, 0);
      for (int i = 0; i < n; ++i) orig_pref[i + 1] = orig_pref[i] + a[i];
      pref.assign(lg, vector<T>(n + 1, 0));
    }

    bv.assign(lg, n);
    vector<T> cur = a, nxt(n);

    for (int h = lg - 1; h >= 0; --h) {
      for (int i = 0; i < n; ++i)
        if ((cur[i] >> h) & 1) bv[h].set(i);

      bv[h].build();

      array<decltype(begin(nxt)), 2> it{begin(nxt), begin(nxt) + bv[h].zeros};
      for (int i = 0; i < n; ++i) *it[bv[h].get(i)]++ = cur[i];

      if constexpr (EnableSum) {
        for (int i = 0; i < n; ++i) pref[h][i + 1] = pref[h][i] + nxt[i];
      }

      swap(cur, nxt);
    }
  }

  void set(u32 i, const T& x) {
    assert(x >= 0);
    a[i] = x;
  }

  inline pair<u32, u32> succ0(int l, int r, int h) const {
    return {bv[h].rank0(l), bv[h].rank0(r)};
  }

  inline pair<u32, u32> succ1(int l, int r, int h) const {
    u32 l0 = bv[h].rank0(l), r0 = bv[h].rank0(r), zeros = bv[h].zeros;
    return {l + zeros - l0, r + zeros - r0};
  }

  T access(u32 k) const {
    T ret = 0;
    for (int h = lg - 1; h >= 0; --h) {
      u32 f = bv[h].get(k);
      ret |= f ? T(1) << h : 0;
      k = f ? bv[h].rank1(k) + bv[h].zeros : bv[h].rank0(k);
    }
    return ret;
  }

  T kth_smallest(u32 l, u32 r, u32 k) const {
    T res = 0;
    for (int h = lg - 1; h >= 0; --h) {
      u32 l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
      if (k < r0 - l0) {
        l = l0, r = r0;
      } else {
        k -= r0 - l0;
        res |= T(1) << h;
        l += bv[h].zeros - l0;
        r += bv[h].zeros - r0;
      }
    }
    return res;
  }

  T kth_largest(int l, int r, int k) const {
    return kth_smallest(l, r, r - l - k - 1);
  }

  int range_freq(int l, int r, T upper) const {
    if (l >= r || upper <= 0) return 0;
    if (upper > max_value) return r - l;

    int ret = 0;
    for (int h = lg - 1; h >= 0; --h) {
      bool f = (upper >> h) & 1;
      u32 l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
      if (f) {
        ret += r0 - l0;
        l += bv[h].zeros - l0;
        r += bv[h].zeros - r0;
      } else {
        l = l0, r = r0;
      }
    }
    return ret;
  }

  int range_freq(int l, int r, T lower, T upper) const {
    return range_freq(l, r, upper) - range_freq(l, r, lower);
  }

  T range_sum(int l, int r, T upper) const {
    static_assert(EnableSum, "range_sum requires WaveletMatrix<T, true>");
    if (l >= r || upper <= 0) return 0;
    if (upper > max_value) return orig_pref[r] - orig_pref[l];

    T ret = 0;
    for (int h = lg - 1; h >= 0; --h) {
      bool f = (upper >> h) & 1;
      u32 l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
      if (f) {
        ret += pref[h][r0] - pref[h][l0];
        l += bv[h].zeros - l0;
        r += bv[h].zeros - r0;
      } else {
        l = l0, r = r0;
      }
    }
    return ret;
  }

  T range_sum(int l, int r, T lower, T upper) const {
    static_assert(EnableSum, "range_sum requires WaveletMatrix<T, true>");
    return range_sum(l, r, upper) - range_sum(l, r, lower);
  }

  T topk_sum(int l, int r, int k) const {
    static_assert(EnableSum, "topk_sum requires WaveletMatrix<T, true>");
    if (k <= 0 || l >= r) return T(0);
    if (k >= r - l) return orig_pref[r] - orig_pref[l];

    T ret = 0, value = 0;
    for (int h = lg - 1; h >= 0; --h) {
      u32 l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
      u32 one_cnt = (r - l) - (r0 - l0);
      u32 l1 = l + bv[h].zeros - l0;
      u32 r1 = r + bv[h].zeros - r0;

      if ((int)one_cnt <= k) {
        ret += pref[h][r1] - pref[h][l1];
        k -= one_cnt;
        l = l0, r = r0;
      } else {
        value |= T(1) << h;
        l = l1, r = r1;
      }
    }
    return ret + value * T(k);
  }

  T prev_value(int l, int r, T upper) const {
    int cnt = range_freq(l, r, upper);
    return cnt == 0 ? T(-1) : kth_smallest(l, r, cnt - 1);
  }

  T next_value(int l, int r, T lower) const {
    int cnt = range_freq(l, r, lower);
    return cnt == r - l ? T(-1) : kth_smallest(l, r, cnt);
  }
};
