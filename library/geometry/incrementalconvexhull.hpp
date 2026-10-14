#pragma once
#include "library/geometry/pointbase.hpp"
// 下側凸包
template <typename T, bool strict = true>
struct IncrementalConvexHull_Lower {
  using P = PointBase<T>;
  set<P> S;

  IncrementalConvexHull_Lower() {}

  int size() { return S.size(); }

  template <typename ADD_V, typename RM_V, typename ADD_E, typename RM_E>
  void add(PointBase<T> p, ADD_V add_v, RM_V rm_v, ADD_E add_e, RM_E rm_e) {
    int s = side(p);
    if (strict && s >= 0) return;
    if (!strict && s > 0) return;
    // 登録済みの点では頂点・辺のコールバックも呼ばない。
    if (!strict && s == 0 && S.find(p) != S.end()) return;

    // 点追加
    add_v(p);
    S.insert(p);

    vector<P> left;
    {
      auto it = S.find(p);
      while (it != S.begin()) {
        --it;
        if (left.empty()) {
          left.emplace_back(*it);
          continue;
        }
        auto a = *it;
        auto b = left.back();
        auto det = cross(difference(b, a), difference(p, a));
        if (strict && det > 0) break;
        if (!strict && det >= 0) break;
        left.emplace_back(a);
      }
    }

    vector<P> right;
    {
      auto it = S.find(p);
      while (1) {
        ++it;
        if (it == S.end()) break;
        if (right.empty()) {
          right.emplace_back(*it);
          continue;
        }
        auto a = right.back();
        auto b = *it;
        auto det = cross(difference(a, p), difference(b, p));
        if (strict && det > 0) break;
        if (!strict && det >= 0) break;
        right.emplace_back(b);
      }
    }

    // 点削除
    if (left.size() > 1) { S.erase(next(S.find(left.back())), S.find(p)); }
    if (right.size() > 1) { S.erase(next(S.find(p)), S.find(right.back())); }
    for(int i = 0;i < (int)left.size() - 1;++i) rm_v(left[i]);
    for(int i = 0;i < (int)right.size() - 1;++i) rm_v(right[i]);

    // 辺削除
    if (!left.empty() && !right.empty()) {
      auto a = left[0], b = right[0];
      rm_e(a, b);
    }
    for(int i = 0;i < (int)left.size() - 1;++i) {
      auto a = left[i + 1], b = left[i];
      rm_e(a, b);
    }
    for(int i = 0;i < (int)right.size() - 1;++i) {
      auto a = right[i], b = right[i + 1];
      rm_e(a, b);
    }
    // 辺追加
    if (!left.empty()) { add_e(left.back(), p); }
    if (!right.empty()) { add_e(p, right.back()); }
  }

  // 中：1, 境界：0, 外：-1
  int side(PointBase<T> p) {
    auto r = S.lower_bound(p);
    if (r == S.begin()) {
      // 全部 p 以上
      if (!S.empty() && (*r) == p) return 0;
      return -1;
    }
    if (r == S.end()) {
      // p は max より大きい
      return -1;
    }
    auto l = prev(r);
    auto p1 = *l, p2 = *r;
    auto det = cross(difference(p, p1), difference(p2, p1));
    if (det == 0) return 0;
    return (det > 0 ? -1 : 1);
  }
};

template <typename T, bool strict = true>
struct Incremental_ConvexHull {
  using P = PointBase<T>;
  using C = typename P::C;
  IncrementalConvexHull_Lower<T, strict> LOWER, UPPER;
  int cnt_E;
  C det_sum;
  P area_origin;
  C area_correction = 0;
  bool is_empty;

  Incremental_ConvexHull() : cnt_E(0), det_sum(0), is_empty(1) {}

  int size() { return cnt_E; }

  bool empty() { return is_empty; }

  template <typename REAL>
  REAL area() {
    return static_cast<REAL>(area_2()) / REAL(2);
  }
  C area_2() { return det_sum + area_correction; }

  C area_cross(P a, P b, bool upper) {
    if constexpr (is_floating_point_v<T>) {
      // 初めに追加した点を原点にする。上側凸包では座標が反転している。
      P origin = upper ? -area_origin : area_origin;
      a -= origin;
      b -= origin;
      C p = a.y * b.x;
      return fma(a.x, b.y, -p) + fma(-a.y, b.x, p);
    } else {
      return cross(a, b);
    }
  }

  void add_area(C value) {
    if constexpr (is_floating_point_v<T>) {
      // 辺の追加・削除を繰り返したときの加算誤差も補償する。
      C next = det_sum + value;
      area_correction += abs(det_sum) >= abs(value)
          ? (det_sum - next) + value : (value - next) + det_sum;
      det_sum = next;
    } else {
      det_sum += value;
    }
  }

  template <typename ADD_V, typename RM_V, typename ADD_E, typename RM_E>
  void add(PointBase<T> p, ADD_V add_v, RM_V rm_v, ADD_E add_e, RM_E rm_e) {
    if (is_empty) area_origin = p;
    is_empty = 0;
    LOWER.add(
        p, add_v, rm_v,
        [&](PointBase<T> a, PointBase<T> b) {
          add_e(a, b);
          ++cnt_E;
          add_area(area_cross(a, b, false));
        },
        [&](PointBase<T> a, PointBase<T> b) {
          rm_e(a, b);
          --cnt_E;
          add_area(-area_cross(a, b, false));
        });
    UPPER.add(
        -p, [&](PointBase<T> p) { add_v(-p); }, [&](PointBase<T> p) { rm_v(-p); },
        [&](PointBase<T> a, PointBase<T> b) {
          add_e(-a, -b);
          ++cnt_E;
          add_area(area_cross(a, b, true));
        },
        [&](PointBase<T> a, PointBase<T> b) {
          rm_e(-a, -b);
          --cnt_E;
          add_area(-area_cross(a, b, true));
        });
  }
  void add(PointBase<T> p) {
    add(
        p, [](PointBase<T> p) {}, [](PointBase<T> p) {}, [](PointBase<T> s, PointBase<T> t) {},
        [](PointBase<T> s, PointBase<T> t) {});
  }

  // 中：1、境界：0、外：-1
  int side(PointBase<T> p) {
    int a = LOWER.side(p);
    int b = UPPER.side(-p);
    if (a == 0 || b == 0) return 0;
    return min(a, b);
  }
};
