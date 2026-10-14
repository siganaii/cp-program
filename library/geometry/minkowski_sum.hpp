#pragma once
#include "library/geometry/geometry-base.hpp"

template <typename T> int lower_or_upper(const PointBase<T> &p) {
  if (p.y != 0)
    return (p.y > 0 ? 1 : -1);
  if (p.x > 0)
    return -1;
  if (p.x < 0)
    return 1;
  return 0;
}
// A,Bは凸多角形の点を順番に列挙したもの。O(|A|+|B|)。
// 整数座標では結果の各座標がTの範囲に収まることを要求する。
template <typename T>
vector<PointBase<T>> minkowski_sum(vector<PointBase<T>> A, vector<PointBase<T>> B) {
    if (A.empty() || B.empty()) return {};
    using P = PointBase<T>;
    auto turn = [](const auto& u, const auto& v) {
        auto c = cross(u, v);
        if constexpr (is_floating_point_v<T>) {
            // 平行な辺を丸め誤差で逆順にしたり、余分な頂点を残したりしない。
            auto tol = EPS + 16 * numeric_limits<T>::epsilon()
                * (abs(u.x * v.y) + abs(u.y * v.x));
            return c > tol ? 1 : c < -tol ? -1 : 0;
        } else {
            return c > 0 ? 1 : c < 0 ? -1 : 0;
        }
    };
    auto normalize = [&](vector<P>& V) {
        V.erase(unique(V.begin(), V.end()), V.end());
        if (V.size() > 1 && V.front() == V.back()) V.pop_back();
        int n = V.size();
        for (int i = 0; i < n; i++) {
            int c = turn(difference(V[(i + 1) % n], V[i]),
                         difference(V[(i + 2) % n], V[(i + 1) % n]));
            if (c == 0) continue;
            if (c < 0) reverse(V.begin(), V.end());
            break;
        }
        rotate(V.begin(), min_element(V.begin(), V.end()), V.end());
    };
    normalize(A);
    normalize(B);
    if (A.size() > B.size()) swap(A, B);
    if (A.size() == 1) {
        for (auto &p : B) {
            p += A[0];
        }
        return B;
    }
    int NA = A.size(), NB = B.size();
    using D = decltype(difference(A[0], A[0]));
    vector<D> DA(NA), DB(NB);
    for(int i = 0;i < NA;i++) DA[i] = difference(A[(i + 1) % NA], A[i]);
    for(int i = 0;i < NB;i++) DB[i] = difference(B[(i + 1) % NB], B[i]);
    auto comp = [&](D L, D R) -> bool {
        int a = lower_or_upper(L), b = lower_or_upper(R);
        if (a != b) return a < b;
        return turn(L, R) > 0;
    };
    int start_a = min_element(DA.begin(), DA.end(), comp) - DA.begin();
    int start_b = min_element(DB.begin(), DB.end(), comp) - DB.begin();

    int a = 0, b = 0;
    vector<P> C;
    C.reserve(NA + NB);
    while (a < NA || b < NB) {
        int ia = (start_a + a) % NA, ib = (start_b + b) % NB;
        // 辺を累積加算せず、元の頂点の和から座標を求める。
        C.push_back(A[ia] + B[ib]);
        if (a == NA || (b < NB && comp(DB[ib], DA[ia]))) {
            ++b;
        } else if (b < NB && !comp(DA[ia], DB[ib])) {
            ++a;
            ++b;
        } else {
            ++a;
        }
    }
    rotate(C.begin(), min_element(C.begin(),C.end()), C.end());
    C.emplace_back(C[0]);
    vector<P> ANS;
    for (P p : C) {
        while (ANS.size() >= 2) {
            P a = ANS[ANS.size() - 2];
            P b = ANS[ANS.size() - 1];
            if (turn(difference(b, a), difference(p, a)) != 0) break;
            auto d = dot(difference(b, a), difference(p, b));
            if (d < 0 && !equals(d, 0)) break;
            ANS.pop_back();
        }
        ANS.emplace_back(p);
    }
    assert(equals(ANS[0], ANS.back()));
    ANS.pop_back();
    return ANS;
}
