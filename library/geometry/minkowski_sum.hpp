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
    auto normalize = [&](vector<P>& V) {
        V.erase(unique(V.begin(), V.end()), V.end());
        if (V.size() > 1 && V.front() == V.back()) V.pop_back();
        int n = V.size();
        for (int i = 0; i < n; i++) {
            int c = geometry_detail::orientation_sign(V[i], V[(i + 1) % n], V[(i + 2) % n]);
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
    auto comp = [&](const vector<P>& L, int i, const vector<P>& R, int j) {
        const P& a = L[i];
        const P& b = L[(i + 1) % L.size()];
        const P& c = R[j];
        const P& d = R[(j + 1) % R.size()];
        int x = lower_or_upper(difference(b, a)), y = lower_or_upper(difference(d, c));
        if (x != y) return x < y;
        // 丸めた辺同士ではなく、元の端点から差分の誤差も含めて比較する。
        return geometry_detail::difference_cross(a, b, c, d) > 0;
    };
    int start_a = 0, start_b = 0;
    for (int i = 1; i < NA; i++) if (comp(A, i, A, start_a)) start_a = i;
    for (int i = 1; i < NB; i++) if (comp(B, i, B, start_b)) start_b = i;

    int a = 0, b = 0;
    vector<P> C;
    C.reserve(NA + NB);
    while (a < NA || b < NB) {
        int ia = (start_a + a) % NA, ib = (start_b + b) % NB;
        // 辺を累積加算せず、元の頂点の和から座標を求める。
        C.push_back(A[ia] + B[ib]);
        if (a == NA || (b < NB && comp(B, ib, A, ia))) {
            ++b;
        } else if (b < NB && !comp(A, ia, B, ib)) {
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
            int c = geometry_detail::orientation_sign(a, b, p);
            if (c > 0) break;
            // 頂点和の丸めで生じた時計回りの折れ返しも取り除く。
            if (c == 0 && dot(difference(b, a), difference(p, b)) < 0) break;
            ANS.pop_back();
        }
        ANS.emplace_back(p);
    }
    assert(equals(ANS[0], ANS.back()));
    ANS.pop_back();
    return ANS;
}
