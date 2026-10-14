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
//A,Bは凸多角形の点を順番に列挙したもの
template <typename T>
vector<PointBase<T>> minkowski_sum(vector<PointBase<T>> A, vector<PointBase<T>> B) {
    if (A.empty() || B.empty()) return {};
    using P = PointBase<T>;
    auto normalize = [](vector<P>& V) {
        int n = V.size();
        for (int i = 0; i < n; i++) {
            auto c = cross(V[(i + 1) % n] - V[i],
                           V[(i + 2) % n] - V[(i + 1) % n]);
            if (equals(c, 0)) continue;
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
    P p0 = A[0] + B[0];
    int NA = A.size(), NB = B.size();
    A.emplace_back(A[0]), B.emplace_back(B[0]);
    vector<P> DA(NA), DB(NB);
    for(int i = 0;i < NA;i++) DA[i] = A[i + 1] - A[i];
    for(int i = 0;i < NB;i++) DB[i] = B[i + 1] - B[i];
    A.pop_back();
    B.pop_back();
    auto comp = [&](P L, P R) -> bool {
        int a = lower_or_upper(L), b = lower_or_upper(R);
        if (a != b) return a < b;
        return cross(L,R) > 0;
    };
    rotate(DA.begin(), min_element(DA.begin(),DA.end(), comp), DA.end());
    rotate(DB.begin(), min_element(DB.begin(),DB.end(), comp), DB.end());

    int a = 0, b = 0;
    vector<P> C(NA + NB + 1);
        for(int i = 0;i < NA + NB;i++) {
        if (a == NA || (b < NB && comp(DB[b], DA[a]))) {
            C[i + 1] = C[i] + DB[b++];
        } else {
            C[i + 1] = C[i] + DA[a++];
        }
    }
    assert(equals(C[NA + NB], P(0, 0)));
    C.pop_back();
    rotate(C.begin(), min_element(C.begin(),C.end()), C.end());
    P add = p0 - C[0];
    for (auto &x : C) x += add;
    C.emplace_back(C[0]);
    vector<P> ANS;
    for (P p : C) {
        while (ANS.size() >= 2) {
            P a = ANS[ANS.size() - 2];
            P b = ANS[ANS.size() - 1];
            if (!equals(cross(b - a, p - a), 0)) break;
            auto d = dot(b - a, p - b);
            if (d < 0 && !equals(d, 0)) break;
            ANS.pop_back();
        }
        ANS.emplace_back(p);
    }
    assert(equals(ANS[0], ANS.back()));
    ANS.pop_back();
    return ANS;
}