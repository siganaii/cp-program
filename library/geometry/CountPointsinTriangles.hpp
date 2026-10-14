#pragma once
#include "library/misc/RNG.hpp"
#include "library/geometry/pointbase.hpp"
#include "library/data-structure/FenwickTree01.hpp"
#include "library/geometry/anglesort.hpp"
struct Count_Points_In_Triangles {
    using P = PointBase<long long>;
    const int LIM = 1000'000'000 + 10;
    vector<P> A, B;
    vector<int> new_idx; // O から見た偏角ソート順を管理
    vector<int> point;   // A[i] と一致する B[j] の数え上げ
    vector<vector<int>> seg;    // 線分 A[i]A[j] 内にある B[k] の数え上げ
    vector<vector<int>> tri;    // OA[i]A[j] 内部にある B[k] の数え上げ
    Count_Points_In_Triangles(const vector<P>& A, const vector<P>& B) : A(A), B(B) {
        for (auto&& p: A) assert(max(abs(p.x), abs(p.y)) < LIM);
        for (auto&& p: B) assert(max(abs(p.x), abs(p.y)) < LIM);
        build();
    }
    int query(int i, int j, int k) {
        i = new_idx[i], j = new_idx[j], k = new_idx[k];
        if (i > j) swap(i, j);
        if (j > k) swap(j, k);
        if (i > j) swap(i, j);
        assert(i <= j && j <= k);
        ll d = cross(A[j] - A[i],A[k] - A[i]);
        if (d == 0) return 0;
        if (d > 0) { return tri[i][j] + tri[j][k] - tri[i][k] - seg[i][k]; }
        int x = tri[i][k] - tri[i][j] - tri[j][k];
        return x - seg[i][j] - seg[j][k] - point[j];
    }
    private:
    P take_origin() {
        int N = A.size(),M = B.size();
        while (1) {
            P O = P{-LIM, RNG(-LIM, LIM)};
            bool ok = 1;
            for(int i = 0;i < N;i++) for(int j = 0;j < N;j++) {
                if (A[i] == A[j]) continue;
                if (cross(A[i] - O,A[j] - O) == 0) ok = 0;
            }
            for(int i = 0;i < N;i++) for(int j = 0;j < M;j++) {
                if (A[i] == B[j]) continue;
                if (cross(A[i] - O,B[j] - O) == 0) ok = 0;
            }
            if (ok) return O;
        }
        return P{};
    }
    void build() {
        P O = take_origin();
        for (auto&& p: A) p = p - O;
        for (auto&& p: B) p = p - O;
        int N = A.size(), M = B.size();
        vector<int> I = angle_sort(A);
        vector<P> tmp(N);
        for(int i = 0;i < N;i++) tmp[i] = A[I[i]];
        swap(tmp,A);
        new_idx.resize(N);
        for(int i = 0;i < N;i++) new_idx[I[i]] = i;
        I = angle_sort(B);
        tmp.resize(M);
        for(int i = 0;i < M;i++) tmp[i] = B[I[i]];
        swap(tmp,B);
        point.assign(N, 0);
        seg.assign(N, vector<int>(N));
        tri.assign(N, vector<int>(N));
        // point
        for(int i = 0;i < N;i++) {
            for(int j = 0;j < M;j++) {
                if (A[i] == B[j]) ++point[i];
            }
        }
        int m = 0;
        for(int j = 0;j < N;j++) {
            // OA[i]A[j], B[k]
            while (m < M && cross(A[j],B[m]) < 0) ++m;
            vector<P> C(m);
            for(int k = 0;k < m;k++) C[k] = B[k] - A[j];
            vector<int> I(m);
            for(int i = 0;i < m;i++) I[i] = i;
            sort(all(I),
                [&](auto& a, auto& b) -> bool { return cross(C[a],C[b]) > 0; });
            vector<P> ctmp(m);
            for(int i = 0;i < m;i++) ctmp[i] = C[I[i]];
            swap(ctmp,C);
            vector<int> rk(m);
            for(int k = 0;k < m;k++) rk[I[k]] = k;
            FenwickTree_01 bit(m);
            int k = m;
            for(int i = j - 1;i >= 0;i--) {
                while (k > 0 && cross(A[i],B[k-1]) > 0) { bit.add(rk[--k], 1); }
                P p = A[i] - A[j];
                int ok = 0;
                int ng = m + 1;
                while(ng - ok > 1) {
                    int mid = (ng + ok) / 2;
                    if(cross(C[mid-1],p) > 0) ok = mid;
                    else ng = mid;
                }
                int lb = ok;
                ok = 0;
                ng = m + 1;
                while(ng - ok > 1) {
                    int mid = (ng + ok) / 2;
                    if(cross(C[mid-1],p) >= 0) ok = mid;
                    else ng = mid;
                }
                int ub = ok;
                seg[i][j] += bit.prod(lb, ub), tri[i][j] += bit.prod(lb-1);
            }
        }
    }
};
