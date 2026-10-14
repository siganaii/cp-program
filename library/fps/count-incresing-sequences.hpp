#pragma once
#include "library/modulo/binomial.hpp"
#include "library/ntt/ntt.hpp"
NTT<mint> ntt;
template<class mint>
pair<vector<mint>, vector<mint>> count_square(vector<mint> L, vector<mint> D,Binomial<mint> &C){
    assert(!L.empty() && !D.empty());
    int N = L.size();
    int M = D.size();
    if (min(N, M) <= 400){
        int sw = 0;
        if (N > M) swap(N, M), swap(L, D), sw = 1;
        vector<mint> R(N);
        for (int i = 0; i < N; i++){
            D[0] += L[i];
            for (int j = 1; j < M; j++) D[j] += D[j - 1];
            R[i] = D.back();
        }
        if (sw) swap(R, D);
        return {R, D};
    }
    auto cyclic_convolution = [&](vector<mint> f,vector<mint> g) {
        ntt.ntt(f);
        ntt.ntt(g);
        for(int i = 0;i < (int)f.size();i++) f[i] *= g[i];
        ntt.intt(f);
        return f;
    };
    vector<mint> R(N), U(M);
    int z = 0;
    while ((1 << z) < (N + M - 1)) z++;
    // 左から右
    {
        vector<mint> tmp(N);
        for (int i = 0; i < N; i++) tmp[i] = C.C(M - 1 + i, i);
        tmp = ntt.multiply(tmp, L);
        for (int i = 0; i < N; i++) R[i] += tmp[i];
    }
    // 左から上
    {
        vector<mint> tmp(1 << z);
        for (int i = 0; i < N; i++) L[i] *= C.finv(N - 1 - i);
        for (int i = 0; i < N + M - 1; i++) tmp[i] = C.fac(i);
        L.resize(1 << z,0);
        tmp = cyclic_convolution(tmp,L);
        for (int i = 0; i < M; i++) U[i] += tmp[N - 1 + i] * C.finv(i);
    }
    // 下から上
    {
        vector<mint> tmp(M);
        for (int i = 0; i < M; i++) tmp[i] = C.C(N - 1 + i, i);
        tmp = ntt.multiply(tmp, D);
        for (int i = 0; i < M; i++) U[i] += tmp[i];
    }
    // 下から右
    {
        vector<mint> tmp(1 << z);
        for (int i = 0; i < M; i++) D[i] *= C.finv(M - i - 1);
        for (int i = 0; i < N + M - 1; i++) tmp[i] = C.fac(i);
        D.resize(1 << z,0);
        tmp = cyclic_convolution(tmp, D);
        for (int i = 0; i < N; i++) R[i] += tmp[M - 1 + i] * C.finv(i);
    }
    return {R, U};
}
template<class mint>
/*
 * g(A, x) を
 * 0 <= B[i] < A[i] かつ B[i] = x を満たす
 * 広義単調増加列 B の数とする
 * res[x] = sum C[i] * g(A[i:N], x)
 * を返す
 */
vector<mint> count_increase_sequences_with_upper_bounds(Binomial<mint> &BC, vector<int> A, vector<mint> C){
    int N = A.size();
    assert((int)C.size() == N);
    assert(N);
    for (int i = (int)(A.size()) - 1; i > 0; i--) A[i - 1] = min(A[i - 1], A[i]);
    if (A.back() == 0) return {};
    if (min(A.back(), N) <= 400){
        vector<mint> dp(0);
        dp.reserve(A.back());
        for (int i = 0; i < N; i++){
            dp.resize(A[i], 0);
            if (A[i]) dp[0] += C[i];
            for (int j = 1; j < (int)dp.size(); j++){
                dp[j] += dp[j - 1];
            }
        }
        return dp;
    }
    if (N == 1){
        vector<mint> res(A[0]);
        for (int i = 0; i < A[0]; i++) res[i] = C[0];
        return res;
    }
    int m = N / 2;
    vector<int> LA(m), RA(N - m);
    vector<mint> LC(m), RC(N - m);
    for (int i = 0; i < m; i++){
        LA[i] = A[i];
        LC[i] = C[i];
    }
    for (int i = 0; i < N - m; i++){
        RA[i] = A[i + m] - A[m - 1];
        RC[i] = C[i + m];
    }
    vector<mint> res;
    res.reserve(A.back());
    auto L = count_increase_sequences_with_upper_bounds(BC,LA,LC);
    if (!L.empty()){
        auto [R, U] = count_square(L,RC,BC);
        for (int i = 0; i < (int)R.size(); i++) res.push_back(R[i]);
        swap(U, RC);
    }
    auto R = count_increase_sequences_with_upper_bounds(BC,RA,RC);
    for (auto x : R) res.push_back(x);
    return res;
}

template<class mint>
/*
 * f(a, b) を X[0] = a, X[N - 1] = b であるような、A, B に挟まれたものとする
 * 長さ B[N - 1] - A[N - 1] を返す
 * res[b - A.back()] = sum C[a - A[0]] * f(a, b)
 * A, B は広義単調増加が嬉しい
 * C は空ならば、全て 1 であるとする。
 * そうでないなら、|C| = B[0] - A[0] でないといけない
 */
vector<mint> count_increase_sequences_with_upper_lower_bounds(Binomial<mint> &BC,vector<int> A, vector<int> B,vector<mint> C = {}){
    int N = A.size();
    assert(A.size() == B.size());
    assert(N > 0);
    assert(C.empty() || (int)C.size() == B[0] - A[0]);
    for (int i = 0; i < N - 1; i++){
        A[i + 1] = max(A[i], A[i + 1]);
    }
    for (int i = N - 1; i > 0; i--){
        B[i - 1] = min(B[i], B[i - 1]);
    }
    if (A.back() >= B.back()) return {};
    // A[0] == 0 にする
    vector<mint> res(B.back() - A.back(), 0);
    {
        int tmp = A[0];
        for (int i = 0; i < N; i++){
            A[i] -= tmp;
            B[i] -= tmp;
            if (A[i] >= B[i]) return res;
        }
    }
    if (C.empty()){
        C.resize(B[0] - A[0], 1);
    }
    // Tightening the upper bounds can remove impossible starting values.
    else C.resize(B[0] - A[0]);
    int l = 0;
    while (B[l] <= A.back()){
        for (int i = (int)(C.size()) - 1; i > 0; i--) C[i] -= C[i - 1];
        int nl = l;
        while (A[nl] < B[l]) nl++;
        vector<int> tmp(B[l] - A[l]);
        tmp[0] = 1;
        for (int i = l; i < nl; i++){
            tmp[A[i] - A[l]]++;
        }
        for (int i = 1; i < B[l] - A[l]; i++) tmp[i] += tmp[i - 1];
        auto X = count_increase_sequences_with_upper_bounds(BC,tmp,C);
        vector<int> nB(nl - l + 1);
        for (int i = l; i <= nl; i++){
            nB[i - l] = B[i] - B[l];
        }
        auto Y = count_increase_sequences_with_upper_bounds(BC,nB,X);
        C.resize(B[nl] - A[nl]);
        for (int i = 0; i < B[nl] - A[nl]; i++){
            C[i] = Y[i + A[nl] - B[l]];
        }
        l = nl;
    }
    // A を揃えてしまえ
    {
        int a = A[l];
        for (int i = l; i < N; i++){
            A[i] -= a;
            B[i] -= a;
        }
    }
    for (int i = (int)(C.size()) - 1; i > 0; i--) C[i] -= C[i - 1];
    vector<mint> D(N - l, 0);
    if (A.back() != 0){
        vector<mint> L(A.back());
        for (int i = 0; i < (int)L.size(); i++) L[i] = C[i];
        vector<int> tmp(L.size());
        tmp[0] = 1;
        for (int i = l; i < N; i++){
            if (A[i] < (int)tmp.size()) tmp[A[i]]++;
        }
        for (int i = 1; i < (int)tmp.size(); i++){
            tmp[i] += tmp[i - 1];
        }
        auto nD = count_increase_sequences_with_upper_bounds(BC,tmp,L);
        for (int i = 0; i < (int)nD.size(); i++) D[i] = nD[i];
    }
    for (int i = A.back(); i < B[l]; i++) C[i - A.back()] = C[i];
    C.resize(B[l] - A.back());
    auto [R, U] = count_square(C,D,BC);
    res = R;
    vector<int> nB(N - l);
    for (int i = 0; i < N - l; i++) nB[i] = B[i + l] - B[l];
    R = count_increase_sequences_with_upper_bounds(BC,nB,U);
    for (auto x : R) res.push_back(x);
    return res;
}
