#pragma once
#include "library/modulo/binomial.hpp"
template <typename mint>
mint lagrange_interpolation(const vector<mint>& y, long long x,Binomial<mint>& C) {
    int N = (int)y.size() - 1;
    if (0 <= x && x <= N) return y[x];
    mint ret = 0;
    vector<mint> dp(N + 1, 1), pd(N + 1, 1);
    mint a = x;
    for (int i = 0; i < N; i++) dp[i + 1] = dp[i] * a, --a;
    for (int i = N; i > 0; i--) pd[i - 1] = pd[i] * a, ++a;
    for (int i = 0; i <= N; i++) {
        mint tmp = y[i] * dp[i] * pd[i] * C.finv(i) * C.finv(N - i);
        ret += ((N - i) & 1) ? -tmp : tmp;
    }
    return ret;
}
