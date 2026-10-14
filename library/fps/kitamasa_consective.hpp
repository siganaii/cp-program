#pragma once
#include "mod_pow.hpp"

// 添字 i が範囲内なら A[i]、範囲外なら 0 を返す
template <typename mint>
mint fps_get(const FormalPowerSeries<mint>& A, int i) {
    if (i < 0 || i >= (int)A.size()) return mint(0);
    return A[i];
}

// A(x) -> A_even(x) = A_0 + A_2 x + A_4 x^2 + ...
template <typename mint>
FormalPowerSeries<mint> fps_even(const FormalPowerSeries<mint>& A) {
    FormalPowerSeries<mint> B;
    for (int i = 0; i < (int)A.size(); i += 2) {
        B.push_back(A[i]);
    }
    B.shrink();
    return B;
}

// A(x) -> A_odd(x) = A_1 + A_3 x + A_5 x^2 + ...
template <typename mint>
FormalPowerSeries<mint> fps_odd(const FormalPowerSeries<mint>& A) {
    FormalPowerSeries<mint> B;
    for (int i = 1; i < (int)A.size(); i += 2) {
        B.push_back(A[i]);
    }
    B.shrink();
    return B;
}

// A(x) -> A(-x)
template <typename mint>
FormalPowerSeries<mint> fps_neg_x(FormalPowerSeries<mint> A) {
    for (int i = 1; i < (int)A.size(); i += 2) {
        A[i] = -A[i];
    }
    return A;
}

// r が小さい場合の愚直な展開
// P(x) / Q(x) の x^l から x^{r-1} までの係数を返す
// 0 <= l <= r を仮定する
template <typename mint>
FormalPowerSeries<mint> NaiveRangeRational(
    FormalPowerSeries<mint> P,
    FormalPowerSeries<mint> Q,
    long long l,
    long long r
) {
    vector<mint> coef(r);
    mint inv_q0 = Q[0].inverse();
    for (int i = 0; i < r; i++) {
        mint s = fps_get(P, i);
        for (int j = 1; j < (int)Q.size() && j <= i; j++) {
            s -= Q[j] * coef[i - j];
        }
        coef[i] = s * inv_q0;
    }
    FormalPowerSeries<mint> res;
    for (long long i = l; i < r; i++) {
        res.push_back(coef[i]);
    }
    return res;
}

// deg(P) < deg(Q) を仮定した区間係数取得
// P(x) / Q(x) の x^l から x^{r-1} までの係数を返す
// l >= 0, d = deg(Q), w = r-l のとき O(M(d) log(l+1) + M(d+w))
// M(n) は畳み込みの計算量
template <typename mint>
FormalPowerSeries<mint> RangeRationalProper(
    FormalPowerSeries<mint> P,
    FormalPowerSeries<mint> Q,
    long long l,
    long long r
) {
    if (l >= r) return {};

    // 負次数の係数は 0
    if (r <= 0) {
        return FormalPowerSeries<mint>(r - l, mint(0));
    }

    if (l < 0) {
        FormalPowerSeries<mint> res(-l, mint(0));
        auto right = RangeRationalProper(P, Q, 0, r);
        res.insert(res.end(), right.begin(), right.end());
        return res;
    }

    P.shrink();
    Q.shrink();

    if (P.empty()) {
        return FormalPowerSeries<mint>(r - l, mint(0));
    }

    assert(!Q.empty());
    assert(Q[0] != mint(0));

    static constexpr long long THRESHOLD = 64;

    if (r <= THRESHOLD) {
        return NaiveRangeRational(P, Q, l, r);
    }

    using fps = FormalPowerSeries<mint>;
    assert(P.size() < Q.size());
    const int d = (int)Q.size() - 1;
    const int width = (int)(r - l);

    // 開始位置が小さい場合は、必要な位置まで直接展開する
    if (l <= d) {
        auto coef = (P.pre((int)r) * Q.inv((int)r)).pre((int)r);
        return fps(coef.begin() + l, coef.end());
    }

    // h(x) = x^l mod rev(Q) とすると、漸化式から
    // a[l+i] = sum_j h[j] * a[i+j] が成り立つ。これらの和を畳み込みでまとめて求める
    auto h = mod_pow(l, fps{0, 1}, Q.rev());
    h.resize(d);
    const int terms = d + width - 1;
    auto coef = (P * Q.inv(terms)).pre(terms);
    auto product = coef * h.rev();
    return fps(product.begin() + d - 1, product.begin() + terms);
}

// P(x) / Q(x) の x^l から x^{r-1} までの係数を返す
template <typename mint>
FormalPowerSeries<mint> LinearRecurrenceRange(
    long long l,
    long long r,
    FormalPowerSeries<mint> Q,
    FormalPowerSeries<mint> P
) {
    if (l >= r) return {};

    Q.shrink();
    P.shrink();

    assert(!Q.empty());
    assert(Q[0] != mint(0));

    FormalPowerSeries<mint> ret(r - l, mint(0));

    // 分子の次数が分母以上の場合は、多項式部分を先に処理する
    if ((int)P.size() >= (int)Q.size()) {
        auto R = P / Q;
        P -= R * Q;
        P.shrink();

        for (long long i = l; i < r; i++) {
            if (0 <= i && i < (long long)R.size()) {
                ret[i - l] += R[i];
            }
        }
    }

    if (P.empty()) {
        return ret;
    }

    auto frac = RangeRationalProper(P, Q, l, r);

    for (int i = 0; i < (int)frac.size(); i++) {
        ret[i] += frac[i];
    }

    return ret;
}


